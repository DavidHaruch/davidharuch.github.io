import sys
import time
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets
from labjack import ljm
from scipy import signal
import csv

class ThermalController(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        # --- Controller Configuration ---
        self.setpoint =28.0
        self.t_amb_init = 22.10
        self.k_thermal = 6.13
        self.k_thermal = 19.49
        
        # Normal Operating Gains (Scheduled)
        self.high_gains = {'kp': 0.3, 'ki': 0.005, 'kd': 0.0}
        self.low_gains = {'kp': 0.2, 'ki': 0.005, 'kd': 0.0}
        
        self.high_gains = {'kp': 0.4, 'ki': 0.015, 'kd': 0.0}
        self.low_gains = {'kp': 0.4, 'ki': 0.015, 'kd': 0.0}
        
        self.high_gains = {'kp': 0.8, 'ki': 0.015, 'kd': 0.0}
        self.low_gains = {'kp': 0.8, 'ki': 0.015, 'kd': 0.0}
        
        
        # Identification "Lazy" Gains
        self.id_gains = {'kp': 0.05, 'ki': 0.01, 'kd': 0.0}
        
        # NTC Constants
        self.V_ref, self.R_fixed, self.B_ntc = 1.8, 10000.0, 3950.0
        self.T0, self.R0 = 298.15, 10000.0

        # Bode Config
        self.test_duration = 600   
        self.test_duration = 2400   
        self.f_min, self.f_max = 0.0005, 1.0
        self.num_frequencies, self.injection_amp, self.fs = 65, 0.3, 40.0             
        
        # --- State Variables ---
        self.is_identifying = False
        self.integral, self.last_error = 0.0, 0.0
        self.start_time = time.time()
        self.start_bode_time = 0
        
        # Buffers for Live Plotting
        self.max_points = 800 
        self.time_data = np.zeros(self.max_points)
        self.y_data = np.full(self.max_points, self.t_amb_init)
        self.amb_data = np.full(self.max_points, self.t_amb_init)
        self.u_data = np.zeros(self.max_points)
        
        self.u_history, self.y_history, self.d_history = [], [], []
        
        self.handle = None
        self.init_hardware()
        self.generate_multisine()
        self.init_ui()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.control_loop)
        if self.handle:
            self.timer.start(25) # 40Hz

    def init_hardware(self):
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            
            # --- Mixed Range Differential Setup ---
            # AIN0/2 (Pos) = 1V, AIN1/3 (Neg/Ref) = 10mV
            configs = [
                "AIN0_NEGATIVE_CH", 1, "AIN0_RANGE", 1.0, "AIN0_RESOLUTION_INDEX", 8,
                "AIN0_SETTLING_US", 5000,
                "AIN1_RANGE", 0.01, "AIN1_RESOLUTION_INDEX", 8,
                
                "AIN2_NEGATIVE_CH", 3, "AIN2_RANGE", 1.0, "AIN2_RESOLUTION_INDEX", 8,
                "AIN2_SETTLING_US", 5000,
                "AIN3_RANGE", 0.01, "AIN3_RESOLUTION_INDEX", 8,
                
                "DAC0", self.V_ref
            ]
            ljm.eWriteNames(self.handle, len(configs)//2, configs[::2], configs[1::2])
            
            # PWM Setup (DIO0)
            self.roll_value = int((80000000 / 256) / 5.0)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_DIVISOR", 256)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ROLL_VALUE", self.roll_value)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 1)
            ljm.eWriteName(self.handle, "DIO0_EF_INDEX", 0)
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 1)
            print("Hardware Initialized: Swapped Feedback to AIN2-3, Ambient to AIN0-1.")
        except Exception as e:
            print(f"HW Init Error: {e}")
            self.handle = None

    def voltage_to_temp(self, v_out):
        v_out = max(0.001, min(v_out, self.V_ref - 0.001))
        r_ntc = (v_out * self.R_fixed) / (self.V_ref - v_out)
        inv_t = (1.0 / self.T0) + (1.0 / self.B_ntc) * np.log(r_ntc / self.R0)
        return (1.0 / inv_t) - 273.15

    def control_loop(self):
        if not self.handle: return
        try:
            # READINGS SWAPPED: 
            # res[0] is AIN0 (Ambient), res[1] is AIN2 (Feedback/Control)
            res = ljm.eReadNames(self.handle, 2, ["AIN0", "AIN2"])
            y_amb = self.voltage_to_temp(res[0])
            y_sec = self.voltage_to_temp(res[1])
            
            now = time.time()
            dt, error = 0.025, self.setpoint - y_sec
            abs_err = abs(error)

            # Gain Selection
            if self.is_identifying:
                kp, ki, kd = self.id_gains['kp'], self.id_gains['ki'], self.id_gains['kd']
            else:
                gains = self.high_gains if abs_err < 0.2 else self.low_gains
                kp, ki, kd = gains['kp'], gains['ki'], gains['kd']
            
            # Integrator logic
            if abs_err < 4.0: 
                self.integral = np.clip(self.integral + error * dt, -100, 100)
            else:
                self.integral = 0.0
            
            # Control Calculation
            u_base = (kp * error) + (ki * self.integral) + ((self.setpoint - y_amb)/self.k_thermal)

            if self.is_identifying:
                elapsed = now - self.start_bode_time
                sig = np.sum(np.cos(2 * np.pi * self.target_freqs * elapsed + self.phases))
                inj = (sig / np.sqrt(len(self.target_freqs))) * self.injection_amp
                final_u = np.clip(u_base + inj, 0, 1)
                
                self.u_history.append(final_u); self.y_history.append(y_sec); self.d_history.append(y_amb)
                self.progress_bar.setValue(int((elapsed / self.test_duration) * 100))
                if elapsed >= self.test_duration: self.finalize_scan()
            else:
                final_u = np.clip(u_base, 0, 1)

            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", int(final_u * self.roll_value))

            # UI Update
            self.time_data = np.roll(self.time_data, -1); self.time_data[-1] = now - self.start_time
            self.y_data = np.roll(self.y_data, -1); self.y_data[-1] = y_sec
            self.amb_data = np.roll(self.amb_data, -1); self.amb_data[-1] = y_amb
            self.u_data = np.roll(self.u_data, -1); self.u_data[-1] = final_u * 100
            
            self.curve_sec.setData(self.time_data, self.y_data)
            self.curve_amb.setData(self.time_data, self.amb_data)
            self.curve_u.setData(self.time_data, self.u_data)

        except Exception as e:
            print(f"Loop Error: {e}")
            self.cleanup_hardware()

    def generate_multisine(self):
        freqs = np.logspace(np.log10(self.f_min), np.log10(self.f_max), self.num_frequencies)
        self.target_freqs = np.unique(np.round(freqs * self.test_duration) / self.test_duration)
        N = len(self.target_freqs)
        self.phases = np.array([np.pi * k**2 / N for k in range(N)])

    def toggle_bode(self):
        if not self.is_identifying:
            self.u_history, self.y_history, self.d_history = [], [], []
            self.start_bode_time = time.time()
            self.is_identifying = True
            self.btn_bode.setText("Cancel Scan (ID Gains Active)")
        else:
            self.is_identifying = False
            self.btn_bode.setText("Start Bode Scan")

    def init_ui(self):
        self.setWindowTitle("Thermal Controller - Differential (Mixed Ranges)")
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        layout = QtWidgets.QVBoxLayout(self.central_widget)
        self.win = pg.GraphicsLayoutWidget()
        layout.addWidget(self.win)
        
        self.p1 = self.win.addPlot(title="Feedback Temp (AIN2-3)"); self.curve_sec = self.p1.plot(pen='c')
        self.p1.addItem(pg.InfiniteLine(pos=self.setpoint, angle=0, pen=pg.mkPen('g', style=QtCore.Qt.DashLine)))
        self.win.nextRow()
        self.p_amb = self.win.addPlot(title="Ambient Temp (AIN0-1)"); self.curve_amb = self.p_amb.plot(pen='m')
        self.p_amb.setXLink(self.p1)
        self.win.nextRow()
        self.p_u = self.win.addPlot(title="Heater Duty %"); self.curve_u = self.p_u.plot(pen='y')
        self.p_u.setYRange(0, 100); self.p_u.setXLink(self.p1)

        btn_box = QtWidgets.QHBoxLayout()
        self.btn_bode = QtWidgets.QPushButton("Start Bode Scan")
        self.btn_bode.clicked.connect(self.toggle_bode)
        self.progress_bar = QtWidgets.QProgressBar()
        btn_box.addWidget(self.btn_bode); btn_box.addWidget(self.progress_bar)
        layout.addLayout(btn_box)

    def finalize_scan(self):
        self.is_identifying = False
        self.btn_bode.setText("Start Bode Scan")
        u, y, d = np.array(self.u_history), np.array(self.y_history), np.array(self.d_history)
        ts = int(time.time())
        filename = f"bode_data_{ts}.csv"
        with open(filename, 'w', newline='') as f:
            csv.writer(f).writerow(['u', 'y', 'd'])
            csv.writer(f).writerows(zip(u, y, d))
        print(f"Scan Saved: {filename}")

    def cleanup_hardware(self):
        if self.handle:
            print("\nSAFETY SHUTDOWN: Setting heater to 0%...")
            try:
                self.timer.stop()
                ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", 0)
                ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 0)
                ljm.close(self.handle)
                self.handle = None
                print("Hardware Released.")
            except:
                pass

    def closeEvent(self, event):
        self.cleanup_hardware()
        event.accept()

def exception_handler(exctype, value, traceback):
    # This ensures that even if Python crashes (SyntaxError, etc), we try to stop the heater
    if hasattr(ctrl, 'cleanup_hardware'):
        ctrl.cleanup_hardware()
    sys.__excepthook__(exctype, value, traceback)

if __name__ == '__main__':
    sys.excepthook = exception_handler # Catch internal crashes
    app = QtWidgets.QApplication(sys.argv)
    ctrl = ThermalController()
    try:
        ctrl.show()
        sys.exit(app.exec_())
    finally:
        ctrl.cleanup_hardware()