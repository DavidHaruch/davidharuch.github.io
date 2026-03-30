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
        
        # --- Configs ---
        self.setpoint = 26.0
        self.t_amb = 23.0
        self.k_thermal = 22.0
        self.threshold = 0.2
        self.hysteresis = 0.02
        self.high_gains = {'kp': 0.4, 'ki': 0.01, 'kd': 0.0}
        self.low_gains = {'kp': 0.15, 'ki': 0.005, 'kd': 0.0}
        
        self.current_mode = "LOW"
        self.kp, self.ki, self.kd = self.low_gains['kp'], self.low_gains['ki'], self.low_gains['kd']
        
        # NTC Constants
        self.V_ref, self.R_fixed, self.B_ntc = 1.8, 10000.0, 3950.0
        self.T0, self.R0 = 298.15, 10000.0

        # Bode Config
        self.test_duration = 600   
        self.f_min, self.f_max = 0.001, 1.0
        self.num_frequencies, self.injection_amp, self.fs = 50, 0.2, 40.0             
        
        # State & Buffers
        self.is_identifying = False
        self.integral, self.last_error = 0.0, 0.0
        self.start_time = time.time()
        self.max_points = 800 
        self.time_data = np.zeros(self.max_points)
        self.y_data = np.full(self.max_points, self.t_amb)
        self.amb_data = np.full(self.max_points, self.t_amb)
        self.u_data = np.zeros(self.max_points)
        self.u_history, self.y_history, self.d_history = [], [], []
        
        self.handle = None
        self.init_hardware()
        self.generate_multisine()
        self.init_ui()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.control_loop)
        if self.handle is not None:
            self.timer.start(25) # 40Hz
        else:
            print("CRITICAL: Device not found.")

    def init_hardware(self):
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            # Only AIN3 (Control) and AIN2 (Ambient)
            configs = [
                "AIN3_NEGATIVE_CH", 199, "AIN3_RANGE", 1.0, "AIN3_RESOLUTION_INDEX", 8,
                "AIN2_NEGATIVE_CH", 199, "AIN2_RANGE", 1.0, "AIN2_RESOLUTION_INDEX", 8,
                "DAC0", self.V_ref
            ]
            ljm.eWriteNames(self.handle, len(configs)//2, configs[::2], configs[1::2])
            
            # PWM Setup
            self.roll_value = int((80000000 / 256) / 200.0)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_DIVISOR", 256)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ROLL_VALUE", self.roll_value)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 1)
            ljm.eWriteName(self.handle, "DIO0_EF_INDEX", 0)
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 1)
            print("LabJack Connected (AIN3 + AIN2).")
        except Exception as e:
            print(f"Init Hardware Failed: {e}")
            self.handle = None

    def generate_multisine(self):
        freqs = np.logspace(np.log10(self.f_min), np.log10(self.f_max), self.num_frequencies)
        self.target_freqs = np.unique(np.round(freqs * self.test_duration) / self.test_duration)
        N = len(self.target_freqs)
        self.phases = np.array([np.pi * k**2 / N for k in range(N)])

    def voltage_to_temp(self, v_out):
        if v_out <= 0.005 or v_out >= (self.V_ref - 0.005): return float(self.t_amb)
        try:
            r_ntc = (v_out * self.R_fixed) / (self.V_ref - v_out)
            inv_t = (1.0 / self.T0) + (1.0 / self.B_ntc) * np.log(r_ntc / self.R0)
            return float((1.0 / inv_t) - 273.15)
        except: return float(self.t_amb)

    def init_ui(self):
        self.setWindowTitle("Thermal Bode Analyzer (AIN3 + Ambient)")
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QtWidgets.QVBoxLayout(self.central_widget)
        self.win = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.win)
        
        self.p1 = self.win.addPlot(title="Control Temperature (AIN3)")
        self.curve_sec = self.p1.plot(pen='c')
        self.p1.addItem(pg.InfiniteLine(pos=self.setpoint, angle=0, pen=pg.mkPen('g', style=QtCore.Qt.DashLine)))
        
        self.win.nextRow()
        self.p_amb = self.win.addPlot(title="Ambient Temperature (AIN2)")
        self.curve_amb = self.p_amb.plot(pen='m')
        self.p_amb.setXLink(self.p1)
        
        self.win.nextRow()
        self.p_u = self.win.addPlot(title="Heater Duty Cycle %")
        self.curve_u = self.p_u.plot(pen='y')
        self.p_u.setYRange(0, 100)
        self.p_u.setXLink(self.p1)

        ctrl_box = QtWidgets.QHBoxLayout()
        self.btn_bode = QtWidgets.QPushButton("Start Bode Scan")
        self.btn_bode.clicked.connect(self.toggle_bode)
        self.progress_bar = QtWidgets.QProgressBar()
        ctrl_box.addWidget(self.btn_bode); ctrl_box.addWidget(self.progress_bar)
        self.layout.addLayout(ctrl_box)

    def toggle_bode(self):
        if not self.is_identifying:
            self.u_history, self.y_history, self.d_history = [], [], []
            self.start_bode_time = time.time()
            self.is_identifying = True
            self.btn_bode.setText("Cancel Scan")
        else:
            self.is_identifying = False
            self.btn_bode.setText("Start Bode Scan")

    def control_loop(self):
        if self.handle is None:
            self.timer.stop()
            return
        try:
            # Only read AIN3 and AIN2
            res = ljm.eReadNames(self.handle, 2, ["AIN3", "AIN2"])
            y_sec = self.voltage_to_temp(res[0])
            y_amb = self.voltage_to_temp(res[1])
            
            now = time.time()
            dt, error = 0.025, self.setpoint - y_sec
            abs_err = abs(error)
            
            # Gain Switching Logic
            gains = self.high_gains if abs_err < self.threshold else self.low_gains
            self.kp, self.ki, self.kd = gains['kp'], gains['ki'], gains['kd']
            
            if abs_err < 4.0: self.integral = np.clip(self.integral + error * dt, -90, 90)
            u_total = (self.kp * error) + (self.ki * self.integral) + ((self.setpoint - y_amb)/self.k_thermal)

            if self.is_identifying:
                elapsed = now - self.start_bode_time
                sig = np.sum(np.cos(2 * np.pi * self.target_freqs * elapsed + self.phases))
                inj = (sig / np.sqrt(len(self.target_freqs))) * self.injection_amp
                final_u = np.clip(u_total + inj, 0, 1)
                self.u_history.append(final_u); self.y_history.append(y_sec); self.d_history.append(y_amb)
                self.progress_bar.setValue(int((elapsed / self.test_duration) * 100))
                if elapsed >= self.test_duration: self.finalize_scan()
            else:
                final_u = np.clip(u_total, 0, 1)

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
            self.handle = None
            self.timer.stop()

    def finalize_scan(self):
        self.is_identifying = False
        self.btn_bode.setText("Start Bode Scan")
        u, y, d = np.array(self.u_history), np.array(self.y_history), np.array(self.d_history)
        
        # Save Raw Data
        ts = int(time.time())
        with open(f"bode_time_{ts}.csv", 'w', newline='') as f:
            csv.writer(f).writerow(['u_duty', 'y_temp', 'd_amb'])
            csv.writer(f).writerows(zip(u, y, d))

        # Frequency Analysis
        y_clean = y - ((self.setpoint - d) / self.k_thermal)
        nper = int(self.fs * 150)
        f, Pxy = signal.csd(u, y_clean, fs=self.fs, nperseg=nper)
        _, Pxx = signal.welch(u, fs=self.fs, nperseg=nper)
        H = Pxy / Pxx
        coh = np.abs(Pxy)**2 / (Pxx * signal.welch(y_clean, fs=self.fs, nperseg=nper)[1])
        
        self.show_results(f, H, coh)

    def show_results(self, f, H, coh):
        self.res_win = pg.GraphicsLayoutWidget(title="Bode Result")
        p1 = self.res_win.addPlot(title="Mag (dB)"); p1.setLogMode(True, False)
        p1.plot(f[1:], 20*np.log10(np.abs(H[1:])), pen='y')
        self.res_win.nextRow()
        p2 = self.res_win.addPlot(title="Coherence"); p2.plot(f[1:], coh[1:], pen='m')
        self.res_win.show()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ctrl = ThermalController(); ctrl.show(); sys.exit(app.exec_())