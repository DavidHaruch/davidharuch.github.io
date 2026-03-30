import sys
import time
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets
from labjack import ljm
import csv
from collections import deque

class ThermalController(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        # --- 1. Identified Plant Parameters ---
        self.Kp_plant = 10.7977
        self.tau_plant = 30.0
        self.L_plant = 1.76
        self.fs = 40.0
        self.dt = 1.0 / self.fs

        # --- 2. Controller & Ramping Configuration ---
        self.target_setpoint = 30.0    
        self.ramp_rate = 10.0          # Celsius per second
        self.t_amb_init = 23.0
        self.k_thermal = 19.49 
        
        self.current_setpoint = self.t_amb_init 
        
        # Gains
        self.sp_gains = {'kp': 5.55674, 'ki': 0.18522}
        # self.std_gains = {'kp': 5.55674, 'ki': 0.18522}
        self.std_gains = {'kp': 0.8523, 'ki': 0.0284}
        self.id_gains = {'kp': 0.05, 'ki': 0.01}
        
        # --- 3. Smith Predictor Math (Tustin Transform) ---
        den = (2 * self.tau_plant + self.dt)
        self.a1 = (2 * self.tau_plant - self.dt) / den
        self.b0 = (self.Kp_plant * self.dt) / den
        self.b1 = self.b0
        
        self.buffer_len = int(np.ceil(self.L_plant * self.fs))
        self.u_model_buffer = deque([self.t_amb_init] * self.buffer_len, maxlen=self.buffer_len)
        self.y_model_instant = self.t_amb_init
        self.u_last_shadow = 0.0

        # Hardware/NTC Constants
        self.V_ref, self.R_fixed, self.B_ntc = 1.8, 10000.0, 3950.0
        self.T0, self.R0 = 298.15, 10000.0

        # Sine Sweep Config
        self.test_duration = 2400
        self.f_start, self.f_end = 0.0005, 0.1
        self.injection_amp = 0.2    
        
        # --- 4. Closed-Loop Step Response Config ---
        self.is_logging_step = False
        self.step_log_duration = 300 # 5 minutes
        self.step_start_time = 0
        
        # --- State Variables ---
        self.is_identifying = False
        self.use_smith = True
        self.integral = 0.0
        self.start_time = time.time()
        self.start_bode_time = 0
        
        # Buffers for Live Plotting
        self.max_points = 1000 
        self.time_data = np.zeros(self.max_points)
        self.y_data = np.full(self.max_points, self.t_amb_init)
        self.amb_data = np.full(self.max_points, self.t_amb_init)
        self.u_data = np.zeros(self.max_points)
        self.sp_line_data = np.full(self.max_points, self.t_amb_init)
        
        # Data Logging Buffers
        self.log_history = [] # Stores [time, setpoint, y_sec, y_amb, u, error, integral, is_smith]
        
        self.handle = None
        self.init_hardware()
        self.init_ui()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.control_loop)
        if self.handle:
            self.timer.start(25) # 40Hz

    def init_hardware(self):
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            configs = [
                "AIN0_NEGATIVE_CH", 1, "AIN0_RANGE", 1.0, "AIN0_RESOLUTION_INDEX", 8,
                "AIN0_SETTLING_US", 5000,
                "AIN2_NEGATIVE_CH", 3, "AIN2_RANGE", 1.0, "AIN2_RESOLUTION_INDEX", 8,
                "AIN2_SETTLING_US", 5000,
                "DAC0", self.V_ref
            ]
            ljm.eWriteNames(self.handle, len(configs)//2, configs[::2], configs[1::2])
            
            self.roll_value = int((80000000 / 256) / 5.0)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_DIVISOR", 256)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ROLL_VALUE", self.roll_value)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 1)
            ljm.eWriteName(self.handle, "DIO0_EF_INDEX", 0)
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 1)
            print("Hardware Initialized.")
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
            res = ljm.eReadNames(self.handle, 2, ["AIN0", "AIN2"])
            y_amb = self.voltage_to_temp(res[0])
            y_sec = self.voltage_to_temp(res[1])
            now = time.time()
            t_run = now - self.start_time
            
            # --- 1. Setpoint Handling ---
            if self.is_logging_step:
                # Instant jump, no ramping
                self.current_setpoint = self.target_setpoint
                t_step = now - self.step_start_time
                self.progress_bar.setValue(int((t_step / self.step_log_duration) * 100))
                if t_step >= self.step_log_duration:
                    self.finalize_log("step_response")
            else:
                # Standard Ramping
                if self.current_setpoint < self.target_setpoint:
                    self.current_setpoint = min(self.target_setpoint, 
                                                self.current_setpoint + (self.ramp_rate * self.dt))
                elif self.current_setpoint > self.target_setpoint:
                    self.current_setpoint = max(self.target_setpoint, 
                                                self.current_setpoint - (self.ramp_rate * self.dt))

            # --- 2. Smith Predictor Math ---
            y_instant_new = (self.a1 * self.y_model_instant) + \
                            (self.b0 * self.u_last_shadow) + \
                            (self.b1 * self.u_last_shadow)
            
            self.u_model_buffer.append(y_instant_new)
            y_delayed = self.u_model_buffer[0]
            y_pred_correction = y_instant_new - y_delayed
            self.y_model_instant = y_instant_new

            # --- 3. Controller Logic ---
            error_raw = self.current_setpoint - y_sec

            if self.is_identifying:
                kp, ki = self.id_gains['kp'], self.id_gains['ki']
                error = error_raw
            elif self.use_smith:
                kp, ki = self.sp_gains['kp'], self.sp_gains['ki']
                error = self.current_setpoint - (y_sec + y_pred_correction)
            else:
                kp, ki = self.std_gains['kp'], self.std_gains['ki']
                error = error_raw

            if abs(error_raw) < 4.0: 
                self.integral = np.clip(self.integral + error * self.dt, -100, 100)
            else:
                self.integral = 0.0
            
            u_base = (kp * error) + (ki * self.integral) + ((self.current_setpoint - y_amb)/self.k_thermal)

            # Sine Sweep Injection
            inj = 0.0
            if self.is_identifying:
                t_bode = now - self.start_bode_time
                if t_bode < self.test_duration:
                    exponent = t_bode / self.test_duration
                    freq_ratio = self.f_end / self.f_start
                    phi = 2 * np.pi * self.f_start * (self.test_duration / np.log(freq_ratio)) * (freq_ratio**exponent - 1)
                    inj = self.injection_amp * np.sin(phi)
                    self.progress_bar.setValue(int((t_bode / self.test_duration) * 100))
                else:
                    self.finalize_log("sine_sweep")

            final_u = np.clip(u_base + inj, 0, 1)
            self.u_last_shadow = final_u 
            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", int(final_u * self.roll_value))

            # --- 4. Data Logging ---
            # Capture state for CSV
            if self.is_identifying or self.is_logging_step:
                self.log_history.append([
                    t_run, self.current_setpoint, y_sec, y_amb, final_u, 
                    error, self.integral, int(self.use_smith), kp, ki
                ])

            # --- 5. UI Updates ---
            self.time_data = np.roll(self.time_data, -1); self.time_data[-1] = t_run
            self.y_data = np.roll(self.y_data, -1); self.y_data[-1] = y_sec
            self.amb_data = np.roll(self.amb_data, -1); self.amb_data[-1] = y_amb
            self.u_data = np.roll(self.u_data, -1); self.u_data[-1] = final_u * 100
            self.sp_line_data = np.roll(self.sp_line_data, -1); self.sp_line_data[-1] = self.current_setpoint
            
            self.curve_sec.setData(self.time_data, self.y_data)
            self.curve_amb.setData(self.time_data, self.amb_data)
            self.curve_u.setData(self.time_data, self.u_data)
            self.setpoint_line.setData(self.time_data, self.sp_line_data)

        except Exception as e:
            print(f"Loop Error: {e}"); self.cleanup_hardware()

    def init_ui(self):
        self.setWindowTitle("Thermal Controller - Advanced Diagnostics")
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        layout = QtWidgets.QVBoxLayout(self.central_widget)
        self.win = pg.GraphicsLayoutWidget()
        layout.addWidget(self.win)
        
        self.p1 = self.win.addPlot(title="Feedback Temp (AIN2-3)")
        self.curve_sec = self.p1.plot(pen='c', name="Feedback")
        self.setpoint_line = self.p1.plot(pen=pg.mkPen('g', style=QtCore.Qt.DashLine), name="Setpoint")
        
        self.win.nextRow()
        self.p_amb = self.win.addPlot(title="Ambient Temp (AIN0-1)"); self.curve_amb = self.p_amb.plot(pen='m')
        self.p_amb.setXLink(self.p1)
        self.win.nextRow()
        self.p_u = self.win.addPlot(title="Heater Duty %"); self.curve_u = self.p_u.plot(pen='y')
        self.p_u.setYRange(0, 100); self.p_u.setXLink(self.p1)

        # Control Panel
        controls = QtWidgets.QHBoxLayout()
        
        # ID Group
        id_group = QtWidgets.QGroupBox("Identification & Steps")
        id_layout = QtWidgets.QHBoxLayout()
        
        self.btn_bode = QtWidgets.QPushButton("Start Sine Sweep")
        self.btn_bode.clicked.connect(self.toggle_bode)
        
        self.spin_step_target = QtWidgets.QDoubleSpinBox()
        self.spin_step_target.setRange(20, 60); self.spin_step_target.setValue(35.0)
        
        self.btn_step = QtWidgets.QPushButton("Execute Step SP")
        self.btn_step.clicked.connect(self.trigger_step_response)
        
        id_layout.addWidget(self.btn_bode)
        id_layout.addWidget(QtWidgets.QLabel("Step To:"))
        id_layout.addWidget(self.spin_step_target)
        id_layout.addWidget(self.btn_step)
        id_group.setLayout(id_layout)

        # Settings Group
        set_group = QtWidgets.QGroupBox("Controller Settings")
        set_layout = QtWidgets.QHBoxLayout()
        self.chk_smith = QtWidgets.QCheckBox("Use Smith Predictor"); self.chk_smith.setChecked(True)
        self.chk_smith.stateChanged.connect(self.update_smith_toggle)
        set_layout.addWidget(self.chk_smith)
        set_group.setLayout(set_layout)

        controls.addWidget(id_group)
        controls.addWidget(set_group)
        layout.addLayout(controls)
        
        self.progress_bar = QtWidgets.QProgressBar()
        layout.addWidget(self.progress_bar)

    def trigger_step_response(self):
        if not self.is_logging_step:
            self.log_history = []
            self.target_setpoint = self.spin_step_target.value()
            self.current_setpoint = self.target_setpoint # Snap to target
            self.step_start_time = time.time()
            self.is_logging_step = True
            self.btn_step.setText("Cancel Step")
            print(f"Recording closed-loop step to {self.target_setpoint} C")
        else:
            self.is_logging_step = False
            self.btn_step.setText("Execute Step SP")

    def toggle_bode(self):
        if not self.is_identifying:
            self.log_history = []
            self.start_bode_time = time.time(); self.is_identifying = True
            self.btn_bode.setText("Cancel Sweep")
        else:
            self.is_identifying = False; self.btn_bode.setText("Start Sine Sweep")

    def update_smith_toggle(self, state):
        self.use_smith = (state == QtCore.Qt.Checked)
        self.integral = 0.0

    def finalize_log(self, prefix):
        self.is_identifying = False
        self.is_logging_step = False
        self.btn_bode.setText("Start Sine Sweep")
        self.btn_step.setText("Execute Step SP")
        
        if not self.log_history: return
        
        ts = int(time.time())
        filename = f"{prefix}_{ts}.csv"
        header = ['time_s', 'setpoint', 'temp_sec', 'temp_amb', 'duty', 'error', 'integral', 'is_smith', 'kp', 'ki']
        
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(self.log_history)
        
        print(f"Log Saved: {filename}")
        self.log_history = []

    def cleanup_hardware(self):
        if self.handle:
            try:
                self.timer.stop()
                ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", 0)
                ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 0)
                ljm.close(self.handle); self.handle = None
            except: pass

    def closeEvent(self, event):
        self.cleanup_hardware(); event.accept()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ctrl = ThermalController()
    try:
        ctrl.show()
        sys.exit(app.exec_())
    finally:
        ctrl.cleanup_hardware()