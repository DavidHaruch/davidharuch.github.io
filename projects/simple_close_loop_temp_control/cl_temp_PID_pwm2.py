import sys
import time
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets
from labjack import ljm
from scipy.ndimage import uniform_filter1d
import csv
import os

class ThermalController(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        # --- Controller Configuration ---
        self.setpoint = 30.0  
        self.t_amb = 23.0     
        self.k_thermal = 22.0 
        
        self.use_secondary_sensor = True  
        self.show_primary_plot = False 
        
        # --- Gain Scheduling + Hysteresis Configuration ---
        self.threshold = 0.2    # Base threshold (°C)
        self.hysteresis = 0.02  # Hysteresis band (°C)
        
        # High Gains (Precision Mode)
        self.high_gains = {'kp': 0.6, 'ki': 0.015, 'kd': 0.001}
        # Low Gains (Approach Mode)
        self.low_gains = {'kp': 0.2, 'ki': 0.005, 'kd': 0.0}

        self.current_mode = "LOW" 
        self.kp, self.ki, self.kd = self.low_gains['kp'], self.low_gains['ki'], self.low_gains['kd']
        
        # NTC Constants
        self.R_fixed = 10000.0  
        self.V_ref = 1.8        
        self.B_ntc = 3950.0     
        self.T0 = 298.15        
        self.R0 = 10000.0       

        self.integral = 0.0
        self.last_error = 0.0
        self.window_size = 1
        self.max_points = 400 

        # --- Logging Setup ---
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        self.log_filename = f"thermal_log_{timestamp}.csv"
        self.init_logger()
        
        # --- UI Setup ---
        self.setWindowTitle("Ultra-Precision Thermal Controller")
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QtWidgets.QVBoxLayout(self.central_widget)

        self.win = pg.GraphicsLayoutWidget(show=True)
        self.layout.addWidget(self.win)
        
        self.p1 = self.win.addPlot(title="Temperature Control Loop")
        self.p1.setLabel('left', 'Temp', units='°C')
        self.p1.showGrid(x=True, y=True)
        self.p1.addLegend()
        
        self.curve1 = self.p1.plot(pen=pg.mkPen('r', width=2), name="Primary (AIN0)")
        self.curve_ntc = self.p1.plot(pen=pg.mkPen('c', width=1.5), name="Secondary (AIN3)")
        
        if not self.show_primary_plot:
            self.curve1.hide()
        
        self.setpoint_line = pg.InfiniteLine(pos=self.setpoint, angle=0, pen=pg.mkPen('g', style=QtCore.Qt.DashLine))
        self.p1.addItem(self.setpoint_line)
        
        self.stats_text = pg.TextItem(anchor=(0,0))
        self.p1.addItem(self.stats_text)

        self.win.nextRow()
        self.p2 = self.win.addPlot(title="Heater Output (PWM Duty %)")
        self.p2.setLabel('left', 'Power', units='%')
        self.p2.setLabel('bottom', 'Time', units='s')
        self.p2.setYRange(0, 100)
        self.p2.showGrid(x=True, y=True)
        self.p2.setXLink(self.p1)
        self.curve2 = self.p2.plot(pen=pg.mkPen('y', width=1.5))

        self.time_data = np.zeros(self.max_points)
        self.temp_data_raw = np.full(self.max_points, self.t_amb)
        self.ntc_temp_data = np.full(self.max_points, self.t_amb) 
        self.duty_data = np.zeros(self.max_points)
        self.start_time = time.time()
        self.last_update_time = time.time()
        
        self.handle = None
        self.init_hardware()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.control_loop)
        if self.handle is not None:
            self.timer.start(50) 

    def init_logger(self):
        """Creates the CSV file and writes the header."""
        with open(self.log_filename, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Timestamp_Unix', 'Time_Elapsed_s', 'Primary_Temp_C', 
                'Secondary_Temp_C', 'Setpoint_C', 'Duty_Cycle_Pct', 
                'Control_Mode', 'KP', 'KI', 'KD'
            ])

    def init_hardware(self):
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            ljm.eWriteName(self.handle, "DAC0", self.V_ref)
            configs = [
                "AIN0_NEGATIVE_CH", 1, "AIN0_RANGE", 0.01, "AIN0_RESOLUTION_INDEX", 8,
                "AIN0_SETTLING_US", 10000, "AIN0_EF_INDEX", 22, "AIN0_EF_CONFIG_A", 1, "AIN0_EF_CONFIG_B", 60052,
                "AIN3_NEGATIVE_CH", 199, "AIN3_RANGE", 1.0, "AIN3_RESOLUTION_INDEX", 8
            ]
            ljm.eWriteNames(self.handle, len(configs)//2, configs[::2], configs[1::2])
            
            self.pwm_hz = 200.0
            self.base_clock = 80000000
            self.divisor = 256
            self.roll_value = int((self.base_clock / self.divisor) / self.pwm_hz)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_DIVISOR", self.divisor)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ROLL_VALUE", self.roll_value)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 1)
            ljm.eWriteName(self.handle, "DIO0_EF_INDEX", 0) 
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 1)
        except ljm.LJMError as e:
            print(f"Hardware Error: {e}")
            self.handle = None

    def voltage_to_temp(self, v_out):
        if v_out <= 0.05 or v_out >= (self.V_ref - 0.05): return float(self.t_amb)
        try:
            r_ntc = (v_out * self.R_fixed) / (self.V_ref - v_out)
            inv_t = (1.0 / self.T0) + (1.0 / self.B_ntc) * np.log(r_ntc / self.R0)
            return float((1.0 / inv_t) - 273.15)
        except Exception: return float(self.t_amb)

    def control_loop(self):
        if self.handle is None: return
        try:
            results = ljm.eReadNames(self.handle, 2, ["AIN0_EF_READ_A", "AIN3"])
            temp_raw = float(results[0])
            secondary_temp = self.voltage_to_temp(float(results[1]))
            
            now = time.time()
            dt = now - self.last_update_time
            self.last_update_time = now

            self.temp_data_raw = np.roll(self.temp_data_raw, -1)
            self.temp_data_raw[-1] = temp_raw
            primary_filtered_arr = uniform_filter1d(self.temp_data_raw, size=self.window_size, mode='nearest')
            primary_filtered_val = primary_filtered_arr[-1]
            
            self.ntc_temp_data = np.roll(self.ntc_temp_data, -1)
            self.ntc_temp_data[-1] = secondary_temp
            
            current_temp = secondary_temp if self.use_secondary_sensor else primary_filtered_val
            stability_data = self.ntc_temp_data if self.use_secondary_sensor else primary_filtered_arr
            abs_err = abs(float(self.setpoint - current_temp))
            derivative = (float(self.setpoint - current_temp) - self.last_error) / dt

            # --- Hysteresis-Aware Gain Mode Switching ---
            new_mode = self.current_mode
            if self.current_mode == "LOW":
                if abs_err < (self.threshold - self.hysteresis):
                    new_mode = "HIGH"
            else: 
                if abs_err > (self.threshold + self.hysteresis):
                    new_mode = "LOW"

            # --- Bumpless Transfer ---
            if new_mode != self.current_mode:
                old_output = (self.kp * (self.setpoint - current_temp)) + (self.ki * self.integral) + (self.kd * derivative)
                gains = self.high_gains if new_mode == "HIGH" else self.low_gains
                self.kp, self.ki, self.kd = gains['kp'], gains['ki'], gains['kd']
                if self.ki != 0:
                    self.integral = (old_output - (self.kp * (self.setpoint - current_temp)) - (self.kd * derivative)) / self.ki
                self.current_mode = new_mode

            # --- PID Calculation ---
            error = float(self.setpoint - current_temp)
            if abs(error) < 4.0:
                self.integral += error * dt
            else:
                self.integral = 0.0 
            
            self.integral = np.clip(self.integral, -90, 90) 
            self.last_error = error
            
            ff_output = max(0.0, (self.setpoint - self.t_amb) / self.k_thermal)
            output = (self.kp * error) + (self.ki * self.integral) + (self.kd * derivative) + ff_output
            duty_percent = np.clip(output * 100, 0.0, 100.0)
            
            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", int((duty_percent / 100.0) * self.roll_value))

            # --- Logging Data ---
            elapsed = now - self.start_time
            with open(self.log_filename, mode='a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    now, f"{elapsed:.3f}", f"{primary_filtered_val:.6f}", 
                    f"{secondary_temp:.6f}", self.setpoint, f"{duty_percent:.4f}", 
                    self.current_mode, self.kp, self.ki, self.kd
                ])

            # --- UI Buffers and Plotting ---
            self.time_data = np.roll(self.time_data, -1)
            self.time_data[-1] = elapsed
            self.duty_data = np.roll(self.duty_data, -1)
            self.duty_data[-1] = duty_percent

            if self.show_primary_plot:
                self.curve1.setData(self.time_data.copy(), primary_filtered_arr.copy())
            self.curve_ntc.setData(self.time_data.copy(), self.ntc_temp_data.copy())
            self.curve2.setData(self.time_data.copy(), self.duty_data.copy())
            
            pk_pk = (np.max(stability_data) - np.min(stability_data)) * 1000
            self.stats_text.setHtml(f'<div style="color: #00FF00; font-size: 11pt;">'
                                    f'Setpoint: {self.setpoint}°C | Mode: {self.current_mode}<br>'
                                    f'AIN0: {primary_filtered_val:.4f}°C<br>'
                                    f'AIN3 (NTC): {secondary_temp:.2f}°C<br>'
                                    f'Pk-Pk (Active): {pk_pk:.1f} m°C</div>')
            self.stats_text.setPos(self.time_data[0], current_temp)

        except Exception as e:
            print(f"Loop Error: {e}")
            self.timer.stop()

    def closeEvent(self, event):
        if self.handle:
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 0)
            ljm.eWriteName(self.handle, "DAC0", 0)
            ljm.close(self.handle)
        event.accept()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ctrl = ThermalController()
    ctrl.show()
    sys.exit(app.exec_())