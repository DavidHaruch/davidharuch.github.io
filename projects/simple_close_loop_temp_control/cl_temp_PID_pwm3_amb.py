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
        
        # --- Gain Scheduling + Hysteresis ---
        self.threshold = 0.2
        self.hysteresis = 0.02
        self.high_gains = {'kp': 0.6, 'ki': 0.015, 'kd': 0.001}
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
        
        # PLOT 1: Control Loop
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
        # PLOT 2: Power Output
        self.p2 = self.win.addPlot(title="Heater Output (PWM Duty %)")
        self.p2.setLabel('left', 'Power', units='%')
        self.p2.setYRange(0, 100)
        self.p2.showGrid(x=True, y=True)
        self.p2.setXLink(self.p1)
        self.curve2 = self.p2.plot(pen=pg.mkPen('y', width=1.5))

        self.win.nextRow()
        # PLOT 3: Ambient Monitor (New)
        self.p3 = self.win.addPlot(title="Ambient Temperature (AIN2)")
        self.p3.setLabel('left', 'Amb', units='°C')
        self.p3.setLabel('bottom', 'Time', units='s')
        self.p3.showGrid(x=True, y=True)
        self.p3.setXLink(self.p1) # Sync zoom/pan with other plots
        self.curve_amb = self.p3.plot(pen=pg.mkPen('m', width=1.5), name="Ambient")

        self.time_data = np.zeros(self.max_points)
        self.temp_data_raw = np.full(self.max_points, self.t_amb)
        self.ntc_temp_data = np.full(self.max_points, self.t_amb) 
        self.amb_temp_data = np.full(self.max_points, self.t_amb) 
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
        with open(self.log_filename, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Timestamp_Unix', 'Elapsed_s', 'Primary_C', 'Secondary_C', 'Ambient_C', 'Setpoint_C', 'Duty_Pct', 'Mode'])

    def init_hardware(self):
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            ljm.eWriteName(self.handle, "DAC0", self.V_ref)
            configs = [
                "AIN0_NEGATIVE_CH", 1, "AIN0_RANGE", 0.01, "AIN0_RESOLUTION_INDEX", 8,
                "AIN3_NEGATIVE_CH", 199, "AIN3_RANGE", 1.0, "AIN3_RESOLUTION_INDEX", 8,
                "AIN2_NEGATIVE_CH", 199, "AIN2_RANGE", 1.0, "AIN2_RESOLUTION_INDEX", 8
            ]
            ljm.eWriteNames(self.handle, len(configs)//2, configs[::2], configs[1::2])
            
            # PWM Setup
            self.pwm_hz, self.divisor = 200.0, 256
            self.roll_value = int((80000000 / self.divisor) / self.pwm_hz)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_DIVISOR", self.divisor)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ROLL_VALUE", self.roll_value)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 1)
            ljm.eWriteName(self.handle, "DIO0_EF_INDEX", 0) 
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 1)
        except Exception as e:
            print(f"HW Error: {e}"); self.handle = None

    def voltage_to_temp(self, v_out):
        if v_out <= 0.05 or v_out >= (self.V_ref - 0.05): return float(self.t_amb)
        try:
            r_ntc = (v_out * self.R_fixed) / (self.V_ref - v_out)
            inv_t = (1.0 / self.T0) + (1.0 / self.B_ntc) * np.log(r_ntc / self.R0)
            return float((1.0 / inv_t) - 273.15)
        except: return float(self.t_amb)

    def control_loop(self):
        if self.handle is None: return
        try:
            results = ljm.eReadNames(self.handle, 3, ["AIN0_EF_READ_A", "AIN3", "AIN2"])
            temp_raw, secondary_temp, ambient_temp = float(results[0]), self.voltage_to_temp(results[1]), self.voltage_to_temp(results[2])
            
            now = time.time()
            dt, self.last_update_time = now - self.last_update_time, now
            elapsed = now - self.start_time

            # Shift buffers
            self.temp_data_raw = np.roll(self.temp_data_raw, -1); self.temp_data_raw[-1] = temp_raw
            self.ntc_temp_data = np.roll(self.ntc_temp_data, -1); self.ntc_temp_data[-1] = secondary_temp
            self.amb_temp_data = np.roll(self.amb_temp_data, -1); self.amb_temp_data[-1] = ambient_temp
            self.time_data = np.roll(self.time_data, -1); self.time_data[-1] = elapsed
            
            primary_filtered_val = uniform_filter1d(self.temp_data_raw, size=self.window_size, mode='nearest')[-1]
            current_temp = secondary_temp if self.use_secondary_sensor else primary_filtered_val
            
            # --- PID & Mode Switching ---
            error = float(self.setpoint - current_temp)
            abs_err = abs(error)
            derivative = (error - self.last_error) / dt
            
            if self.current_mode == "LOW" and abs_err < (self.threshold - self.hysteresis):
                self.current_mode = "HIGH"
            elif self.current_mode == "HIGH" and abs_err > (self.threshold + self.hysteresis):
                self.current_mode = "LOW"

            gains = self.high_gains if self.current_mode == "HIGH" else self.low_gains
            self.kp, self.ki, self.kd = gains['kp'], gains['ki'], gains['kd']
            
            if abs_err < 4.0: self.integral = np.clip(self.integral + error * dt, -90, 90)
            else: self.integral = 0.0
            
            ff_output = max(0.0, (self.setpoint - ambient_temp) / self.k_thermal)
            duty_percent = np.clip(((self.kp * error) + (self.ki * self.integral) + (self.kd * derivative) + ff_output) * 100, 0, 100)
            
            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", int((duty_percent / 100.0) * self.roll_value))
            self.last_error = error

            # --- Plotting ---
            self.duty_data = np.roll(self.duty_data, -1); self.duty_data[-1] = duty_percent
            self.curve_ntc.setData(self.time_data, self.ntc_temp_data)
            self.curve2.setData(self.time_data, self.duty_data)
            self.curve_amb.setData(self.time_data, self.amb_temp_data)
            
            self.stats_text.setHtml(f'<div style="color: #00FF00; font-size: 10pt;">'
                                    f'AIN3: {secondary_temp:.3f}°C | AMB: {ambient_temp:.2f}°C<br>'
                                    f'Mode: {self.current_mode} | Duty: {duty_percent:.1f}%</div>')
            self.stats_text.setPos(self.time_data[0], current_temp)

        except Exception as e:
            print(f"Loop Error: {e}"); self.timer.stop()

    def closeEvent(self, event):
        if self.handle:
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 0)
            ljm.close(self.handle)
        event.accept()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ctrl = ThermalController(); ctrl.show()
    sys.exit(app.exec_())