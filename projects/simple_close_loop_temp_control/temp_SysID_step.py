import sys
import time
import numpy as np
import pyqtgraph as pg
import pandas as pd
from pyqtgraph.Qt import QtCore, QtWidgets
from labjack import ljm

class ThermalStepID(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        # --- Constants from Bode Code ---
        self.V_ref, self.R_fixed, self.B_ntc = 1.8, 10000.0, 3950.0
        self.T0, self.R0 = 298.15, 10000.0
        self.setpoint = 30.0 
        
        self.is_identifying = False
        self.step_magnitude = 0.35  
        self.sys_id_buffer = []
        self.max_points = 1200 
        
        self.handle = None
        self.init_hardware() # Call before UI to catch errors early
        self.init_ui()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.control_loop)
        if self.handle:
            self.timer.start(25) 

    def init_hardware(self):
        """Mixed Range Differential Setup - Robust Version"""
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            
            # Increase timeout for high-resolution differential reads
            ljm.writeLibraryConfigS("LJM_SEND_RECEIVE_TIMEOUT_MS", 5000)

            configs = [
                # AIN0-1: Ambient
                "AIN0_NEGATIVE_CH", 1, 
                "AIN0_RANGE", 1.0, 
                "AIN0_RESOLUTION_INDEX", 8,
                "AIN0_SETTLING_US", 10000, # Increased settling for differential
                "AIN1_RANGE", 0.01,
                
                # AIN2-3: Feedback
                "AIN2_NEGATIVE_CH", 3, 
                "AIN2_RANGE", 1.0, 
                "AIN2_RESOLUTION_INDEX", 8,
                "AIN2_SETTLING_US", 10000, 
                "AIN3_RANGE", 0.01,
                
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
            
            print(f"Hardware Connected. DAC0 set to {self.V_ref}V.")
        except Exception as e:
            print(f"CRITICAL HARDWARE ERROR: {e}")
            self.handle = None

    def voltage_to_temp(self, v_out):
        """Returns 25C if voltage is out of bounds to prevent math crashes"""
        if v_out <= 0.001 or v_out >= self.V_ref - 0.001:
            return 25.0 
        try:
            r_ntc = (v_out * self.R_fixed) / (self.V_ref - v_out)
            inv_t = (1.0 / self.T0) + (1.0 / self.B_ntc) * np.log(r_ntc / self.R0)
            return (1.0 / inv_t) - 273.15
        except:
            return 25.0

    def init_ui(self):
        self.setWindowTitle("Thermal System ID - Step Response")
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        layout = QtWidgets.QVBoxLayout(self.central_widget)
        
        self.btn_id = QtWidgets.QPushButton(f"Apply {self.step_magnitude*100}% Step")
        self.btn_id.clicked.connect(self.toggle_id)
        if not self.handle:
            self.btn_id.setEnabled(False)
            self.btn_id.setText("HARDWARE NOT FOUND")
        
        layout.addWidget(self.btn_id)

        self.win = pg.GraphicsLayoutWidget()
        layout.addWidget(self.win)
        self.p1 = self.win.addPlot(title="Step Response (C)"); self.p1.addLegend()
        self.curve_fb = self.p1.plot(pen='c', name="Feedback (AIN2)")
        self.curve_amb = self.p1.plot(pen='m', name="Ambient (AIN0)")
        self.win.nextRow()
        self.p2 = self.win.addPlot(title="Heater Input (%)"); self.curve_u = self.p2.plot(pen='y')
        
        self.time_data = np.zeros(self.max_points)
        self.y_data = np.full(self.max_points, 25.0) # Init to room temp
        self.a_data = np.full(self.max_points, 25.0)
        self.u_data = np.zeros(self.max_points)
        self.start_time = time.time()

    def toggle_id(self):
        if not self.is_identifying:
            self.sys_id_buffer = []
            self.is_identifying = True
            self.btn_id.setText("Recording... (Click to Stop & Save)")
        else:
            self.is_identifying = False
            self.btn_id.setText(f"Apply {self.step_magnitude*100}% Step")
            self.save_data()

    def save_data(self):
        if self.sys_id_buffer:
            df = pd.DataFrame(self.sys_id_buffer)
            filename = f"step_test_{int(time.time())}.csv"
            df.to_csv(filename, index=False)
            print(f"Step Response saved to {filename}")

    def control_loop(self):
        if not self.handle: return
        try:
            # Read both in one call. If it hangs here, it's a timeout/timing issue.
            res = ljm.eReadNames(self.handle, 2, ["AIN0", "AIN2"])
            y_amb = self.voltage_to_temp(res[0])
            y_fb = self.voltage_to_temp(res[1])
            now = time.time()

            if self.is_identifying:
                u_final = self.step_magnitude
                self.sys_id_buffer.append({'time': now, 'y': y_fb, 'd': y_amb, 'u': u_final})
            else:
                # Idle: very low power just to verify DAC/Thermistor circuit works
                u_final = 0.00 if y_fb < self.setpoint else 0.0

            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", int(u_final * self.roll_value))

            # Update Plot Arrays
            self.time_data = np.roll(self.time_data, -1); self.time_data[-1] = now - self.start_time
            self.y_data = np.roll(self.y_data, -1); self.y_data[-1] = y_fb
            self.a_data = np.roll(self.a_data, -1); self.a_data[-1] = y_amb
            self.u_data = np.roll(self.u_data, -1); self.u_data[-1] = u_final * 100
            
            self.curve_fb.setData(self.time_data, self.y_data)
            self.curve_amb.setData(self.time_data, self.a_data)
            self.curve_u.setData(self.time_data, self.u_data)
        except Exception as e:
            print(f"Runtime Loop Error: {e}")

    def closeEvent(self, event):
        if self.handle:
            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", 0); ljm.close(self.handle)
        event.accept()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ctrl = ThermalStepID()
    ctrl.show()
    sys.exit(app.exec_())