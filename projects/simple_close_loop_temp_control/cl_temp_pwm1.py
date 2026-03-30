import sys
import time
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets
from labjack import ljm
from scipy.ndimage import uniform_filter1d

class ThermalMonitor(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        # --- Filter & Stats Configuration ---
        self.window_size = 10  # Moving average window
        self.max_points = 500  # Number of points shown on screen
        
        # --- UI Setup ---
        self.setWindowTitle("Precision Thermal Metrology")
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QtWidgets.QVBoxLayout(self.central_widget)

        self.win = pg.GraphicsLayoutWidget(show=True)
        self.layout.addWidget(self.win)
        
        # Top Plot: Temperature
        self.p1 = self.win.addPlot(title="Filtered Temperature Stability")
        self.p1.setLabel('left', 'Temp', units='°C')
        self.p1.showGrid(x=True, y=True)
        self.curve1 = self.p1.plot(pen=pg.mkPen('r', width=2))
        
        # Stats Label (Top Left of Plot)
        self.stats_text = pg.TextItem(html='<div style="color: #FFF; background-color: rgba(0,0,0,150); padding: 5px;">'
                                           'Initialising...</div>', anchor=(0,0))
        self.p1.addItem(self.stats_text)

        self.win.nextRow()

        # Bottom Plot: PWM Monitor
        self.p2 = self.win.addPlot(title="PWM Load Monitor")
        self.p2.setLabel('left', 'Voltage', units='V')
        self.p2.setLabel('bottom', 'Time', units='s')
        self.p2.showGrid(x=True, y=True)
        self.p2.setXLink(self.p1)
        self.curve2 = self.p2.plot(pen=pg.mkPen('b', width=1.5))

        # --- Data Buffers ---
        self.time_data = np.zeros(self.max_points)
        self.temp_data_raw = np.zeros(self.max_points)
        self.pwm_data = np.zeros(self.max_points)
        self.start_time = time.time()
        
        # --- Hardware Initialization ---
        self.handle = None
        self.init_hardware()

        # --- Timer Setup (30Hz) ---
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_plot)
        if self.handle is not None:
            self.timer.start(33) 
        else:
            print("Check LabJack connection.")

    def init_hardware(self):
        try:
            self.handle = ljm.openS("T7", "ANY", "ANY")
            configs = [
                "AIN0_NEGATIVE_CH", 1,
                "AIN0_RANGE", 0.01,
                "AIN0_RESOLUTION_INDEX", 8,
                "AIN0_SETTLING_US", 10000, # Increased settling for metrology
                "AIN2_RANGE", 10.0,
                "AIN0_EF_INDEX", 22,
                "AIN0_EF_CONFIG_A", 1,        
                "AIN0_EF_CONFIG_B", 60052     
            ]
            ljm.eWriteNames(self.handle, len(configs)//2, configs[::2], configs[1::2])
            
            # PWM Setup (0.5Hz)
            roll = int((80000000 / 256) / 0.5)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_DIVISOR", 256)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ROLL_VALUE", roll)
            ljm.eWriteName(self.handle, "DIO_EF_CLOCK0_ENABLE", 1)
            ljm.eWriteName(self.handle, "DIO0_EF_INDEX", 0)
            ljm.eWriteName(self.handle, "DIO0_EF_CONFIG_A", int(roll * 0.5))
            ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 1)
            
        except ljm.LJMError as e:
            print(f"Hardware Init Error: {e}")
            self.handle = None

    def update_plot(self):
        if self.handle is None: return

        try:
            # 1. Acquire
            temp_raw = ljm.eReadName(self.handle, "AIN0_EF_READ_A")
            pwm = ljm.eReadName(self.handle, "AIN2")
            now = time.time() - self.start_time

            # 2. Update Buffers
            self.time_data = np.roll(self.time_data, -1)
            self.time_data[-1] = now
            self.temp_data_raw = np.roll(self.temp_data_raw, -1)
            self.temp_data_raw[-1] = temp_raw
            self.pwm_data = np.roll(self.pwm_data, -1)
            self.pwm_data[-1] = pwm

            # 3. Filter (Uniform Filter handles edges by reflecting data)
            filtered_temp = uniform_filter1d(self.temp_data_raw, size=self.window_size, mode='nearest')

            # 4. Calculate Peak-to-Peak Noise (on filtered data)
            # Only calculate on the visible window (ignore the initial zeros)
            valid_data = filtered_temp[self.window_size:]
            pk_pk = (np.max(valid_data) - np.min(valid_data)) * 1000 # Convert to mC
            
            # Update Stats Label
            self.stats_text.setHtml(
                f'<div style="color: #00FF00; font-family: monospace; font-size: 12pt;">'
                f'Temp: {filtered_temp[-1]:.3f} °C<br>'
                f'Pk-Pk: {pk_pk:.1f} m°C</div>'
            )
            self.stats_text.setPos(self.time_data[0], np.max(filtered_temp))

            # 5. Push to GPU
            self.curve1.setData(self.time_data, filtered_temp)
            self.curve2.setData(self.time_data, self.pwm_data)
            
        except ljm.LJMError as e:
            if e.errorCode == 1224: self.timer.stop()

    def closeEvent(self, event):
        if self.handle:
            try:
                ljm.eWriteName(self.handle, "DIO0_EF_ENABLE", 0)
                ljm.close(self.handle)
            except: pass
        event.accept()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    monitor = ThermalMonitor()
    monitor.show()
    sys.exit(app.exec_())