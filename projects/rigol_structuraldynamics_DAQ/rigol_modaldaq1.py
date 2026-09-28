import sys
import numpy as np
import scipy.signal as signal
import pyvisa

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QDoubleSpinBox,
    QGroupBox, QTextEdit, QProgressBar, QFileDialog, QMessageBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


# ==============================================================================
# 1. HARDWARE COMMUNICATION ENGINE
# ==============================================================================
class RigolMHO934:
    def __init__(self, ip_address):
        self.rm = pyvisa.ResourceManager()
        self.resource_name = f"TCPIP0::{ip_address}::INSTR"
        self.scope = self.rm.open_resource(self.resource_name)
        self.scope.timeout = 5000
        self.scope.write("*CLS")
        
    def fetch_channel_data(self, channel_num):
        self.scope.write(f":WAVeform:SOURce CHANnel{channel_num}")
        self.scope.write(":WAVeform:MODE RAW")
        self.scope.write(":WAVeform:FORMat BYTE")
        
        preamble = self.scope.query(":WAVeform:PREamble?").split(',')
        x_incr = float(preamble[4])
        x_orig = float(preamble[5])
        y_incr = float(preamble[7])
        y_orig = float(preamble[8])
        y_ref = float(preamble[9])
        
        raw_data = self.scope.query_binary_values(":WAVeform:DATA?", datatype='B', container=np.array)
        volts = (raw_data - y_orig - y_ref) * y_incr
        time = np.arange(len(volts)) * x_incr + x_orig
        fs = 1.0 / x_incr if x_incr > 0 else 1.0
        return time, volts, fs

    def close(self):
        try:
            self.scope.close()
        except Exception:
            pass


# ==============================================================================
# 2. SIGNAL PROCESSING ENGINE
# ==============================================================================
class ModalSignalProcessor:
    def __init__(self, fs, exp_decay_rate=0.05):
        self.fs = fs
        self.exp_decay_rate = exp_decay_rate

    def apply_windows(self, force, resp):
        N = len(force)
        t = np.arange(N) / self.fs
        
        peak_idx = np.argmax(np.abs(force))
        force_win = np.zeros(N)
        win_len = int(0.02 * self.fs)
        start = max(0, peak_idx - int(0.002 * self.fs))
        end = min(N, start + win_len)
        force_win[start:end] = 1.0
        
        resp_win = np.exp(-self.exp_decay_rate * t)
        return force * force_win, resp * resp_win

    def compute_average_frf(self, force_runs, resp_runs):
        Sff_total = 0.0
        Sxf_total = 0.0
        Sxx_total = 0.0
        num_avg = len(force_runs)
        
        for i in range(num_avg):
            f, Sff = signal.welch(force_runs[i], fs=self.fs, window='boxcar', nperseg=len(force_runs[i]), scaling='spectrum')
            _, Sxx = signal.welch(resp_runs[i], fs=self.fs, window='boxcar', nperseg=len(resp_runs[i]), scaling='spectrum')
            _, Sxf = signal.csd(force_runs[i], resp_runs[i], fs=self.fs, window='boxcar', nperseg=len(force_runs[i]), scaling='spectrum')
            
            Sff_total += Sff
            Sxx_total += Sxx
            Sxf_total += Sxf

        Sff_avg = Sff_total / num_avg
        Sxx_avg = Sxx_total / num_avg
        Sxf_avg = Sxf_total / num_avg

        H1 = Sxf_avg / Sff_avg
        coherence = (np.abs(Sxf_avg)**2) / (Sff_avg * Sxx_avg)
        return f, H1, coherence


# ==============================================================================
# 3. DATA MODEL & UFF58 EXPORTER
# ==============================================================================
class RovingModalSession:
    def __init__(self, impact_node=1):
        self.impact_node = impact_node
        self.measurements = {}

    def add_point_data(self, resp_node, freq, H1_complex, coherence):
        self.measurements[str(resp_node)] = {
            "freq": freq.tolist(),
            "real": np.real(H1_complex).tolist(),
            "imag": np.imag(H1_complex).tolist(),
            "mag_dB": (20 * np.log10(np.abs(H1_complex) + 1e-12)).tolist(),
            "phase_deg": np.angle(H1_complex, deg=True).tolist(),
            "coherence": coherence.tolist()
        }

    def export_to_unv58(self, filename):
        with open(filename, "w") as f:
            for resp_node, data in self.measurements.items():
                freqs = data["freq"]
                reals = data["real"]
                imags = data["imag"]
                num_pts = len(freqs)
                f_start = freqs[0]
                f_inc = freqs[1] - freqs[0] if num_pts > 1 else 0.0

                f.write("    -1\n")
                f.write("    58\n")
                f.write(f"Roving Modal Test Pt {resp_node}\n")
                f.write("NONE\n")
                f.write("NONE\n")
                f.write(f"         1         1{int(self.impact_node):10d}         3         1{int(resp_node):10d}         3\n")
                f.write(f"         2        18         2{num_pts:10d}         1{f_start:13.5e}{f_inc:13.5e} 0.00000e+00\n")
                f.write("         0         0         0         0         0         0\n")
                
                for r, i in zip(reals, imags):
                    f.write(f"{r:13.5e}{i:13.5e}\n")
                f.write("    -1\n")


# ==============================================================================
# 4. VISA WORKER THREAD
# ==============================================================================
class FetchWorker(QThread):
    finished = pyqtSignal(np.ndarray, np.ndarray, np.ndarray, float)
    error = pyqtSignal(str)

    def __init__(self, ip, use_sim, parent=None):
        super().__init__(parent)
        self.ip = ip
        self.use_sim = use_sim

    def run(self):
        if self.use_sim:
            # Synthetic voltage signals assuming default ~100 mV/g and ~10 mV/N sensors
            fs = 10000.0
            N = 8192
            t = np.arange(N) / fs
            force_volts = 0.5 * np.exp(-((t - 0.005) / 0.0008)**2) + 0.002 * np.random.randn(N)
            resp_volts = 0.08 * np.exp(-12 * t) * np.sin(2 * np.pi * 180 * t) + \
                         0.04 * np.exp(-25 * t) * np.sin(2 * np.pi * 420 * t) + \
                         0.001 * np.random.randn(N)
            self.finished.emit(t, resp_volts, force_volts, fs)
        else:
            try:
                scope = RigolMHO934(self.ip)
                t, resp_volts, fs = scope.fetch_channel_data(1)
                _, force_volts, _ = scope.fetch_channel_data(4)
                scope.close()
                self.finished.emit(t, resp_volts, force_volts, fs)
            except Exception as e:
                self.error.emit(str(e))


# ==============================================================================
# 5. MAIN GUI APPLICATION
# ==============================================================================
class ModalAnalysisApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Rigol MHO934 Roving Modal Analysis Software (PyQt5)")
        self.resize(1400, 950)

        self.session = RovingModalSession(impact_node=1)
        self.current_runs_force = []
        self.current_runs_resp = []
        self.latest_fs = 10000.0

        self.init_ui()

    def init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # ----------------------------------------------------------------------
        # LEFT CONTROL PANEL
        # ----------------------------------------------------------------------
        left_panel = QVBoxLayout()
        main_layout.addLayout(left_panel, stretch=1)

        # Connection Box
        conn_box = QGroupBox("1. Scope Connection (MHO934)")
        conn_layout = QGridLayout(conn_box)
        conn_layout.addWidget(QLabel("IP Address:"), 0, 0)
        self.ip_edit = QLineEdit("10.0.0.66")
        conn_layout.addWidget(self.ip_edit, 0, 1)
        self.sim_check = QPushButton("Mode: Hardware")
        self.sim_check.setCheckable(True)
        self.sim_check.toggled.connect(self.toggle_sim_mode)
        conn_layout.addWidget(self.sim_check, 1, 0, 1, 2)
        left_panel.addWidget(conn_box)

        # Modal Configuration Box
        config_box = QGroupBox("2. Test Configuration & Calibration")
        config_layout = QGridLayout(config_box)
        
        config_layout.addWidget(QLabel("Fixed Impact Node (Ch4):"), 0, 0)
        self.impact_node_spin = QSpinBox()
        self.impact_node_spin.setValue(1)
        self.impact_node_spin.valueChanged.connect(self.update_impact_node)
        config_layout.addWidget(self.impact_node_spin, 0, 1)

        config_layout.addWidget(QLabel("Roving Resp Node (Ch1):"), 1, 0)
        self.resp_node_spin = QSpinBox()
        self.resp_node_spin.setRange(1, 9999)
        self.resp_node_spin.setValue(101)
        config_layout.addWidget(self.resp_node_spin, 1, 1)

        config_layout.addWidget(QLabel("Target Averages:"), 2, 0)
        self.target_avg_spin = QSpinBox()
        self.target_avg_spin.setValue(3)
        config_layout.addWidget(self.target_avg_spin, 2, 1)

        config_layout.addWidget(QLabel("Exp Decay Rate:"), 3, 0)
        self.exp_decay_spin = QDoubleSpinBox()
        self.exp_decay_spin.setSingleStep(0.01)
        self.exp_decay_spin.setValue(0.05)
        config_layout.addWidget(self.exp_decay_spin, 3, 1)

        # Sensor Sensitivity Adjustments
        config_layout.addWidget(QLabel("Accel Gain (mV/g):"), 4, 0)
        self.accel_gain_spin = QDoubleSpinBox()
        self.accel_gain_spin.setRange(0.001, 100000.0)
        self.accel_gain_spin.setValue(100.0)  # Default 100 mV/g
        self.accel_gain_spin.setDecimals(3)
        config_layout.addWidget(self.accel_gain_spin, 4, 1)

        config_layout.addWidget(QLabel("Hammer Gain (mV/N):"), 5, 0)
        self.hammer_gain_spin = QDoubleSpinBox()
        self.hammer_gain_spin.setRange(0.001, 100000.0)
        self.hammer_gain_spin.setValue(10.0)  # Default 10 mV/N
        self.hammer_gain_spin.setDecimals(3)
        config_layout.addWidget(self.hammer_gain_spin, 5, 1)

        left_panel.addWidget(config_box)

        # Acquisition Action Box
        acq_box = QGroupBox("3. Acquisition")
        acq_layout = QVBoxLayout(acq_box)
        
        self.btn_capture = QPushButton("Trigger / Capture Impact")
        self.btn_capture.setStyleSheet("font-weight: bold; background-color: #007ACC; color: white; height: 35px;")
        self.btn_capture.clicked.connect(self.capture_trigger)
        acq_layout.addWidget(self.btn_capture)

        self.btn_clear = QPushButton("Clear Current Point Averages")
        self.btn_clear.clicked.connect(self.clear_averages)
        acq_layout.addWidget(self.btn_clear)

        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        acq_layout.addWidget(self.progress_bar)
        left_panel.addWidget(acq_box)

        # Export Box
        export_box = QGroupBox("4. Save & Export")
        export_layout = QVBoxLayout(export_box)
        self.btn_next_point = QPushButton("Store Point & Advance Node")
        self.btn_next_point.setEnabled(False)
        self.btn_next_point.clicked.connect(self.store_point)
        export_layout.addWidget(self.btn_next_point)

        self.btn_export = QPushButton("Export Session to UNV 58")
        self.btn_export.setStyleSheet("font-weight: bold; background-color: #28A745; color: white; height: 30px;")
        self.btn_export.clicked.connect(self.export_unv)
        export_layout.addWidget(self.btn_export)
        left_panel.addWidget(export_box)

        # Log Window
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        left_panel.addWidget(self.log_text, stretch=1)

        # ----------------------------------------------------------------------
        # RIGHT PLOTTING DISPLAY (3x2 Subplot Layout)
        # ----------------------------------------------------------------------
        right_panel = QVBoxLayout()
        main_layout.addLayout(right_panel, stretch=3)

        self.fig = Figure(figsize=(10, 8))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        right_panel.addWidget(self.toolbar)
        right_panel.addWidget(self.canvas)

        # Setup 3x2 Grid Subplots
        self.ax_accel_time = self.fig.add_subplot(3, 2, 1)
        self.ax_hammer_time = self.fig.add_subplot(3, 2, 2)
        self.ax_accel_fft = self.fig.add_subplot(3, 2, 3)
        self.ax_hammer_fft = self.fig.add_subplot(3, 2, 4)
        self.ax_frf_mag = self.fig.add_subplot(3, 2, 5)
        self.ax_frf_phase = self.fig.add_subplot(3, 2, 6)

        self.all_axes = [
            self.ax_accel_time, self.ax_hammer_time,
            self.ax_accel_fft, self.ax_hammer_fft,
            self.ax_frf_mag, self.ax_frf_phase
        ]

        self.fig.tight_layout()
        self.log("Application initialized. Connect your scope or run in Simulation Mode.")

    # ----------------------------------------------------------------------
    # LOGIC & LOGGING
    # ----------------------------------------------------------------------
    def log(self, msg):
        self.log_text.append(f"> {msg}")

    def toggle_sim_mode(self, checked):
        if checked:
            self.sim_check.setText("Mode: SIMULATION (Mock Data)")
            self.sim_check.setStyleSheet("background-color: #D9534F; color: white;")
        else:
            self.sim_check.setText("Mode: Hardware")
            self.sim_check.setStyleSheet("")

    def update_impact_node(self, val):
        self.session.impact_node = val

    def capture_trigger(self):
        target_avg = self.target_avg_spin.value()
        if len(self.current_runs_force) >= target_avg:
            self.log("Target averages reached. Click 'Store Point & Advance Node' to save.")
            return

        self.btn_capture.setEnabled(False)
        self.log(f"Fetching trigger from scope (Ch1/Ch4)... [Avg {len(self.current_runs_force) + 1}/{target_avg}]")

        self.worker = FetchWorker(self.ip_edit.text(), self.sim_check.isChecked())
        self.worker.finished.connect(self.on_data_fetched)
        self.worker.error.connect(self.on_fetch_error)
        self.worker.start()

    def on_data_fetched(self, t, resp_volts, force_volts, fs):
        self.btn_capture.setEnabled(True)
        self.latest_fs = fs

        # Apply Sensor Gains (mV/g -> V/g and mV/N -> V/N)
        accel_sensitivity_v_per_g = self.accel_gain_spin.value() / 1000.0
        hammer_sensitivity_v_per_n = self.hammer_gain_spin.value() / 1000.0

        resp_g = resp_volts / accel_sensitivity_v_per_g
        force_n = force_volts / hammer_sensitivity_v_per_n

        # Apply Windows to Scaled Signals
        proc = ModalSignalProcessor(fs=fs, exp_decay_rate=self.exp_decay_spin.value())
        f_win, r_win = proc.apply_windows(force_n, resp_g)

        self.current_runs_force.append(f_win)
        self.current_runs_resp.append(r_win)

        target_avg = self.target_avg_spin.value()
        curr_avg = len(self.current_runs_force)
        self.progress_bar.setValue(int((curr_avg / target_avg) * 100))

        if curr_avg >= target_avg:
            self.btn_next_point.setEnabled(True)
            self.log("Averaging complete for this point.")

        self.plot_all(t, force_n, resp_g, f_win, r_win)

    def on_fetch_error(self, err_msg):
        self.btn_capture.setEnabled(True)
        self.log(f"ERROR: {err_msg}")
        QMessageBox.critical(self, "Connection Error", f"Failed to reach scope:\n{err_msg}")

    def clear_averages(self):
        self.current_runs_force.clear()
        self.current_runs_resp.clear()
        self.progress_bar.setValue(0)
        self.btn_next_point.setEnabled(False)
        for ax in self.all_axes:
            ax.cla()
        self.canvas.draw()
        self.log("Cleared current point buffers.")

    def plot_all(self, t, force_raw, resp_raw, f_win, r_win):
        # Clear all axes
        for ax in self.all_axes:
            ax.cla()

        N = len(force_raw)
        
        # ----------------------------------------------------------------------
        # 1. ACCELEROMETER TIME DOMAIN (Scaled to g)
        # ----------------------------------------------------------------------
        self.ax_accel_time.plot(t * 1000, resp_raw, 'b-', alpha=0.4, label='Raw')
        self.ax_accel_time.plot(t * 1000, r_win, 'b-', label='Windowed')
        self.ax_accel_time.set_title("Accelerometer Time Domain")
        self.ax_accel_time.set_xlabel("Time [ms]")
        self.ax_accel_time.set_ylabel("Acceleration [g]")
        self.ax_accel_time.legend(loc='upper right', fontsize='x-small')
        self.ax_accel_time.grid(True)

        # ----------------------------------------------------------------------
        # 2. IMPACT HAMMER TIME DOMAIN (Scaled to N)
        # ----------------------------------------------------------------------
        peak_idx = np.argmax(np.abs(force_raw))
        t_centered_ms = (t - t[peak_idx]) * 1000.0  # Align peak at t = 0 ms

        self.ax_hammer_time.plot(t_centered_ms, force_raw, 'r-', alpha=0.4, label='Raw')
        self.ax_hammer_time.plot(t_centered_ms, f_win, 'r-', label='Windowed')
        self.ax_hammer_time.set_xlim([-15, 20])  # Force -15ms to 20ms range
        self.ax_hammer_time.set_title("Impact Hammer Time Domain")
        self.ax_hammer_time.set_xlabel("Time relative to impact [ms]")
        self.ax_hammer_time.set_ylabel("Force [N]")
        self.ax_hammer_time.legend(loc='upper right', fontsize='x-small')
        self.ax_hammer_time.grid(True)

        # ----------------------------------------------------------------------
        # FREQUENCY DOMAIN COMPUTATIONS (15 Hz to 1.5 kHz Mask)
        # ----------------------------------------------------------------------
        freqs_raw = np.fft.rfftfreq(N, d=1/self.latest_fs)
        fft_accel = np.abs(np.fft.rfft(resp_raw))
        fft_hammer = np.abs(np.fft.rfft(force_raw))

        # Filter Mask: 15 Hz - 1500 Hz
        freq_mask_fft = (freqs_raw >= 15) & (freqs_raw <= 1500)

        # ----------------------------------------------------------------------
        # 3. RAW ACCELEROMETER FREQUENCY DOMAIN
        # ----------------------------------------------------------------------
        # self.ax_accel_fft.plot(freqs_raw[freq_mask_fft], fft_accel[freq_mask_fft], 'b-')
        self.ax_accel_fft.semilogy(freqs_raw[freq_mask_fft], fft_accel[freq_mask_fft], 'b-')
        self.ax_accel_fft.set_title("Accelerometer Spectrum (15 Hz - 1.5 kHz)")
        self.ax_accel_fft.set_xlabel("Frequency [Hz]")
        self.ax_accel_fft.set_ylabel("Magnitude [g]")
        self.ax_accel_fft.grid(True)

        # ----------------------------------------------------------------------
        # 4. RAW HAMMER FREQUENCY DOMAIN
        # ----------------------------------------------------------------------
        self.ax_hammer_fft.plot(freqs_raw[freq_mask_fft], fft_hammer[freq_mask_fft], 'r-')
        self.ax_hammer_fft.set_title("Hammer Spectrum (15 Hz - 1.5 kHz)")
        self.ax_hammer_fft.set_xlabel("Frequency [Hz]")
        self.ax_hammer_fft.set_ylabel("Magnitude [N]")
        self.ax_hammer_fft.grid(True)

        # ----------------------------------------------------------------------
        # FRF & PHASE COMPUTATIONS
        # ----------------------------------------------------------------------
        proc = ModalSignalProcessor(fs=self.latest_fs, exp_decay_rate=self.exp_decay_spin.value())
        freqs_frf, H1, _ = proc.compute_average_frf(self.current_runs_force, self.current_runs_resp)
        
        freq_mask_frf = (freqs_frf >= 15) & (freqs_frf <= 1500)
        mag_dB = 20 * np.log10(np.abs(H1) + 1e-12)
        phase_deg = np.angle(H1, deg=True)

        # ----------------------------------------------------------------------
        # 5. FRF MAGNITUDE PLOT
        # ----------------------------------------------------------------------
        self.ax_frf_mag.plot(freqs_frf[freq_mask_frf], mag_dB[freq_mask_frf], 'k-')
        self.ax_frf_mag.set_title(f"Averaged FRF H1 Magnitude (Avg: {len(self.current_runs_force)}/{self.target_avg_spin.value()})")
        self.ax_frf_mag.set_xlabel("Frequency [Hz]")
        self.ax_frf_mag.set_ylabel("Magnitude [dB re 1 g/N]")
        self.ax_frf_mag.grid(True)

        # ----------------------------------------------------------------------
        # 6. FRF PHASE PLOT
        # ----------------------------------------------------------------------
        self.ax_frf_phase.plot(freqs_frf[freq_mask_frf], phase_deg[freq_mask_frf], 'g-')
        self.ax_frf_phase.set_title("FRF H1 Phase Angle (15 Hz - 1.5 kHz)")
        self.ax_frf_phase.set_xlabel("Frequency [Hz]")
        self.ax_frf_phase.set_ylabel("Phase [deg]")
        self.ax_frf_phase.set_ylim([-180, 180])
        self.ax_frf_phase.grid(True)

        self.fig.tight_layout()
        self.canvas.draw()

    def store_point(self):
        proc = ModalSignalProcessor(fs=self.latest_fs, exp_decay_rate=self.exp_decay_spin.value())
        freqs, H1, coherence = proc.compute_average_frf(self.current_runs_force, self.current_runs_resp)

        node = self.resp_node_spin.value()
        self.session.add_point_data(resp_node=node, freq=freqs, H1_complex=H1, coherence=coherence)
        self.log(f"Stored Point #{node} in session dataset.")

        # Advance point
        self.resp_node_spin.setValue(node + 1)
        self.clear_averages()

    def export_unv(self):
        if not self.session.measurements:
            QMessageBox.warning(self, "Export Warning", "No measurement points saved yet.")
            return

        file_path, _ = QFileDialog.getSaveFileName(self, "Save UFF 58 File", "", "Universal File Format (*.unv)")
        if file_path:
            self.session.export_to_unv58(file_path)
            self.log(f"Successfully exported {len(self.session.measurements)} point(s) to: {file_path}")
            QMessageBox.information(self, "Export Successful", f"Saved UNV dataset to:\n{file_path}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ModalAnalysisApp()
    window.show()
    sys.exit(app.exec_())