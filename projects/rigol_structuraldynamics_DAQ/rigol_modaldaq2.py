import sys
import csv
import numpy as np
import scipy.signal as signal
import pyvisa

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QDoubleSpinBox,
    QGroupBox, QTextEdit, QProgressBar, QFileDialog, QMessageBox, QCheckBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


# ==============================================================================
# 1. HARDWARE COMMUNICATION ENGINE (RIGOL MHO934 16-BIT LOW-NOISE FRONTEND)
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
        
        # Enforce 16-bit WORD transmission mode for Hi-Res Scope ADC
        self.scope.write(":WAVeform:FORMat WORD")
        
        preamble = self.scope.query(":WAVeform:PREamble?").split(',')
        x_incr = float(preamble[4])
        x_orig = float(preamble[5])
        y_incr = float(preamble[7])
        y_orig = float(preamble[8])
        y_ref = float(preamble[9])
        
        # Fetch 16-bit unsigned integers ('H' datatype)
        raw_data = self.scope.query_binary_values(":WAVeform:DATA?", datatype='H', container=np.array)
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

    def process_signals(self, force, resp, apply_hammer_win=True, apply_accel_win=True, remove_hammer_dc=True):
        N = len(force)
        t = np.arange(N) / self.fs
        
        # Feature 4: Hammer DC Offset Removal
        if remove_hammer_dc:
            # Estimate pre-trigger baseline noise floor from first 1% of time series
            pre_trigger_samples = max(10, int(0.01 * N))
            dc_offset = np.mean(force[:pre_trigger_samples])
            force_proc = force - dc_offset
        else:
            force_proc = force.copy()
            
        resp_proc = resp.copy()

        # Feature 3: Selective Windowing Options
        if apply_hammer_win:
            peak_idx = np.argmax(np.abs(force_proc))
            force_win = np.zeros(N)
            win_len = int(0.02 * self.fs)
            start = max(0, peak_idx - int(0.002 * self.fs))
            end = min(N, start + win_len)
            force_win[start:end] = 1.0
            force_proc = force_proc * force_win

        if apply_accel_win:
            resp_win = np.exp(-self.exp_decay_rate * t)
            resp_proc = resp_proc * resp_win

        return force_proc, resp_proc

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
        coherence = (np.abs(Sxf_avg)**2) / (Sff_avg * Sxx_avg + 1e-15)
        return f, H1, coherence


# ==============================================================================
# 3. VISA WORKER THREAD
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
            fs = 10000.0
            N = 8192
            t = np.arange(N) / fs
            # Synthetic 16-bit equivalent time domain signals
            force_volts = 0.5 * np.exp(-((t - 0.005) / 0.0008)**2) + 0.005 + 0.001 * np.random.randn(N)
            resp_volts = 0.08 * np.exp(-12 * t) * np.sin(2 * np.pi * 180 * t) + \
                         0.04 * np.exp(-25 * t) * np.sin(2 * np.pi * 420 * t) + \
                         0.0005 * np.random.randn(N)
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
# 4. MAIN GUI APPLICATION
# ==============================================================================
class ModalAnalysisApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Rigol MHO934 Roving Modal Analysis Software")
        self.resize(1500, 1000)

        self.current_runs_force = []
        self.current_runs_resp = []
        self.latest_t = None
        self.latest_force_raw = None
        self.latest_resp_raw = None
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

        # 1. Connection Box
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

        # 2. Configuration & Signal Conditioning Box
        config_box = QGroupBox("2. Test Config & Signal Conditioning")
        config_layout = QGridLayout(config_box)

        config_layout.addWidget(QLabel("Fixed Impact Node (Ch4):"), 0, 0)
        self.impact_node_spin = QSpinBox()
        self.impact_node_spin.setValue(1)
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

        # Signal Processing Toggles & Windows
        self.chk_hammer_dc = QCheckBox("Hammer DC Offset Removal")
        self.chk_hammer_dc.setChecked(True)
        config_layout.addWidget(self.chk_hammer_dc, 3, 0, 1, 2)

        self.chk_hammer_win = QCheckBox("Apply Force Window (Hammer)")
        self.chk_hammer_win.setChecked(True)
        config_layout.addWidget(self.chk_hammer_win, 4, 0, 1, 2)

        self.chk_accel_win = QCheckBox("Apply Exp Window (Accel)")
        self.chk_accel_win.setChecked(True)
        config_layout.addWidget(self.chk_accel_win, 5, 0, 1, 2)

        config_layout.addWidget(QLabel("Exp Decay Rate:"), 6, 0)
        self.exp_decay_spin = QDoubleSpinBox()
        self.exp_decay_spin.setSingleStep(0.01)
        self.exp_decay_spin.setValue(0.05)
        config_layout.addWidget(self.exp_decay_spin, 6, 1)

        config_layout.addWidget(QLabel("Accel Gain (mV/g):"), 7, 0)
        self.accel_gain_spin = QDoubleSpinBox()
        self.accel_gain_spin.setRange(0.001, 100000.0)
        self.accel_gain_spin.setValue(100.0)
        self.accel_gain_spin.setDecimals(3)
        config_layout.addWidget(self.accel_gain_spin, 7, 1)

        config_layout.addWidget(QLabel("Hammer Gain (mV/N):"), 8, 0)
        self.hammer_gain_spin = QDoubleSpinBox()
        self.hammer_gain_spin.setRange(0.001, 100000.0)
        self.hammer_gain_spin.setValue(10.0)
        self.hammer_gain_spin.setDecimals(3)
        config_layout.addWidget(self.hammer_gain_spin, 8, 1)

        left_panel.addWidget(config_box)

        # 3. Acquisition Box
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

        # 4. CSV Export Box
        export_box = QGroupBox("4. Export CSV Data")
        export_layout = QVBoxLayout(export_box)

        self.btn_export_time_csv = QPushButton("Export Time Domain Data (CSV)")
        self.btn_export_time_csv.clicked.connect(self.export_time_csv)
        export_layout.addWidget(self.btn_export_time_csv)

        self.btn_export_freq_csv = QPushButton("Export Frequency Domain / FRF (CSV)")
        self.btn_export_freq_csv.setStyleSheet("font-weight: bold; background-color: #28A745; color: white; height: 32px;")
        self.btn_export_freq_csv.clicked.connect(self.export_freq_csv)
        export_layout.addWidget(self.btn_export_freq_csv)

        left_panel.addWidget(export_box)

        # Log Window
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        left_panel.addWidget(self.log_text, stretch=1)

        # ----------------------------------------------------------------------
        # RIGHT PLOTTING DISPLAY (4x2 Subplot Layout)
        # ----------------------------------------------------------------------
        right_panel = QVBoxLayout()
        main_layout.addLayout(right_panel, stretch=3)

        self.fig = Figure(figsize=(10, 10))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        right_panel.addWidget(self.toolbar)
        right_panel.addWidget(self.canvas)

        # Setup 4x2 Grid Subplots
        self.ax_accel_time = self.fig.add_subplot(4, 2, 1)
        self.ax_hammer_time = self.fig.add_subplot(4, 2, 2)
        self.ax_accel_fft = self.fig.add_subplot(4, 2, 3)
        self.ax_hammer_fft = self.fig.add_subplot(4, 2, 4)
        self.ax_frf_mag = self.fig.add_subplot(4, 2, 5)
        self.ax_frf_phase = self.fig.add_subplot(4, 2, 6)
        self.ax_coherence = self.fig.add_subplot(4, 2, 7)
        self.ax_blank = self.fig.add_subplot(4, 2, 8)
        self.ax_blank.axis('off')

        self.all_axes = [
            self.ax_accel_time, self.ax_hammer_time,
            self.ax_accel_fft, self.ax_hammer_fft,
            self.ax_frf_mag, self.ax_frf_phase, self.ax_coherence
        ]

        self.fig.tight_layout()
        self.log("Application initialized. High-Res 16-bit Mode Enabled.")

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

    def capture_trigger(self):
        target_avg = self.target_avg_spin.value()
        if len(self.current_runs_force) >= target_avg:
            self.log("Target averages reached for this point.")
            return

        self.btn_capture.setEnabled(False)
        self.log(f"Fetching 16-bit raw waveform... [Avg {len(self.current_runs_force) + 1}/{target_avg}]")

        self.worker = FetchWorker(self.ip_edit.text(), self.sim_check.isChecked())
        self.worker.finished.connect(self.on_data_fetched)
        self.worker.error.connect(self.on_fetch_error)
        self.worker.start()

    def on_data_fetched(self, t, resp_volts, force_volts, fs):
        self.btn_capture.setEnabled(True)
        self.latest_fs = fs
        self.latest_t = t

        accel_sensitivity_v_per_g = self.accel_gain_spin.value() / 1000.0
        hammer_sensitivity_v_per_n = self.hammer_gain_spin.value() / 1000.0

        self.latest_resp_raw = resp_volts / accel_sensitivity_v_per_g
        self.latest_force_raw = force_volts / hammer_sensitivity_v_per_n

        proc = ModalSignalProcessor(fs=fs, exp_decay_rate=self.exp_decay_spin.value())
        
        f_proc, r_proc = proc.process_signals(
            force=self.latest_force_raw,
            resp=self.latest_resp_raw,
            apply_hammer_win=self.chk_hammer_win.isChecked(),
            apply_accel_win=self.chk_accel_win.isChecked(),
            remove_hammer_dc=self.chk_hammer_dc.isChecked()
        )

        self.current_runs_force.append(f_proc)
        self.current_runs_resp.append(r_proc)

        target_avg = self.target_avg_spin.value()
        curr_avg = len(self.current_runs_force)
        self.progress_bar.setValue(int((curr_avg / target_avg) * 100))

        self.plot_all(t, self.latest_force_raw, self.latest_resp_raw, f_proc, r_proc)

    def on_fetch_error(self, err_msg):
        self.btn_capture.setEnabled(True)
        self.log(f"ERROR: {err_msg}")
        QMessageBox.critical(self, "Connection Error", f"Failed to reach scope:\n{err_msg}")

    def clear_averages(self):
        self.current_runs_force.clear()
        self.current_runs_resp.clear()
        self.progress_bar.setValue(0)
        for ax in self.all_axes:
            ax.cla()
        self.canvas.draw()
        self.log("Cleared current point buffers.")

    def plot_all(self, t, force_raw, resp_raw, f_proc, r_proc):
        for ax in self.all_axes:
            ax.cla()

        N = len(force_raw)

        # 1. Accelerometer Time Domain
        self.ax_accel_time.plot(t * 1000, resp_raw, 'b-', alpha=0.35, label='Raw')
        self.ax_accel_time.plot(t * 1000, r_proc, 'b-', label='Processed')
        self.ax_accel_time.set_title("Accelerometer Time Domain")
        self.ax_accel_time.set_xlabel("Time [ms]")
        self.ax_accel_time.set_ylabel("Accel [g]")
        self.ax_accel_time.legend(loc='upper right', fontsize='x-small')
        self.ax_accel_time.grid(True)

        # 2. Impact Hammer Time Domain
        peak_idx = np.argmax(np.abs(force_raw))
        t_centered_ms = (t - t[peak_idx]) * 1000.0

        self.ax_hammer_time.plot(t_centered_ms, force_raw, 'r-', alpha=0.35, label='Raw')
        self.ax_hammer_time.plot(t_centered_ms, f_proc, 'r-', label='Processed')
        self.ax_hammer_time.set_xlim([-15, 20])
        self.ax_hammer_time.set_title("Impact Hammer Time Domain")
        self.ax_hammer_time.set_xlabel("Time rel. to impact [ms]")
        self.ax_hammer_time.set_ylabel("Force [N]")
        self.ax_hammer_time.legend(loc='upper right', fontsize='x-small')
        self.ax_hammer_time.grid(True)

        # Frequency Calculations
        freqs_raw = np.fft.rfftfreq(N, d=1/self.latest_fs)
        fft_accel = np.abs(np.fft.rfft(resp_raw))
        fft_hammer = np.abs(np.fft.rfft(force_raw))

        freq_mask_fft = (freqs_raw >= 15) & (freqs_raw <= 1500)

        # 3. Accel Frequency Spectrum
        self.ax_accel_fft.semilogy(freqs_raw[freq_mask_fft], fft_accel[freq_mask_fft], 'b-')
        self.ax_accel_fft.set_title("Accelerometer Spectrum (15 Hz - 1.5 kHz)")
        self.ax_accel_fft.set_xlabel("Frequency [Hz]")
        self.ax_accel_fft.set_ylabel("Mag [g]")
        self.ax_accel_fft.grid(True)

        # 4. Hammer Frequency Spectrum
        self.ax_hammer_fft.plot(freqs_raw[freq_mask_fft], fft_hammer[freq_mask_fft], 'r-')
        self.ax_hammer_fft.set_title("Hammer Spectrum (15 Hz - 1.5 kHz)")
        self.ax_hammer_fft.set_xlabel("Frequency [Hz]")
        self.ax_hammer_fft.set_ylabel("Mag [N]")
        self.ax_hammer_fft.grid(True)

        # FRF & Coherence Computation
        proc = ModalSignalProcessor(fs=self.latest_fs, exp_decay_rate=self.exp_decay_spin.value())
        freqs_frf, H1, coherence = proc.compute_average_frf(self.current_runs_force, self.current_runs_resp)
        
        # Feature 5: Note Frequency Step Size (delta_f) on Plot
        delta_f = freqs_frf[1] - freqs_frf[0] if len(freqs_frf) > 1 else 0.0

        freq_mask_frf = (freqs_frf >= 15) & (freqs_frf <= 1500)
        mag_dB = 20 * np.log10(np.abs(H1) + 1e-12)
        phase_deg = np.angle(H1, deg=True)

        # 5. FRF Magnitude Plot
        self.ax_frf_mag.plot(freqs_frf[freq_mask_frf], mag_dB[freq_mask_frf], 'k-')
        self.ax_frf_mag.set_title(f"FRF Magnitude | Avg: {len(self.current_runs_force)}/{self.target_avg_spin.value()} | df = {delta_f:.3f} Hz")
        self.ax_frf_mag.set_xlabel("Frequency [Hz]")
        self.ax_frf_mag.set_ylabel("Mag [dB re 1 g/N]")
        self.ax_frf_mag.grid(True)

        # 6. FRF Phase Plot
        self.ax_frf_phase.plot(freqs_frf[freq_mask_frf], phase_deg[freq_mask_frf], 'g-')
        self.ax_frf_phase.set_title(f"FRF Phase Angle | df = {delta_f:.3f} Hz")
        self.ax_frf_phase.set_xlabel("Frequency [Hz]")
        self.ax_frf_phase.set_ylabel("Phase [deg]")
        self.ax_frf_phase.set_ylim([-180, 180])
        self.ax_frf_phase.grid(True)

        # 7. Coherence Plot
        self.ax_coherence.plot(freqs_frf[freq_mask_frf], coherence[freq_mask_frf], 'm-')
        self.ax_coherence.set_title(f"Coherence gamma^2 | df = {delta_f:.3f} Hz")
        self.ax_coherence.set_xlabel("Frequency [Hz]")
        self.ax_coherence.set_ylabel("Coherence [0 - 1]")
        self.ax_coherence.set_ylim([0, 1.05])
        self.ax_coherence.grid(True)

        self.fig.tight_layout()
        self.canvas.draw()

    # ----------------------------------------------------------------------
    # FEATURE 2: CSV EXPORTERS
    # ----------------------------------------------------------------------
    def export_freq_csv(self):
        if len(self.current_runs_force) == 0:
            QMessageBox.warning(self, "Export Warning", "No acquired runs available to compute FRF for export.")
            return

        impact_node = self.impact_node_spin.value()
        accel_node = self.resp_node_spin.value()
        default_filename = f"FRF_H1_Hammer_Node{impact_node}_Accel_Node{accel_node}.csv"

        file_path, _ = QFileDialog.getSaveFileName(self, "Save FRF Data CSV", default_filename, "CSV Files (*.csv)")
        if file_path:
            proc = ModalSignalProcessor(fs=self.latest_fs, exp_decay_rate=self.exp_decay_spin.value())
            freqs, H1, coherence = proc.compute_average_frf(self.current_runs_force, self.current_runs_resp)
            
            # Restrict frequency array to 5 kHz maximum
            mask_5k = (freqs >= 0) & (freqs <= 5000.0)
            freqs_export = freqs[mask_5k]
            H1_export = H1[mask_5k]
            coherence_export = coherence[mask_5k]

            mag_dB = 20 * np.log10(np.abs(H1_export) + 1e-12)
            phase_deg = np.angle(H1_export, deg=True)

            with open(file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                
                # Metadata Header Block (Explicit Node IDs and Test Context)
                writer.writerow(["# MODAL TEST FRF EXPORT"])
                writer.writerow(["# Hammer Impact Node (Force Input)", impact_node])
                writer.writerow(["# Accelerometer Node (Response Output)", accel_node])
                writer.writerow(["# Number of Averages", len(self.current_runs_force)])
                writer.writerow(["# Sample Rate (Hz)", self.latest_fs])
                writer.writerow([])
                
                # Data Header
                writer.writerow(["Frequency_Hz", "FRF_Real_g_N", "FRF_Imag_g_N", "FRF_Mag_dB_re_1g_N", "FRF_Phase_deg", "Coherence"])
                
                for i in range(len(freqs_export)):
                    writer.writerow([
                        f"{freqs_export[i]:.4f}",
                        f"{np.real(H1_export[i]):.6e}",
                        f"{np.imag(H1_export[i]):.6e}",
                        f"{mag_dB[i]:.4f}",
                        f"{phase_deg[i]:.4f}",
                        f"{coherence_export[i]:.6f}"
                    ])

            self.log(f"Exported FRF (H1 Node {accel_node}/{impact_node}) up to 5 kHz to: {file_path}")

            # Automatically advance Accelerometer Node Number and Reset Buffer
            next_node = accel_node + 1
            self.resp_node_spin.setValue(next_node)
            self.log(f"Advanced Accelerometer Node to #{next_node}.")
            self.clear_averages()

    def export_time_csv(self):
        if self.latest_t is None or len(self.current_runs_force) == 0:
            QMessageBox.warning(self, "Export Warning", "No time domain data captured yet.")
            return

        impact_node = self.impact_node_spin.value()
        accel_node = self.resp_node_spin.value()
        default_filename = f"Time_Domain_Hammer_Node{impact_node}_Accel_Node{accel_node}.csv"

        file_path, _ = QFileDialog.getSaveFileName(self, "Save Time Domain CSV", default_filename, "CSV Files (*.csv)")
        if file_path:
            with open(file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["# MODAL TIME DOMAIN WAVEFORM EXPORT"])
                writer.writerow(["# Hammer Impact Node", impact_node])
                writer.writerow(["# Accelerometer Node", accel_node])
                writer.writerow([])
                writer.writerow(["Time_s", "Raw_Hammer_N", "Processed_Hammer_N", "Raw_Accel_g", "Processed_Accel_g"])
                for i in range(len(self.latest_t)):
                    writer.writerow([
                        f"{self.latest_t[i]:.6e}",
                        f"{self.latest_force_raw[i]:.6e}",
                        f"{self.current_runs_force[-1][i]:.6e}",
                        f"{self.latest_resp_raw[i]:.6e}",
                        f"{self.current_runs_resp[-1][i]:.6e}"
                    ])
            self.log(f"Exported raw time domain data to: {file_path}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ModalAnalysisApp()
    window.show()
    sys.exit(app.exec_())