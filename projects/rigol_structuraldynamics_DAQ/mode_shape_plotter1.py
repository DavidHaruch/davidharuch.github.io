import sys
import os
import re
import glob
import pandas as pd
import numpy as np

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QSlider, QLabel, QComboBox, QSplitter,
    QDoubleSpinBox, QLineEdit
)
from PyQt5.QtCore import Qt, QTimer

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D


class ODSVisualizer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ME'Scope Style ODS Visualizer")
        self.resize(1400, 800)

        # Data Containers
        self.nodes_df = None          # DataFrame: ['NODE', 'X', 'Y', 'Z']
        self.lines_list = []          # List of tuples: [(node_start, node_end), ...]
        self.frf_data = {}            # {node_id: dataframe}
        self.frequencies = None       # Complete frequency array
        self.filtered_indices = []    # Indices corresponding to min/max freq range
        self.selected_freq_idx = 0    # Index within self.frequencies

        # Animation state
        self.phase_deg = 0.0
        self.anim_gain = 15.0         # User visual scale factor (mm)
        self.peak_magnitude = 1.0     # Peak FRF magnitude at current frequency for normalization

        self.init_ui()

    def init_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)

        # --- Top Control Toolbar ---
        top_bar = QHBoxLayout()

        self.btn_load_csvs = QPushButton("Load FRF Directory")
        self.btn_load_csvs.clicked.connect(self.load_frf_directory)
        top_bar.addWidget(self.btn_load_csvs)

        self.btn_load_nodes = QPushButton("Load Nodes CSV")
        self.btn_load_nodes.clicked.connect(self.load_nodes_file)
        top_bar.addWidget(self.btn_load_nodes)

        self.btn_load_lines = QPushButton("Load Lines CSV")
        self.btn_load_lines.clicked.connect(self.load_lines_file)
        top_bar.addWidget(self.btn_load_lines)

        top_bar.addWidget(QLabel("Component:"))
        self.combo_component = QComboBox()
        self.combo_component.addItems(["Magnitude", "Imaginary"])
        self.combo_component.currentIndexChanged.connect(self.update_frf_plot)
        top_bar.addWidget(self.combo_component)

        self.btn_animate = QPushButton("Start/Stop Animation")
        self.btn_animate.setCheckable(True)
        self.btn_animate.toggled.connect(self.toggle_animation)
        top_bar.addWidget(self.btn_animate)

        main_layout.addLayout(top_bar)

        # --- Frequency & Connectivity Options ---
        opts_layout = QHBoxLayout()

        opts_layout.addWidget(QLabel("Min Freq (Hz):"))
        self.spin_min_freq = QDoubleSpinBox()
        self.spin_min_freq.setRange(0.0, 50000.0)
        self.spin_min_freq.setValue(5.0)
        self.spin_min_freq.setSingleStep(5.0)
        self.spin_min_freq.valueChanged.connect(self.apply_frequency_range)
        opts_layout.addWidget(self.spin_min_freq)

        opts_layout.addWidget(QLabel("Max Freq (Hz):"))
        self.spin_max_freq = QDoubleSpinBox()
        self.spin_max_freq.setRange(0.0, 50000.0)
        self.spin_max_freq.setValue(1500.0)
        self.spin_max_freq.setSingleStep(50.0)
        self.spin_max_freq.valueChanged.connect(self.apply_frequency_range)
        opts_layout.addWidget(self.spin_max_freq)

        opts_layout.addWidget(QLabel("Anim Gain (mm):"))
        self.spin_gain = QDoubleSpinBox()
        self.spin_gain.setRange(0.1, 500.0)
        self.spin_gain.setValue(15.0)
        self.spin_gain.valueChanged.connect(self.on_gain_changed)
        opts_layout.addWidget(self.spin_gain)

        opts_layout.addWidget(QLabel("Lines (e.g. 101-102,102-103):"))
        self.txt_lines = QLineEdit()
        self.txt_lines.setPlaceholderText("101-102, 102-103")
        self.txt_lines.editingFinished.connect(self.parse_manual_lines)
        opts_layout.addWidget(self.txt_lines)

        main_layout.addLayout(opts_layout)

        # --- Side-by-Side Plots (Horizontal Splitter) ---
        splitter = QSplitter(Qt.Horizontal)

        # Left: 2D FRF Overlay
        self.fig_frf = Figure(figsize=(6, 6))
        self.canvas_frf = FigureCanvas(self.fig_frf)
        self.ax_frf = self.fig_frf.add_subplot(111)
        splitter.addWidget(self.canvas_frf)

        # Right: 3D Animated Structure
        self.fig_ods = Figure(figsize=(6, 6))
        self.canvas_ods = FigureCanvas(self.fig_ods)
        self.ax_ods = self.fig_ods.add_subplot(111, projection='3d')
        splitter.addWidget(self.canvas_ods)

        splitter.setSizes([600, 600])
        main_layout.addWidget(splitter)

        # --- Bottom Slider Bar ---
        bottom_bar = QHBoxLayout()
        bottom_bar.addWidget(QLabel("Frequency Selector:"))
        self.slider_freq = QSlider(Qt.Horizontal)
        self.slider_freq.setMinimum(0)
        self.slider_freq.setMaximum(100)
        self.slider_freq.valueChanged.connect(self.on_slider_change)
        bottom_bar.addWidget(self.slider_freq)

        self.lbl_freq_val = QLabel("Freq: -- Hz")
        bottom_bar.addWidget(self.lbl_freq_val)

        main_layout.addLayout(bottom_bar)

        self.canvas_frf.mpl_connect('button_press_event', self.on_frf_click)

        # Timer for ODS Animation (~30 FPS)
        self.timer = QTimer()
        self.timer.setInterval(33)
        self.timer.timeout.connect(self.update_ods_frame)

    # --- Data Loading Handlers ---
    def load_frf_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select FRF CSV Directory")
        if not dir_path:
            return

        csv_files = glob.glob(os.path.join(dir_path, "*.csv"))
        self.frf_data.clear()

        for file_path in csv_files:
            fname = os.path.basename(file_path)
            match = re.search(r'Accel_Node(\d+)', fname, re.IGNORECASE)
            if match:
                node_id = int(match.group(1))
            else:
                continue

            df = pd.read_csv(file_path, comment='#')
            df.columns = [c.strip() for c in df.columns]
            self.frf_data[node_id] = df

        if self.frf_data:
            sample_node = list(self.frf_data.keys())[0]
            self.frequencies = self.frf_data[sample_node].iloc[:, 0].values

            if self.nodes_df is None:
                self.generate_default_grid()

            self.apply_frequency_range()

    def load_nodes_file(self):
        """ CSV format: Node, X, Y, Z """
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Node Coordinates CSV", "", "CSV Files (*.csv)")
        if not file_path:
            return

        df = pd.read_csv(file_path)
        df.columns = [c.strip().upper() for c in df.columns]
        self.nodes_df = df
        self.update_ods_plot()

    def load_lines_file(self):
        """ CSV format: NodeA, NodeB """
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Line Connections CSV", "", "CSV Files (*.csv)")
        if not file_path:
            return

        df = pd.read_csv(file_path)
        df.columns = [c.strip().upper() for c in df.columns]
        self.lines_list = [(row.iloc[0], row.iloc[1]) for _, row in df.iterrows()]
        
        line_strs = [f"{int(a)}-{int(b)}" for a, b in self.lines_list]
        self.txt_lines.setText(", ".join(line_strs))
        self.update_ods_plot()

    def parse_manual_lines(self):
        text = self.txt_lines.text().strip()
        if not text:
            self.lines_list = []
        else:
            new_lines = []
            segments = text.split(',')
            for seg in segments:
                parts = seg.strip().split('-')
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    new_lines.append((int(parts[0]), int(parts[1])))
            self.lines_list = new_lines
        self.update_ods_plot()

    def generate_default_grid(self):
        nodes = sorted(list(self.frf_data.keys()))
        coords = [{'NODE': n, 'X': i * 50.0, 'Y': 0.0, 'Z': 0.0} for i, n in enumerate(nodes)]
        self.nodes_df = pd.DataFrame(coords)
        self.lines_list = [(nodes[i], nodes[i+1]) for i in range(len(nodes)-1)]

    # --- Frequency & View Logic ---
    def apply_frequency_range(self):
        if self.frequencies is None:
            return

        f_min = self.spin_min_freq.value()
        f_max = self.spin_max_freq.value()

        self.filtered_indices = np.where((self.frequencies >= f_min) & (self.frequencies <= f_max))[0]

        if len(self.filtered_indices) == 0:
            self.filtered_indices = [0]

        self.slider_freq.blockSignals(True)
        self.slider_freq.setMinimum(0)
        self.slider_freq.setMaximum(len(self.filtered_indices) - 1)
        self.slider_freq.setValue(0)
        self.slider_freq.blockSignals(False)

        self.selected_freq_idx = self.filtered_indices[0]
        self.update_peak_normalization()
        self.update_frf_plot()
        self.update_ods_plot()

    def on_slider_change(self, val):
        if len(self.filtered_indices) > val:
            self.selected_freq_idx = self.filtered_indices[val]
            freq_val = self.frequencies[self.selected_freq_idx]
            self.lbl_freq_val.setText(f"Freq: {freq_val:.2f} Hz")
            self.update_peak_normalization()
            self.update_frf_plot()
            self.update_ods_plot()

    def on_gain_changed(self, val):
        self.anim_gain = val
        self.update_ods_plot()

    def update_peak_normalization(self):
        mags = []
        for df in self.frf_data.values():
            real_col = [c for c in df.columns if 'Real' in c][0]
            imag_col = [c for c in df.columns if 'Imag' in c][0]
            r = df[real_col].iloc[self.selected_freq_idx]
            i = df[imag_col].iloc[self.selected_freq_idx]
            mags.append(np.sqrt(r**2 + i**2))
        
        max_m = np.max(mags) if len(mags) > 0 else 1.0
        self.peak_magnitude = max_m if max_m > 1e-12 else 1.0

    # --- Plotting Methods ---
    def update_frf_plot(self):
        if not self.frf_data or self.frequencies is None:
            return

        self.ax_frf.clear()
        comp_type = self.combo_component.currentText()
        f_min, f_max = self.spin_min_freq.value(), self.spin_max_freq.value()

        for node_id, df in self.frf_data.items():
            freq = df.iloc[:, 0]
            if comp_type == "Magnitude":
                col = [c for c in df.columns if 'Mag' in c][0]
            else:
                col = [c for c in df.columns if 'Imag' in c][0]

            self.ax_frf.plot(freq, df[col], alpha=0.6, label=f'Node {node_id}')

        current_freq = self.frequencies[self.selected_freq_idx]
        self.ax_frf.axvline(x=current_freq, color='red', linestyle='--', linewidth=2)

        self.ax_frf.set_xlim(f_min, f_max)
        self.ax_frf.set_title(f"FRF Overlay ({comp_type})")
        self.ax_frf.set_xlabel("Frequency (Hz)")
        self.ax_frf.set_ylabel(comp_type)
        self.ax_frf.grid(True)
        self.canvas_frf.draw()

    def on_frf_click(self, event):
        if event.inaxes == self.ax_frf and self.frequencies is not None:
            clicked_freq = event.xdata
            subset_freqs = self.frequencies[self.filtered_indices]
            sub_idx = (np.abs(subset_freqs - clicked_freq)).argmin()
            self.slider_freq.setValue(sub_idx)

    def update_ods_plot(self):
        if self.nodes_df is None or not self.frf_data:
            return

        elev, azim = self.ax_ods.elev, self.ax_ods.azim
        self.ax_ods.clear()

        # Build node coordinate dictionary
        node_coords = {}
        for _, row in self.nodes_df.iterrows():
            node_coords[int(row['NODE'])] = (row['X'], row['Y'], row['Z'])

        node_ids = np.array(list(node_coords.keys()))
        x_base = np.array([node_coords[n][0] for n in node_ids])
        y_base = np.array([node_coords[n][1] for n in node_ids])
        z_base = np.array([node_coords[n][2] for n in node_ids])

        # Compute complex response displacement
        z_disp = np.zeros(len(node_ids))
        phase_rad = np.radians(self.phase_deg)

        for i, nid in enumerate(node_ids):
            if nid in self.frf_data:
                df = self.frf_data[nid]
                real_col = [c for c in df.columns if 'Real' in c][0]
                imag_col = [c for c in df.columns if 'Imag' in c][0]

                r = df[real_col].iloc[self.selected_freq_idx]
                im = df[imag_col].iloc[self.selected_freq_idx]

                complex_val = complex(r, im)
                z_disp[i] = (complex_val * np.exp(1j * phase_rad)).real

        scaled_z_disp = (z_disp / self.peak_magnitude) * self.anim_gain
        z_deformed = z_base + scaled_z_disp

        deformed_coords = {nid: (x, y, zd) for nid, x, y, zd in zip(node_ids, x_base, y_base, z_deformed)}

        # Plot Undeformed Wireframe Lines
        for na, nb in self.lines_list:
            if na in node_coords and nb in node_coords:
                xa, ya, za = node_coords[na]
                xb, yb, zb = node_coords[nb]
                self.ax_ods.plot([xa, xb], [ya, yb], [za, zb], 'k--', alpha=0.3)

        # Plot Deformed Wireframe Lines
        for na, nb in self.lines_list:
            if na in deformed_coords and nb in deformed_coords:
                xa, ya, za = deformed_coords[na]
                xb, yb, zb = deformed_coords[nb]
                self.ax_ods.plot([xa, xb], [ya, yb], [za, zb], 'r-', linewidth=2)

        # Plot Nodes & Labels
        self.ax_ods.scatter(x_base, y_base, z_base, color='gray', s=15, alpha=0.3)
        self.ax_ods.scatter(x_base, y_base, z_deformed, color='red', s=25)

        for nid, x, y, zd in zip(node_ids, x_base, y_base, z_deformed):
            self.ax_ods.text(x, y, zd, f' {nid}', color='blue', fontsize=8)

        # --- Enforce True 1:1 Spatial Aspect Ratio & Fixed Bounds ---
        max_span = max(
            np.ptp(x_base),
            np.ptp(y_base),
            max(self.anim_gain * 2.0, np.ptp(z_base))
        )
        if max_span == 0:
            max_span = 1.0

        # Center points in 3D space
        x_mid = (np.max(x_base) + np.min(x_base)) / 2.0
        y_mid = (np.max(y_base) + np.min(y_base)) / 2.0
        z_mid = (np.max(z_base) + np.min(z_base)) / 2.0

        half_span = max_span * 0.55  # 10% margin padding around structure

        self.ax_ods.set_xlim(x_mid - half_span, x_mid + half_span)
        self.ax_ods.set_ylim(y_mid - half_span, y_mid + half_span)
        self.ax_ods.set_zlim(z_mid - half_span, z_mid + half_span)

        # Lock physical 1:1 ratio across X, Y, Z axes
        self.ax_ods.set_box_aspect((1, 1, 1))

        self.ax_ods.set_title(f"ODS @ {self.frequencies[self.selected_freq_idx]:.2f} Hz")
        self.ax_ods.set_xlabel("X (mm)")
        self.ax_ods.set_ylabel("Y (mm)")
        self.ax_ods.set_zlabel("Z Deflection (mm)")
        self.ax_ods.view_init(elev=elev, azim=azim)

        self.canvas_ods.draw()

    # --- Animation Control ---
    def toggle_animation(self, checked):
        if checked:
            self.timer.start()
        else:
            self.timer.stop()

    def update_ods_frame(self):
        self.phase_deg = (self.phase_deg + 12) % 360
        self.update_ods_plot()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ODSVisualizer()
    window.show()
    sys.exit(app.exec_())