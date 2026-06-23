import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QGridLayout, QLabel, QSlider)

class AmplifierSimGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LM3886 VCM Precision Current Driver Loop-Shaper")
        self.setGeometry(100, 100, 1400, 850)
        
        # Initialize default baseline values (Updated to your real 4.31mH coil!)
        self.params = {
            'R1': 0.55,     # Coil R
            'L1': 4.31,     # Measured Coil L in mH
            'R2': 0.1,      # Shunt R
            'R10': 670.0,   # Feedback R
            'R11': 1210.0,  # Input R
            'R12': 67.0,    # Comp R in kOhm
            'C1': 22.0,     # Comp C1 in nF
            'C2': 100.0,    # Comp C2 in pF
            'GBW_U3': 10.0,  # Error Amp GBW in MHz
            'GBW_XU2': 8.0, # LM3886 GBW in MHz
        }
        
        # Setup frequency range vector once (High-density 800 points is plenty for graphics)
        self.f = np.logspace(2, 5.5, 50)
        self.w = 2 * np.pi * self.f
        self.s = 1j * self.w

        self.init_ui()
        self.update_plots()

    def init_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        # Left Column: Interactive Control Sliders Panel
        controls_widget = QWidget()
        controls_widget.setFixedWidth(420)
        controls_layout = QVBoxLayout(controls_widget)
        
        title_label = QLabel("Component Tuning Panel")
        title_label.setStyleSheet("font-weight: bold; font-size: 14px; color: #2c3e50;")
        controls_layout.addWidget(title_label)
        
        grid = QGridLayout()
        grid.setSpacing(8)
        
        self.sliders = {}
        self.val_labels = {}
        
        # Config array: (key, Display Name, min_val, max_val, steps, scale_factor, unit)
        slider_config = [
            ('R11', 'R11 (Input R)', 100, 5000, 4900, 1.0, "Ω"),
            ('R10', 'R10 (Feedback R)', 100, 5000, 4900, 1.0, "Ω"),
            ('R12', 'R12 (Comp R)', 1, 200, 199, 1.0, "kΩ"),
            ('C1', 'C1 (Comp Cap)', 1, 200, 199, 1.0, "nF"),
            ('C2', 'C2 (HFX Cap)', 10, 1000, 990, 1.0, "pF"),
            ('L1', 'L1 (Coil Inductance)', 1, 200, 199, 0.1, "mH"),
            ('R1', 'R1 (Coil Resistance)', 1, 50, 49, 0.1, "Ω"),
            ('GBW_U3', 'GBW Error Amp', 1, 20, 19, 1.0, "MHz"),
            ('GBW_XU2', 'GBW LM3886', 1, 30, 29, 1.0, "MHz")
        ]

        for idx, (key, name, p_min, p_max, steps, scale, unit) in enumerate(slider_config):
            lbl = QLabel(name)
            v_lbl = QLabel()
            v_lbl.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            v_lbl.setStyleSheet("font-family: monospace; font-weight: bold; color: #16a085;")
            
            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(0)
            slider.setMaximum(steps)
            
            # Map baseline values precisely to step scales
            if key in ['L1', 'R1']:
                init_step = int((self.params[key] * 10) - p_min)
            else:
                init_step = int(self.params[key] - p_min)
            slider.setValue(init_step)
            slider.valueChanged.connect(self.update_plots)
            
            grid.addWidget(lbl, idx*2, 0)
            grid.addWidget(v_lbl, idx*2, 1)
            grid.addWidget(slider, idx*2+1, 0, 1, 2)
            
            self.sliders[key] = (slider, p_min, scale, key)
            self.val_labels[key] = (v_lbl, unit)

        controls_layout.addLayout(grid)
        controls_layout.addStretch()
        
        self.status_box = QLabel()
        self.status_box.setStyleSheet("border: 1px solid #bdc3c7; background-color: #fdfefe; "
                                     "font-family: monospace; font-size: 11px; padding: 10px;")
        controls_layout.addWidget(self.status_box)
        main_layout.addWidget(controls_widget)

        # Right Column: Matplotlib Subplots Canvas
        self.fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = plt.subplots(2, 2, figsize=(10, 8))
        self.canvas = FigureCanvas(self.fig)
        main_layout.addWidget(self.canvas, stretch=1)

        # --- CRITICAL PERFORMANCE FIXED SETUP ---
        # Initialize the lines ONCE with dummy zeros. We only manipulate line properties later.
        self.line_ol_mag, = self.ax1.semilogx(self.f, np.zeros_like(self.f), color='firebrick', lw=2)
        self.line_ol_phase, = self.ax3.semilogx(self.f, np.zeros_like(self.f), color='firebrick', lw=2)
        self.line_cl_mag, = self.ax2.semilogx(self.f, np.zeros_like(self.f), color='teal', lw=2)
        self.line_cl_phase, = self.ax4.semilogx(self.f, np.zeros_like(self.f), color='teal', lw=2)
        
        # Add crossover tracking marker line
        self.v_marker = self.ax3.axvline(1000, color='purple', linestyle=':', alpha=0.7)

        # Set static limits and labels once so they don't repaint dynamically
        self.ax1.axhline(0, color='black', linestyle='--', alpha=0.5)
        self.ax3.axhline(-180, color='black', linestyle='--', alpha=0.5)
        
        self.ax1.set_title('Open-Loop Return Ratio Gain $|T(f)|$')
        self.ax1.set_ylabel('Gain (dB)')
        self.ax1.set_ylim(-60, 60)
        
        self.ax3.set_title('Open-Loop Phase $\\angle T(f)$')
        self.ax3.set_ylabel('Phase (deg)')
        self.ax3.set_xlabel('Frequency (Hz)')
        self.ax3.set_ylim(-270, 90)
        
        self.ax2.set_title('Closed-Loop Gain $|V_{sense}/V_{cmd}|$')
        self.ax2.set_ylabel('Gain (dB [V/V])')
        self.ax2.set_ylim(-40, 20)
        
        self.ax4.set_title('Closed-Loop Phase Tracking')
        self.ax4.set_ylabel('Phase (deg)')
        self.ax4.set_xlabel('Frequency (Hz)')
        self.ax4.set_ylim(-190, 190)

        for ax in [self.ax1, self.ax2, self.ax3, self.ax4]:
            ax.set_xlim(10, 3e6)
            ax.grid(True, which="both", ls="-", alpha=0.3)
            
        self.fig.tight_layout()

    def read_sliders(self):
        for key, (slider, p_min, scale, _) in self.sliders.items():
            step = slider.value()
            if key in ['L1', 'R1']:
                self.params[key] = (p_min + step) * 0.1
                self.val_labels[key][0].setText(f"{self.params[key]:.2f} {self.val_labels[key][1]}")
            else:
                self.params[key] = p_min + (step * scale)
                self.val_labels[key][0].setText(f"{int(self.params[key])} {self.val_labels[key][1]}")

    def update_plots(self):
        self.read_sliders()
        
        # Parameter Unpacking
        R1 = self.params['R1']
        L1 = self.params['L1'] * 1e-3
        R2 = self.params['R2']
        R3, R4, R5, R6 = 1000.0, 1000.0, 10000.0, 10000.0
        R10 = self.params['R10']
        R11 = self.params['R11']
        R12 = self.params['R12'] * 1e3
        C1 = self.params['C1'] * 1e-9
        C2 = self.params['C2'] * 1e-12
        
        omega_U3 = (self.params['GBW_U3'] * 1e6) * 2 * np.pi
        omega_XU2 = (self.params['GBW_XU2'] * 1e6) * 2 * np.pi
        omega_U1 = 3.0 * 1e6 * 2 * np.pi
        
        # --- Nodal Loop Shape Vector Math ---
        Z_c1 = 1.0 / (self.s * C1)
        Z_c2 = 1.0 / (self.s * C2)
        Z_f = ((R12 + Z_c1) * Z_c2) / (R12 + Z_c1 + Z_c2)
        
        A_u3 = omega_U3 / self.s
        A_u1 = omega_U1 / self.s
        A_xu2 = omega_XU2 / self.s
        
        A_v_xu2 = 1.0 + (14700.0 / 1000.0)
        G_inner = A_xu2 / (1.0 + A_xu2 * (1.0 / A_v_xu2))
        
        G_plant = 1.0 / (self.s * L1 + R1 + R2)
        G_sense = (R2 * (R6 / R3)) * (A_u1 / (1.0 + A_u1))
        G_forward_path = G_inner * G_plant * G_sense
        
        noise_gain = 1.0 + (Z_f / R11) + (Z_f / R10)
        denom_u3 = 1.0 + (noise_gain / A_u3)
        
        G_vcmd_to_u3 = -(Z_f / R11) / denom_u3
        G_vsense_to_u3 = -(Z_f / R10) / denom_u3
        
        T_loop = -G_vsense_to_u3 * G_forward_path
        sys_cl = (G_vcmd_to_u3 * G_forward_path) / (1.0 - G_vsense_to_u3 * G_forward_path)
        
        # Extract Magnitudes and Phases
        mag_ol = 20 * np.log10(np.abs(T_loop))
        phase_ol = np.unwrap(np.degrees(np.angle(T_loop)), period=360)
        if phase_ol[0] > 90: phase_ol -= 180
        if phase_ol[0] < -90: phase_ol += 180
            
        mag_cl = 20 * np.log10(np.abs(sys_cl))
        phase_cl = np.unwrap(np.degrees(np.angle(sys_cl)), period=360)
        
        # Metrics Calculation
        idx_0dB = np.argmin(np.abs(mag_ol))
        f_c = self.f[idx_0dB]
        pm = phase_ol[idx_0dB] - (-180)
        
        # Update text readout box layout without lagging frame triggers
        self.status_box.setText(
            f"=== LOOP SHAPE METRICS ===\n"
            f"Loop Crossover (fc) : {f_c:.1f} Hz\n"
            f"Phase Margin (PM)   : {pm:.2f}°\n"
            f"Target DC Attn floor: {20*np.log10(R10/R11):.2f} dB\n"
            f"Simulated DC Gain   : {mag_cl[0]:.2f} dB\n"
            f"Transconductance    : {10**(mag_cl[0]/20):.3f} A/V"
        )
        
        # --- THE SPEED FIX: Mutate existing line vector properties instantly ---
        self.line_ol_mag.set_ydata(mag_ol)
        self.line_ol_phase.set_ydata(phase_ol)
        self.line_cl_mag.set_ydata(mag_cl)
        self.line_cl_phase.set_ydata(phase_cl)
        
        # Move crossover visual marker line location
        self.v_marker.set_xdata([f_c, f_c])
        
        # Redraw canvas background instantly without reconstruction steps
        self.canvas.draw_idle()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AmplifierSimGUI()
    window.show()
    sys.exit(app.exec_())