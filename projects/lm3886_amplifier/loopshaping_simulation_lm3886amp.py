import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# =====================================================================
# 1. Hardware Inputs & Target Loop Shaping Specifications
# =====================================================================
# Physical Plant Parameters
R_coil = 0.6889925-0.1         # Coil resistance (Ohms)
L_coil = 0.1486e-3        # Coil inductance (Henries, e.g., 2mH = 0.002)
R_shunt = 0.1         # Shunt resistor (Ohms)
R_total = R_coil + R_shunt

# Target Bandwidth selection (e.g., 1500 Hz standard reliable range)
Desired_BW_Hz = 1000.0 
R_in = 33000.0        # Anchor input resistor (10 kOhm)

# Op-Amp Non-Ideal Specifications
GBWP_op27 = 8.0e6     # 8 MHz
A0_op27 = 1.8e6       # 125 dB Open Loop DC Gain
F_c_diffamp = 800e3   # 800 kHz (OP27 DiffAmp stage corner)

Gain_LM3886_dc = 15.0
GBWP_lm3886 = 8.0e6
F_c_lm3886 = GBWP_lm3886 / Gain_LM3886_dc # ~533.3 kHz
H_gain = R_shunt * 10.0 # Feedback Shunt + DiffAmp nominal scalar gain (1.0 V/A)

# =====================================================================
# 2. Automated Component Quantization Math (E96 / E24 Standard Stock)
# =====================================================================
e96_base = np.array([
    1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30, 
    1.33, 1.37, 1.40, 1.43, 1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74, 
    1.78, 1.82, 1.87, 1.91, 1.96, 2.00, 2.05, 2.10, 2.15, 2.21, 2.26, 2.32, 
    2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09, 
    3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12, 
    4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49, 
    5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81, 6.98, 7.15, 7.32, 
    7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76
])

e24_base = np.array([
    1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 
    3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1
])

def find_closest_standard(target, base_values):
    exponent = np.floor(np.log10(target))
    remainder = target / (10 ** exponent)
    idx = np.argmin(np.abs(base_values - remainder))
    return base_values[idx] * (10 ** exponent)

# Calculate pole-zero cancellation targets
omega_c = 2 * np.pi * Desired_BW_Hz
Kp_ideal = (omega_c * L_coil) / (Gain_LM3886_dc * H_gain)
Ki_ideal = Kp_ideal * (R_total / L_coil)

# Standard bench value snapping
R_f = find_closest_standard(Kp_ideal * R_in, e96_base)
C_i = find_closest_standard(1.0 / (Ki_ideal * R_in), e24_base)

print(f"--- Automated Hardware Loop Selection ---")
print(f"  * Targeted Bandwidth:  {Desired_BW_Hz:.1f} Hz")
print(f"  * Selected R_in:       {R_in/1000:.1f} kΩ")
print(f"  * Selected R_f (1%):   {R_f/1000:.3f} kΩ (Ideal: {Kp_ideal*R_in/1000:.3f} kΩ)")
print(f"  * Selected C_i (5%):   {C_i*1e9:.1f} nF (Ideal: {1.0/(Ki_ideal*R_in)*1e9:.1f} nF)")

# =====================================================================
# 3. Non-Ideal Transfer Function Derivations (Using Snapped Components)
# =====================================================================
# --- Non-Ideal OP27 PI Controller C(s) ---
num_Zf = [R_f * C_i, 1.0]
den_Zf = [R_in * C_i, 0.0]

tau_ol = A0_op27 / (2.0 * np.pi * GBWP_op27)
num_As = [A0_op27]
den_As = [tau_ol, 1.0]

# Exact closed loop matrix math: C_real(s) = Zf / (Zin + (Zf + Zin)/A(s))
num_C_real = np.polymul(num_Zf, num_As)
p1 = np.polymul(den_Zf, num_As)
p2 = np.polymul(den_Zf, den_As)
p3 = np.polymul(num_Zf, den_As)
den_C_real = np.polyadd(p1, np.polyadd(p2, p3))
C_real = signal.TransferFunction(num_C_real, den_C_real)

# --- Non-Ideal LM3886 Power Amp ---
num_Amp = [Gain_LM3886_dc]
den_Amp = [1.0 / (2.0 * np.pi * F_c_lm3886), 1.0]
Amp_real = signal.TransferFunction(num_Amp, den_Amp)

# --- Real Plant Combination: P_real(s) = Amp_real(s) * Coil(s) ---
num_P_real = num_Amp
den_P_real = np.polymul(den_Amp, [L_coil, R_total])
P_real = signal.TransferFunction(num_P_real, den_P_real)

# --- Real Feedback Path: H_real(s) ---
tau_diff = 1.0 / (2.0 * np.pi * F_c_diffamp)
H_real = signal.TransferFunction([H_gain], [tau_diff, 1.0])

# =====================================================================
# 4. System Loops Construction
# =====================================================================
# Open Loop: L(s) = C_real * P_real * H_real
num_L = np.polymul(np.polymul(num_C_real, num_P_real), [H_gain])
den_L = np.polymul(np.polymul(den_C_real, den_P_real), [tau_diff, 1.0])
L = signal.TransferFunction(num_L, den_L)

# Closed Loop Tracking: T(s)
num_T = np.polymul(np.polymul(num_C_real, num_P_real), [tau_diff, 1.0])
den_T = np.polyadd(den_L, num_L)
T = signal.TransferFunction(num_T, den_T)

# Sensitivity: S(s)
S = signal.TransferFunction(den_L, den_T)

# =====================================================================
# 5. Frequency Analysis & Plotting
# =====================================================================
w = np.logspace(0, 7, 2500) * 2 * np.pi
freq_hz = w / (2 * np.pi)

def get_bode(tf, w_vec):
    _, mag, phase = signal.bode(tf, w_vec)
    return mag, phase

mag_P, phase_P = get_bode(P_real, w)
mag_C, phase_C = get_bode(C_real, w)
mag_H, phase_H = get_bode(H_real, w)
mag_L, phase_L = get_bode(L, w)
mag_T, phase_T = get_bode(T, w)
mag_S, phase_S = get_bode(S, w)

idx_crossover = np.argmin(np.abs(mag_L))
crossover_freq = freq_hz[idx_crossover]
phase_margin = 180 + phase_L[idx_crossover]

print(f"\n--- Realized System Performance Metrics ---")
print(f"  * Realized 0 dB Crossover Frequency: {crossover_freq:.1f} Hz")
print(f"  * Realized Phase Margin:             {phase_margin:.1f} degrees")

# Generate 6x2 Plotting Grid
fig, axs = plt.subplots(6, 2, figsize=(14, 18), sharex=True)

def setup_plot(ax, freq, data, title, ylabel, is_mag, color, y_lim=None):
    ax.semilogx(freq, data, color=color, linewidth=2)
    ax.grid(True, which="both", linestyle=':', color='gray', alpha=0.6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10, fontweight='bold')
    if y_lim: ax.set_ylim(y_lim)

# --- Subplot Layout Setup ---
setup_plot(axs[0,0], freq_hz, mag_P, 'Real Plant P(s) [LM3886 Pole + Coil]', 'Gain (dB)', True, 'orange')
setup_plot(axs[0,1], freq_hz, phase_P, 'Real Plant Phase', 'Phase (deg)', False, 'orange', [-200, 10])

setup_plot(axs[1,0], freq_hz, mag_C, f'Non-Ideal OP27 PI ($R_f={R_f/1000:.1f}k\\Omega, C_i={C_i*1e9:.0f}nF$)', 'Gain (dB)', True, 'purple')
setup_plot(axs[1,1], freq_hz, phase_C, 'Controller Phase', 'Phase (deg)', False, 'purple', [-100, 10])

setup_plot(axs[2,0], freq_hz, mag_H, 'Real Feedback H(s) [OP27 Diff Amp Bandwidth]', 'Gain (dB)', True, 'teal')
setup_plot(axs[2,1], freq_hz, phase_H, 'Feedback Phase', 'Phase (deg)', False, 'teal', [-100, 10])

setup_plot(axs[3,0], freq_hz, mag_L, 'Open Loop L(s) Magnitude', 'Gain (dB)', True, 'crimson')
axs[3,0].axhline(0.0, color='black', linestyle='--')
setup_plot(axs[3,1], freq_hz, phase_L, 'Open Loop Phase', 'Phase (deg)', False, 'crimson', [-280, 10])
axs[3,1].axhline(-180.0, color='darkred', linestyle='--')
axs[3,1].axvline(crossover_freq, color='gray', linestyle='-.')

setup_plot(axs[4,0], freq_hz, mag_T, 'Closed Loop Tracking T(s)', 'Gain (dB)', True, 'royalblue')
setup_plot(axs[4,1], freq_hz, phase_T, 'Closed Loop Phase', 'Phase (deg)', False, 'royalblue', [-280, 10])

setup_plot(axs[5,0], freq_hz, mag_S, 'Sensitivity S(s) [Disturbance Rejection]', 'Gain (dB)', True, 'forestgreen')
axs[5,0].set_xlabel('Frequency (Hz)')
setup_plot(axs[5,1], freq_hz, phase_S, 'Sensitivity Phase', 'Phase (deg)', False, 'forestgreen')
axs[5,1].set_xlabel('Frequency (Hz)')

plt.tight_layout()
plt.show()