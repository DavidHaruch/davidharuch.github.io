import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --- 1. Parameters (FOPDT Identified + Lambda Tuned) ---
KP_PLANT = 10.7977  
TAU_PLANT = 30.00   
L_PLANT = 1.76      
FS = 40.0           
TS = 1.0 / FS

# Lambda Tuning Parameters
# tau_cl = max(L_PLANT, 1.0)
tau_cl = 1.5
Kc = (1.0 / KP_PLANT) * (TAU_PLANT / (tau_cl + L_PLANT))
Ti = TAU_PLANT

# Parallel Form Gains for Implementation
Kp_reg = Kc
Ki_reg = Kc / Ti

# --- 2. Transfer Function Engine ---
def get_responses(f):
    s = 2j * np.pi * f
    L_total = L_PLANT + (TS / 2.0)
    
    # Plant G(s)
    G = (KP_PLANT * np.exp(-L_total * s)) / (TAU_PLANT * s + 1)
    # Controller C(s)
    C = Kp_reg + (Ki_reg / s)
    # Open Loop L(s)
    L = C * G
    # Closed Loop T(s)
    T = L / (1 + L)
    
    return G, L, T

def wrap_phase(data_rad):
    """Wraps phase to be strictly between 0 and -360 degrees."""
    deg = np.degrees(np.angle(data_rad))
    return (deg + 360) % 360 - 360

# --- 3. Extracting Metrics ---
def find_metrics():
    # Crossover: |L(f)| = 1 (0 dB)
    f_gc = fsolve(lambda f: np.abs(get_responses(f)[1]) - 1.0, 0.05)[0]
    _, L_gc, _ = get_responses(f_gc)
    pm = 180 + np.degrees(np.angle(L_gc))
    
    # Bandwidth: |T(f)| = 0.707 (-3 dB)
    f_bw = fsolve(lambda f: 20*np.log10(np.abs(get_responses(f)[2])) + 3.0, 0.1)[0]
    
    return f_gc, pm, f_bw

# --- 4. Plotting Construction ---
freqs = np.logspace(-3, 0, 1000)
G, L, T = get_responses(freqs)
f_gc, pm, f_bw = find_metrics()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

fig.patch.set_facecolor('#F0F0F0')

# Define Legend Labels with Gains
controller_label = f'PI Gains: Kp={Kp_reg:.4f}, Ki={Ki_reg:.4f}'

# --- Magnitude Plot ---
ax1.semilogx(freqs, 20*np.log10(np.abs(G)), 'gray', alpha=0.5, label='Plant G(s)')
ax1.semilogx(freqs, 20*np.log10(np.abs(L)), 'blue', linewidth=2, label='Open Loop L(s)')
ax1.semilogx(freqs, 20*np.log10(np.abs(T)), 'green', linewidth=2, label='Closed Loop T(s)')

# Annotations (Magnitude)
ax1.axhline(0, color='black', linewidth=1, linestyle='--')
# # ax1.axhline(-3, color='green', linewidth=0.8, linestyle=':', label='-3dB Bandwidth')
# ax1.annotate(f'Crossover: {f_gc:.3f} Hz', xy=(f_gc, 0), xytext=(f_gc*1.5, 10),
#              arrowprops=dict(arrowstyle='->', color='blue'))
# ax1.annotate(f'Bandwidth: {f_bw:.3f} Hz', xy=(f_bw, -3), xytext=(f_bw*1.5, -15),
#              arrowprops=dict(arrowstyle='->', color='green'))

# Summary Box
stats_text = (f'TUNING PARAMETERS\n'
              f'Kp: {Kp_reg:.4f}\n'
              f'Ki: {Ki_reg:.4f}\n'
              f'PM: {pm:.1f}°\n'
              f'BW: {f_bw:.3f} Hz\n'
              f'Crossover: {f_gc:.3f} Hz')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.5, 0.95, stats_text, transform=ax1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

KP_PLANT = 10.7977  
TAU_PLANT = 30.00   
L_PLANT = 1.76      
FS = 40.0           
TS = 1.0 / FS

stats_text2 = (f'MODEL PARAMETERS\n'
              f'KP_PLANT: {KP_PLANT:.1f}degC/%\n'
              f'TAU_PLANT: {TAU_PLANT:.1f}s\n'
              f'L_PLANT: {L_PLANT:.2f}s\n'
              f'FS: {FS:.1f} Hz')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.7, 0.95, stats_text2, transform=ax1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax1.set_ylabel("Magnitude (dB)")
ax1.set_title(f"Thermal Control Analysis | {controller_label}")
ax1.legend(loc='lower left')
ax1.grid(True, which="both", alpha=0.3)

# --- Phase Plot ---
ax2.semilogx(freqs, wrap_phase(G), 'gray', alpha=0.5)
ax2.semilogx(freqs, wrap_phase(L), 'blue', linewidth=2)
ax2.semilogx(freqs, wrap_phase(T), 'green', linewidth=2)

# Annotations (Phase)
ax2.axhline(-180, color='red', linewidth=1, linestyle='--')
ax2.set_ylim(-360, 0)
ax2.set_yticks([0, -90, -180, -270, -360])

# Phase Margin Indicator
ax2.annotate('', xy=(f_gc, -180), xytext=(f_gc, wrap_phase(L[np.argmin(np.abs(freqs-f_gc))])),
             arrowprops=dict(arrowstyle='<->', color='red', lw=2))
ax2.text(f_gc * 0.5, -180 + (pm/2), f'PM: {pm:.1f}°', color='red', fontweight='bold')

ax2.set_ylabel("Phase (degrees)")
ax2.set_xlabel("Frequency (Hz)")
ax2.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()

print(f"Summary Table:")
print(f"Kp (Proportional): {Kp_reg:.6f}")
print(f"Ki (Integral):     {Ki_reg:.6f}")
print(f"Phase Margin (PM): {pm:.2f} deg")
print(f"Bandwidth:         {f_bw:.4f} Hz")