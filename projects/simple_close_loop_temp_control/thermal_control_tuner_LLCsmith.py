import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --- 1. Parameters (FOPDT Identified + Lambda Tuned) ---
KP_PLANT = 10.7977  
TAU_PLANT = 30.00   
L_PLANT = 1.76      
FS = 40.0           
TS = 1.0 / FS

# Lambda Tuning Parameters for Smith Predictor
# With Smith, we can be much more aggressive (smaller tau_cl)
tau_cl = 0.5  # Reduced from 1.5 because delay is compensated
Kc = (1.0 / KP_PLANT) * (TAU_PLANT / tau_cl) # Note: L_PLANT is removed from denominator
Ti = TAU_PLANT

Kp_reg = Kc
Ki_reg = Kc / Ti

# --- 2. Transfer Function Engine ---
def get_responses(f, use_smith=True):
    s = 2j * np.pi * f
    L_total = L_PLANT + (TS / 2.0)
    
    # Actual Plant G(s)
    G = (KP_PLANT * np.exp(-L_total * s)) / (TAU_PLANT * s + 1)
    
    # Internal Model (Delay-free) G_model(s)
    G_model_no_delay = KP_PLANT / (TAU_PLANT * s + 1)
    
    # PI Controller C_pi(s)
    C_pi = Kp_reg + (Ki_reg / s)
    
    if use_smith:
        # Smith Predictor Equivalent Controller
        # C_eq = C_pi / (1 + C_pi * G_model_no_delay * (1 - exp(-Ls)))
        feedback_model = G_model_no_delay * (1 - np.exp(-L_total * s))
        C_eff = C_pi / (1 + C_pi * feedback_model)
    else:
        C_eff = C_pi
        
    # Open Loop L(s) = C_eff * G
    L = C_eff * G
    # Closed Loop T(s)
    T = L / (1 + L)
    
    return G, L, T

def wrap_phase(data_rad):
    deg = np.degrees(np.angle(data_rad))
    return (deg + 360) % 360 - 360# Centered around -180 for stability view

# --- 3. Extracting Metrics ---
def find_metrics(use_smith=True):
    # Crossover: |L(f)| = 1
    try:
        f_gc = fsolve(lambda f: np.abs(get_responses(f, use_smith)[1]) - 1.0, 0.1)[0]
        _, L_gc, _ = get_responses(f_gc, use_smith)
        pm = 180 + np.degrees(np.angle(L_gc))
    except:
        f_gc, pm = 0, 0
    
    # Bandwidth: |T(f)| = 0.707
    try:
        f_bw = fsolve(lambda f: 20*np.log10(np.abs(get_responses(f, use_smith)[2])) + 3.0, 0.2)[0]
    except:
        f_bw = 0
        
    return f_gc, pm, f_bw

# --- 4. Plotting ---
freqs = np.logspace(-3, 0, 1000)
G, L_s, T_s = get_responses(freqs, use_smith=True)
# _, L_pi, T_pi = get_responses(freqs, use_smith=False)

f_gc, pm, f_bw = find_metrics(use_smith=True)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
fig.patch.set_facecolor('#F0F0F0')

# Magnitude
# ax1.semilogx(freqs, 20*np.log10(np.abs(L_pi)), 'b--', alpha=0.4, label='Standard PI Open Loop')
ax1.semilogx(freqs, 20*np.log10(np.abs(L_s)), 'blue', linewidth=2, label='Smith Predictor Open Loop')
ax1.semilogx(freqs, 20*np.log10(np.abs(T_s)), 'green', linewidth=2, label='Smith Predictor Closed Loop')
ax1.axhline(0, color='black', linewidth=1, linestyle='--')

# Phase
# ax2.semilogx(freqs, wrap_phase(L_pi), 'b--', alpha=0.4)
ax2.semilogx(freqs, wrap_phase(L_s), 'blue', linewidth=2)
ax2.semilogx(freqs, wrap_phase(T_s), 'green', linewidth=2)
ax2.axhline(-180, color='red', linewidth=1, linestyle='--')

# Summary Box
stats_text = (f'SMITH PREDICTOR METRICS\n'
              f'$\lambda$: {tau_cl}s\n'
              f'Kp: {Kp_reg:.4f}\n'
              f'Ki: {Ki_reg:.4f}\n'
              f'PM: {pm:.1f}°\n'
              f'BW: {f_bw:.3f} Hz\n'
              f'Crossover: {f_gc:.3f} Hz')
ax1.text(0.05, 0.05, stats_text, transform=ax1.transAxes, bbox=dict(facecolor='wheat', alpha=1))

ax1.set_ylabel("Magnitude (dB)")
ax1.set_title("Bode Analysis: Smith Predictor + PI")
ax1.legend()
ax1.grid(True, which="both", alpha=0.3)
ax2.set_ylabel("Phase (deg)")
ax2.set_xlabel("Frequency (Hz)")
ax2.grid(True, which="both", alpha=0.3)
ax2.set_ylim(-360, 0)
ax2.set_xlim(0.001, 1)

plt.tight_layout()
plt.show()