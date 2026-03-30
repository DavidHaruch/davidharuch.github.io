import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit

# --- 1. CORE MODELS & MATH ---

def fopdt_model(f, K, tau, theta):
    """Frequency response of a First-Order Plus Dead Time plant."""
    s = 2j * np.pi * f
    return (K * np.exp(-s * theta)) / (tau * s + 1)

def fit_objective(f, K, tau, theta):
    """Objective function fitting Magnitude in dB."""
    G = fopdt_model(f, K, tau, theta)
    return 20 * np.log10(np.abs(G))

def analyze_system(filename, kp=20.0, ki=0.015, kd=0.0, fs=40.0):
    # --- 2. DATA LOADING & CLEANING ---
    try:
        data = pd.read_csv(filename)
        u_raw = data['u'].values
        y_delta = data['y'].values - data['d'].values # Ambient subtraction
        u = signal.detrend(u_raw)
        y = signal.detrend(y_delta)
    except Exception as e:
        print(f"Error: {e}"); return

    # --- 3. SPECTRAL ESTIMATION ---
    nperseg = int(fs / 0.001) * 2  
    f, Puu = signal.welch(u, fs, nperseg=nperseg)
    f, Pyy = signal.welch(y, fs, nperseg=nperseg)
    f, Puy = signal.csd(u, y, fs, nperseg=nperseg)
    
    valid = (f > 0.0005) & (f < 1.0)
    f_v = f[valid]
    G_raw = Puy[valid] / Puu[valid]
    coh_v = (np.abs(Puy[valid])**2) / (Puu[valid] * Pyy[valid])
    mag_db_raw = 20 * np.log10(np.abs(G_raw))

    # --- 4. ROBUST MODEL FITTING ---
    # We weight the fit by coherence: trust clear signals more than noise
    weights = coh_v**2 
    
    try:
        # Bounds: K [0, 50], tau [1, 1000], theta [0, 20]
        popt, _ = curve_fit(
            fit_objective, f_v, mag_db_raw, 
            p0=[0.5, 60, 1.5], 
            sigma=1/weights, 
            bounds=([0.01, 1, 0], [100, 2000, 50])
        )
        K_m, tau_m, theta_m = popt
        G_model = fopdt_model(f_v, K_m, tau_m, theta_m)
        print(f"\nSUCCESSFUL FIT:\nK: {K_m:.3f}\ntau: {tau_m:.1f}s\ntheta: {theta_m:.3f}s")
    except Exception as e:
        print(f"Fit failed: {e}"); G_model = None

    # --- 5. LOOP MATH ---
    s = 2j * np.pi * f_v
    C = kp + (ki / s) + (kd * s)
    L_raw = C * G_raw
    L_mod = C * G_model if G_model is not None else None

    # --- 6. FIGURE 1: PLANT MODEL OVERLAY (The "Fit Check") ---
    plt.figure(figsize=(10, 8))
    plt.subplot(2, 1, 1)
    plt.semilogx(f_v, mag_db_raw, 'k.', alpha=0.3, label='Experimental Data')
    if G_model is not None:
        plt.semilogx(f_v, 20*np.log10(np.abs(G_model)), 'r-', lw=2, label='FOPDT Fit')
    plt.ylabel('Plant Magnitude (dB)'); plt.legend(); plt.grid(True, which='both', alpha=0.3)
    plt.title("Figure 1: Plant Model Identification Accuracy")

    plt.subplot(2, 1, 2)
    plt.semilogx(f_v, np.rad2deg(np.unwrap(np.angle(G_raw))), 'k.', alpha=0.3)
    if G_model is not None:
        plt.semilogx(f_v, np.rad2deg(np.unwrap(np.angle(G_model))), 'r-', lw=2)
    plt.ylabel('Plant Phase (deg)'); plt.xlabel('Frequency (Hz)'); plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout(); plt.show()

    # --- 7. FIGURE 2: SYSTEM STABILITY & SENSITIVITY ---
    if L_mod is not None:
        fig2, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
        S_mod = 1 / (1 + L_mod)
        T_mod = L_mod / (1 + L_mod)
        
        ax_mag.semilogx(f_v, 20*np.log10(np.abs(L_mod)), 'cyan', label='Open Loop (L)')
        ax_mag.semilogx(f_v, 20*np.log10(np.abs(S_mod)), 'tomato', lw=2, label='Disturbance Rejection (S)')
        ax_mag.semilogx(f_v, 20*np.log10(np.abs(T_mod)), 'lime', lw=2, label='Tracking (T)')
        ax_mag.axhline(0, color='black', lw=1); ax_mag.set_ylabel('Magnitude (dB)')
        ax_mag.set_title("Figure 2: Modeled Loop Performance"); ax_mag.legend(); ax_mag.grid(True, which='both', alpha=0.3)

        phase_mod = np.rad2deg(np.unwrap(np.angle(L_mod)))
        ax_phase.semilogx(f_v, phase_mod, 'cyan', lw=2)
        ax_phase.axhline(-180, color='red', ls='--'); ax_phase.set_ylabel('Phase (deg)')
        ax_phase.grid(True, which='both', alpha=0.3); ax_phase.set_xlabel('Frequency (Hz)')
        
        # Stability Margins
        idx_u = np.argmin(np.abs(np.abs(L_mod) - 1.0))
        pm = 180 + phase_mod[idx_u]
        ms_peak = np.max(20 * np.log10(np.abs(S_mod)))
        print(f"\n--- PERFORMANCE ---\nCrossover: {f_v[idx_u]:.4f} Hz\nPhase Margin: {pm:.2f}°\nMax Sensitivity (Ms): {ms_peak:.2f} dB\n-------------------")
        
    plt.tight_layout(); plt.show()

if __name__ == "__main__":
    analyze_system("bode_data_1774504507.csv", kp=20, ki=0.015)