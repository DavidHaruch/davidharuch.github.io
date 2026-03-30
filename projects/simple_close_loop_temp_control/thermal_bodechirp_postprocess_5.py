import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import detrend
from scipy.optimize import differential_evolution

# --- 1. Model Definitions (Includes ZOH Delay) ---
def fopdt_model(f, Kp, tau, L, fs=40.0):
    """First-Order Plus Dead Time: G(s) = (Kp * e^-Ls) / (tau*s + 1)"""
    s = 2j * np.pi * f
    total_L = L + (1.0 / (2.0 * fs)) # Add 1/2 sample period ZOH delay
    return (Kp * np.exp(-total_L * s)) / (tau * s + 1)

def sopdt_model(f, Kp, tau1, tau2, L, fs=40.0):
    """Second-Order Plus Dead Time: G(s) = (Kp * e^-Ls) / ((tau1*s+1)(tau2*s+1))"""
    s = 2j * np.pi * f
    total_L = L + (1.0 / (2.0 * fs)) # Add 1/2 sample period ZOH delay
    denominator = (tau1 * s + 1) * (tau2 * s + 1)
    return (Kp * np.exp(-total_L * s)) / denominator

# --- 2. Global Objective Functions ---
def objective_fopdt(params, f_data, actual_combined, fs):
    Kp, tau, L = params
    G = fopdt_model(f_data, Kp, tau, L, fs)
    mag_db = 20 * np.log10(np.clip(np.abs(G), 1e-9, None))
    phase_deg = np.degrees(np.angle(G))
    return np.mean((actual_combined - np.concatenate([mag_db, phase_deg]))**2)

def objective_sopdt(params, f_data, actual_combined, fs):
    Kp, tau1, tau2, L = params
    G = sopdt_model(f_data, Kp, tau1, tau2, L, fs)
    mag_db = 20 * np.log10(np.clip(np.abs(G), 1e-9, None))
    phase_deg = np.degrees(np.angle(G))
    return np.mean((actual_combined - np.concatenate([mag_db, phase_deg]))**2)

# --- 3. Main Analysis Suite ---
def analyze_thermal_performance(filename, f_start=0.0005, f_end=0.2, fs=40.0):
    # --- Data Loading & Prep ---
    try:
        data = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: {filename} not found."); return

    u = data['u'].values  # Heater Duty (0-1)
    y = data['y'].values  # Feedback Sensor
    d = data['d'].values  # Ambient Sensor
    
    t_sec = np.arange(len(u)) / fs
    t_hrs = t_sec / 3600
    y_iso = y - d # Isolated plant response
    
    # Detrend for Bode extraction (removes non-oscillatory thermal drift)
    u_dt = detrend(u)
    y_iso_dt = detrend(y_iso)

    # --- 4. Time Domain Visualization (3 Subplots) ---
    fig_t, (t_ax1, t_ax2, t_ax3) = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    
    t_ax1.plot(t_hrs, y, 'b', label='Feedback Temp (y)', alpha=0.8)
    t_ax1.plot(t_hrs, d, 'r', label='Ambient Temp (d)', alpha=0.6)
    t_ax1.set_ylabel('Temp (°C)')
    t_ax1.set_title(f'Time Domain Signal Integrity: {filename}')
    t_ax1.legend(loc='upper right'); t_ax1.grid(True, alpha=0.3)

    t_ax2.plot(t_hrs, u * 100, 'g', label='Heater Duty Cycle (%)')
    t_ax2.set_ylabel('Duty %')
    t_ax2.legend(loc='upper right'); t_ax2.grid(True, alpha=0.3)

    t_ax3.plot(t_hrs, y_iso, 'purple', label='Isolated Plant (y - d)')
    t_ax3.set_ylabel('Δ Temp (°C)')
    t_ax3.set_xlabel('Time (hours)')
    t_ax3.legend(loc='upper right'); t_ax3.grid(True, alpha=0.3)

    # --- 5. Frequency Domain (Bode) Extraction ---
    freqs = np.logspace(np.log10(f_start), np.log10(f_end), 80)
    mags, phases = [], []

    print(f"Extracting Bode response for {len(freqs)} points...")
    for f in freqs:
        ref = np.exp(-1j * 2 * np.pi * f * t_sec)
        U_f = np.sum(u_dt * ref)
        Y_f = np.sum(y_iso_dt * ref)
        G_f = Y_f / U_f
        mags.append(20 * np.log10(np.abs(G_f)))
        phases.append(np.degrees(np.angle(G_f)))
    
    actual_combined = np.concatenate([mags, phases])

    # --- 6. Dual Global Optimization ---
    print("\n[1/2] Fitting First-Order Model (FOPDT)...")
    # Bounds: [Kp, tau, L]
    res1 = differential_evolution(objective_fopdt, 
                                  [(0.1, 20), (10, 500), (0, 100)], 
                                  args=(freqs, actual_combined, fs), 
                                  polish=True)
    
    print("[2/2] Fitting Second-Order Model (SOPDT)...")
    # Bounds: [Kp, tau1, tau2, L]
    res2 = differential_evolution(objective_sopdt, 
                                  [(0.1, 20), (30, 400), (0.1, 150), (0, 60)], 
                                  args=(freqs, actual_combined, fs), 
                                  polish=True)

    # Generate Model Results
    G1 = fopdt_model(freqs, *res1.x, fs=fs)
    G2 = sopdt_model(freqs, *res2.x, fs=fs)

    # --- 7. Bode Comparison Plot ---
    fig_b, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    
    # Magnitude
    ax1.semilogx(freqs, mags, 'ko', alpha=0.2, label='Measured Data')
    ax1.semilogx(freqs, 20*np.log10(np.abs(G1)), 'r--', label=f'FOPDT: τ={res1.x[1]:.1f}s, L={res1.x[2]:.1f}s')
    ax1.semilogx(freqs, 20*np.log10(np.abs(G2)), 'b-', linewidth=2, label=f'SOPDT: τ1={res2.x[1]:.1f}s, τ2={res2.x[2]:.1f}s, L={res2.x[3]:.1f}s')
    ax1.set_ylabel('Magnitude (dB)')
    ax1.set_title('Bode Model Comparison: Heater to Feedback')
    ax1.grid(True, which="both", alpha=0.2); ax1.legend()

    # Phase
    ax2.semilogx(freqs, phases, 'ko', alpha=0.2)
    ax2.semilogx(freqs, np.degrees(np.angle(G1)), 'r--')
    ax2.semilogx(freqs, np.degrees(np.angle(G2)), 'b-', linewidth=2)
    ax2.set_ylabel('Phase (degrees)')
    ax2.set_xlabel('Frequency (Hz)')
    ax2.grid(True, which="both", alpha=0.2)
    ax2.set_ylim(-270, 0)
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nModel Evaluation:")
    print(f"FOPDT Final MSE: {res1.fun:.4f}")
    print(f"SOPDT Final MSE: {res2.fun:.4f}")

# --- Execution ---
analyze_thermal_performance('sweep_data_1774731446.csv')