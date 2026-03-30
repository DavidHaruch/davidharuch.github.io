import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import detrend
from scipy.optimize import differential_evolution

# --- 1. FOPDT Model Definition ---
def fopdt_model(f, Kp, tau, L, fs=40.0):
    """
    G(s) = (Kp * e^-Ls) / (tau*s + 1)
    tau: Main thermal time constant
    L: Physical transport delay
    """
    s = 2j * np.pi * f
    # Include ZOH delay (Ts/2) + Physical Delay L
    total_L = L + (1.0 / (2.0 * fs)) 
    
    return (Kp * np.exp(-total_L * s)) / (tau * s + 1)

# --- 2. Global Objective Function ---
def objective_fopdt(params, f_data, actual_combined, fs):
    """
    Objective function for differential evolution. 
    params = [Kp, tau, L]
    """
    Kp, tau, L = params
    G = fopdt_model(f_data, Kp, tau, L, fs)
    
    # Magnitude in dB
    mag_db = 20 * np.log10(np.clip(np.abs(G), 1e-9, None))
    # Phase in Degrees
    phase_deg = np.degrees(np.angle(G))
    
    model_combined = np.concatenate([mag_db, phase_deg])
    return np.mean((actual_combined - model_combined)**2)

# --- 3. Main Analysis Suite ---
def analyze_thermal_performance(filename, f_start=0.0005, f_end=0.12, fs=40.0):
    # --- Data Loading ---
    try:
        data = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return

    u = data['u'].values  # Heater Duty (0.0 to 1.0)
    y = data['y'].values  # Feedback Sensor
    d = data['d'].values  # Ambient Sensor
    
    t_sec = np.arange(len(u)) / fs
    t_hrs = t_sec / 3600
    y_isolated = y - d
    
    # Detrend for Bode extraction
    u_dt = detrend(u)
    y_iso_dt = detrend(y_isolated)

    # --- 4. Time Domain Visualization (3 Subplots) ---
    fig_t, (t_ax1, t_ax2, t_ax3) = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    
    fig_t.patch.set_facecolor('#F0F0F0')
    
    t_ax1.plot(t_hrs, y, color='blue', label='Feedback Temp (y)', alpha=0.8)
    t_ax1.plot(t_hrs, d, color='red', label='Ambient Temp (d)', alpha=0.6)
    t_ax1.set_ylabel('Temp (°C)')
    t_ax1.set_title(f'Time Domain Signal Integrity: {filename}')
    t_ax1.legend(loc='best'); t_ax1.grid(True, alpha=0.3)

    t_ax2.plot(t_hrs, u * 100, color='green', label='Heater Duty (%)')
    t_ax2.set_ylabel('Duty %')
    t_ax2.legend(loc='upper right'); t_ax2.grid(True, alpha=0.3)

    t_ax3.plot(t_hrs, y_isolated, color='purple', label='Isolated Plant (y - d)')
    t_ax3.set_ylabel('Δ Temp (°C)')
    t_ax3.set_xlabel('Time (hours)')
    t_ax3.legend(loc='upper right'); t_ax3.grid(True, alpha=0.3)

    # --- 5. Frequency Domain (Bode) Extraction ---
    # Log-spaced frequencies up to the 0.12Hz cutoff
    freqs = np.logspace(np.log10(f_start), np.log10(f_end), 80)
    mags, phases = [], []

    print(f"Extracting Bode response for {len(freqs)} points up to {f_end} Hz...")
    for f in freqs:
        ref = np.exp(-1j * 2 * np.pi * f * t_sec)
        U_f = np.sum(u_dt * ref)
        Y_f = np.sum(y_iso_dt * ref)
        G_f = Y_f / U_f
        
        mags.append(20 * np.log10(np.abs(G_f)))
        phases.append(np.degrees(np.angle(G_f)))

    actual_combined = np.concatenate([mags, phases])

    # --- 6. Global Optimization (Differential Evolution) ---
    print("\nTuning FOPDT Model (Plant + Dead Time)...")
    # Bounds: [Kp, tau, L]
    search_bounds = [
        (0.1, 20.0),    # Kp: Static Gain
        (30.0, 400.0),  # tau: Time Constant
        (0.0, 100.0)    # L: Dead Time (s)
    ]
    
    result = differential_evolution(
        objective_fopdt, 
        bounds=search_bounds, 
        args=(freqs, actual_combined, fs),
        strategy='best1bin', 
        maxiter=1000,
        popsize=15,
        tol=0.001,
        polish=True
    )

    Kp_f, tau_f, L_f = result.x
    G_fit = fopdt_model(freqs, Kp_f, tau_f, L_f, fs)
    mags_fit = 20 * np.log10(np.abs(G_fit))
    phases_fit = np.degrees(np.angle(G_fit))

    # --- 7. Final Plotting ---
    fig_b, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    
    fig_b.patch.set_facecolor('#F0F0F0')
    
    # Magnitude Plot
    ax1.semilogx(freqs, mags, 'ko', alpha=0.2, label='Measured Data')
    ax1.semilogx(freqs, mags_fit, 'r-', linewidth=2, 
                 label=f'FOPDT Fit: Kp={Kp_f:.3f}, τ={tau_f:.1f}s, L={L_f:.1f}s')
    ax1.set_ylabel('Magnitude (dB)')
    ax1.set_title(f'Bode Analysis & FOPDT Fit: {filename}')
    ax1.grid(True, which="both", alpha=0.2); ax1.legend()

    # Phase Plot
    ax2.semilogx(freqs, phases, 'ko', alpha=0.2)
    ax2.semilogx(freqs, phases_fit, 'r-', linewidth=2)
    ax2.set_ylabel('Phase (degrees)')
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylim(-210, 0)
    ax2.grid(True, which="both", alpha=0.2)
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nFinal Identified Parameters:")
    print(f"  Plant Gain (Kp):      {Kp_f:.4f} °C/Unit Duty")
    print(f"  Time Constant (τ):    {tau_f:.2f} seconds")
    print(f"  Dead Time (L):        {L_f:.2f} seconds")
    print(f"  Final MSE:            {result.fun:.4f}")

# --- Execution ---
analyze_thermal_performance('sweep_data_1774731446.csv')