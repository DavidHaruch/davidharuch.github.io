import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# -----------------------------------------------------------------------------
# 1. Safe Dutch Mechatronics Style
# -----------------------------------------------------------------------------
plt.rcParams.update({
    "text.usetex": False,            
    "font.family": "sans-serif",
    "mathtext.fontset": "cm",        
    # "font.serif": ["cmr10", "Computer Modern Serif", "DejaVu Serif"], 
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.grid": True,              
    "grid.linestyle": ":",          
    "grid.linewidth": 0.5,
    "grid.color": "#b0b0b0",
    "lines.linewidth": 1.0,         
    "figure.figsize": (5.5, 4.5),   
    "figure.autolayout": True
})

# -----------------------------------------------------------------------------
# 2. Rigol Specific Parser
# -----------------------------------------------------------------------------
def load_rigol_bode(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Could not find the file: {filename}")
    
    df = pd.read_csv(filename, skiprows=4)
    df.columns = df.columns.str.strip()
    
    freq = pd.to_numeric(df.iloc[:, 0], errors='coerce').values
    mag = pd.to_numeric(df.iloc[:, 1], errors='coerce').values
    phase = pd.to_numeric(df.iloc[:, 2], errors='coerce').values
    
    mask = ~np.isnan(freq) & ~np.isnan(mag) & ~np.isnan(phase)
    return freq[mask], mag[mask], phase[mask]

# -----------------------------------------------------------------------------
# 3. Execution
# -----------------------------------------------------------------------------
csv_filename = r"C:\Users\david\Documents\GitHub\davidharuch.github.io\projects\lm3886_amplifier\2026-06-04 exp results\op27gain10noninv0.csv"

csv_filename = r"C:\Users\david\Documents\GitHub\davidharuch.github.io\projects\lm3886_amplifier\2026-06-04 exp results\hihghbwbode june40.csv"


COLOR_MAG = "blue"    
COLOR_PHASE = "blue"  

plot_title = "OP27 Non-Inverting Current Sense Gain 100 Breadboard\n David J. Haruch 2026-06-07"
plot_title = "Closed Current Loop Bode\n David J. Haruch 2026-06-07"

try:
    freq, mag, phase = load_rigol_bode(csv_filename)
    phase = (phase + 180) % 360 - 180

    # -------------------------------------------------------------------------
    # --- Characterization Block ---
    # -------------------------------------------------------------------------
    # 1. DC Gain & -3dB Frequency
    dc_gain = mag[0] 
    target_3db = dc_gain - 3.0
    
    f_3db = None
    if np.any(mag < target_3db):
        idx_below = np.where(mag < target_3db)[0][0]
        if idx_below > 0:
            idx_above = idx_below - 1
            log_f1 = np.log10(freq[idx_above])
            log_f2 = np.log10(freq[idx_below])
            m1 = mag[idx_above]
            m2 = mag[idx_below]
            log_f_3db = log_f1 + (target_3db - m1) * (log_f2 - log_f1) / (m2 - m1)
            f_3db = 10**log_f_3db

    # 2. Peaking, Damping Ratio, and Phase Margin Estimation
    peak_idx = np.argmax(mag)
    f_peak = freq[peak_idx]
    max_mag = mag[peak_idx]
    peaking_db = max_mag - dc_gain

    zeta_est = None
    pm_est = None

    if peaking_db > 0.05: 
        M_p = 10**(peaking_db / 20.0)
        if M_p > 1.0:
            inside_root = 1.0 - (1.0 / (M_p**2))
            inside_root = max(0.0, inside_root) 
            zeta_est = np.sqrt((1.0 - np.sqrt(inside_root)) / 2.0)
            pm_est = zeta_est * 100.0

    # Printouts to the console terminal
    print("\n" + "="*50)
    print(f"BODE ANALYSIS RESULTS:")
    print(f"  DC Gain:            {dc_gain:.3f} dB")
    
    if f_3db:
        f_unit = "MHz" if f_3db >= 1e6 else ("kHz" if f_3db >= 1e3 else "Hz")
        f_scale = 1e6 if f_3db >= 1e6 else (1e3 if f_3db >= 1e3 else 1.0)
        print(f"  -3dB Frequency:     {f_3db/f_scale:.3f} {f_unit}")
    
    print(f"  Peak Closed-Loop:   {max_mag:.3f} dB at {f_peak/1e3:.2f} kHz")
    print(f"  Relative Peaking:   {peaking_db:.3f} dB")
    
    if zeta_est is not None:
        print(f"  Estimated Damping (ζ): {zeta_est:.3f}")
        print(f"  Estimated PM (ϕm):     {pm_est:.1f}°")
    else:
        print("  Estimated Damping (ζ): > 0.707 (No significant peaking detected)")
    print("="*50 + "\n")

    # -------------------------------------------------------------------------
    # --- Plotting ---
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # Calculate magnitude layout breathing room safely ahead of execution
    mag_margin = (max(mag) - min(mag)) * 0.08 if max(mag) != min(mag) else 5

    # --- Magnitude Subplot ---
    ax1.semilogx(freq, mag, color=COLOR_MAG, linestyle='-', label="Measured")
    ax1.set_ylabel("Magnitude [dB]")
    ax1.set_title(plot_title)
    ax1.set_ylim(min(mag) - mag_margin, max(mag) + mag_margin) 
    
    # -3dB Line Overlay & Text Note Restored
    if f_3db:
        ax1.axhline(target_3db, color='#888888', linestyle='--', linewidth=0.6, alpha=0.5)
        ax1.axvline(f_3db, color='#888888', linestyle='--', linewidth=0.6, alpha=0.5)
        f_str = f"{f_3db/1e3:.1f} kHz" if f_3db < 1e6 else f"{f_3db/1e6:.2f} MHz"
        ax1.text(f_3db * 1.2, target_3db + (mag_margin * 0.1), f"$f_{{-3\\mathrm{{dB}}}}$ = {f_str}", fontsize=8, color='#555555')

    # Peaking Marker & Text Overlay Restored
    if peaking_db > 0.05:
        ax1.plot(f_peak, max_mag, 'o', color=COLOR_MAG, markersize=4)
        
        annot_text = (f"$M_p = {peaking_db:.1f}\\mathrm{{dB}}$\n"
                      f"$\\zeta \\approx {zeta_est:.2f}$\n"
                      f"$\\phi_m \\approx {pm_est:.0f}^\\circ$")
        
        ax1.text(f_peak * 1.3, max_mag - (mag_margin * 0.3), annot_text, 
                 fontsize=8, color='#333333', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # --- Phase Subplot ---
    ax2.semilogx(freq, phase, color=COLOR_PHASE, linestyle='-')
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Phase [deg]")
    
    if f_3db:
        ax2.axvline(f_3db, color='#888888', linestyle='--', linewidth=0.6, alpha=0.5)
    
    ax2.set_yticks([-180, -90, 0,45, 90,135, 180])
    ax2.set_ylim(-0, 190)

    # --- Structural Layout Styling (Dutch Academic Archetype) ---
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)   
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xlim(min(freq), max(freq))

    # --- Safe File Saving ---
    try:
        output_img = csv_filename.replace('.csv', '_thesis_bode.png')
        plt.savefig(output_img, bbox_inches='tight', dpi=300)
        print(f"File saved successfully as: {output_img}")
    except Exception as save_error:
        print(f"Skipping auto-save. (Reason: {save_error})")
    
    print("Opening plot window...")
    plt.show()

except Exception as e:
    print(f"Error reading or plotting data: {e}")