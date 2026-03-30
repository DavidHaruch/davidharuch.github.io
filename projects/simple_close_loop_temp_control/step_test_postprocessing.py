import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def analyze_step_run(filename):
    """ Loads CSV and computes step response metrics. """
    df = pd.read_csv(filename)
    
    # 1. Detect Step Input (where setpoint jumps)
    sp_diff = df['setpoint'].diff().abs()
    step_idx = sp_diff.idxmax()
    t_step = df.loc[step_idx, 'time_s']
    
    # Normalize time so the step occurs at t=0
    df['t_rel'] = df['time_s'] - t_step
    
    # Determine initial and final steady-state values
    # (Using averages to avoid noise-induced errors)
    y_initial = df[df['time_s'] < t_step]['temp_sec'].tail(40).mean()
    y_final = df['temp_sec'].tail(200).mean()
    step_height = y_final - y_initial
    
    # 2. Extract Metrics
    # Filter data to start from the step
    df_after = df[df['t_rel'] >= 0].copy()
    y = df_after['temp_sec'].values
    t = df_after['t_rel'].values
    
    d = df['temp_amb']
    
    # Peak Value & Overshoot
    peak_val = y.max()
    t_peak = t[y.argmax()]
    overshoot_pct = max(0, (peak_val - y_final) / step_height * 100)
    
    # jitter 120s to 300s
    y_jitter = np.max(y[4800:12000])-np.min(y[4800:12000])
    print(y_jitter)
    
    # Rise Time (10% to 90% of the step change)
    y_10 = y_initial + 0.1 * step_height
    y_90 = y_initial + 0.9 * step_height
    try:
        t_10 = t[np.where(y >= y_10)[0][0]]
        t_90 = t[np.where(y >= y_90)[0][0]]
        rise_time = t_90 - t_10
    except IndexError:
        rise_time = np.nan

    # Settling Time (2% error band)
    settle_band = 0.02 * step_height
    settle_band = 0.002
    within_band = np.abs(y - y_final) < settle_band
    # Find the last time the signal was OUTSIDE the band
    outside_indices = np.where(~within_band)[0]
    settling_time = t[outside_indices[-1]] if len(outside_indices) > 0 else 0
    
    # Capture control parameters
    metrics = {
        'File': filename,
        'Kp': df['kp'].iloc[0],
        'Ki': df['ki'].iloc[0],
        'Smith': "Yes" if df['is_smith'].iloc[0] == 1 else "No",
        'Rise Time (s)': round(rise_time, 2),
        'Settling Time (s)': round(settling_time, 2),
        'Overshoot (%)': round(overshoot_pct, 2),
        'Peak Temp (°C)': round(peak_val, 3),
        'Final Error (°C)': round(y_final - df['setpoint'].iloc[-1], 4),
        'Jitter (°C)': y_jitter
    }
    
    return df, metrics

# --- MAIN EXECUTION ---
# List your 3 file names here
# data_files = [r".\saved_data\step_response_1774814889_stdPI.csv", r".\saved_data\step_response_1774815566_smithpred.csv", r".\saved_data\step_response_1774828678_toohighgainPI.csv"]
# data_files = [r".\saved_data\step_response_1774815566_smithpred.csv", r".\saved_data\step_response_1774828678_toohighgainPI.csv"]
data_files = [r".\saved_data\step_response_1774814889_stdPI.csv", r".\saved_data\step_response_1774815566_smithpred.csv"]
data_files = [r".\saved_data\step_response_1774814889_stdPI.csv"]


results_table = []
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
fig.patch.set_facecolor('#F0F0F0')

for f in data_files:
    if not os.path.exists(f):
        print(f"Skipping {f}: File not found.")
        continue
        
    df, metrics = analyze_step_run(f)
    results_table.append(metrics)
    
    # Plotting
    label = f"{metrics['File']} (Kp={metrics['Kp']}, Smith={metrics['Smith']})"
    
    # label = f"{metrics['File']} (Kp={metrics['Kp']}, Smith={metrics['Smith']}, Pk2Pk Jitter={metrics['Jitter (°C)']})"
    
    ax1.plot(df['t_rel'], df['temp_sec'], label=label, linewidth=2)
    ax2.plot(df['t_rel'], df['duty'] * 100, alpha=0.7)

# Formatting Temperature Plot
ax1.axhline(df['setpoint'].iloc[-1], color='black', linestyle='--', alpha=0.5, label='Setpoint')
ax1.set_ylabel('Temperature (°C)', fontsize=12)
ax1.set_title('Step Response Comparison (Normalized Time)', fontsize=14)
ax1.set_xlim(0, 120) # Show 5s before and 120s after step

# ax1.set_xlim(300-120, 300) # Show 5s before and 120s after step
# ax1.set_ylim(31.03, 30.97) # Show 5s before and 120s after step

ax1.grid(True, which='both', alpha=0.3)
ax1.legend(loc='lower right')

# Formatting Duty Cycle Plot
ax2.set_ylabel('Heater Duty (%)', fontsize=12)
ax2.set_xlabel('Time from Step Change (s)', fontsize=12)
ax2.set_ylim(-5, 105)
ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig('step_response_analysis.png')

# Display Summary Table
summary_df = pd.DataFrame(results_table)
print("\n--- CONTROL METRICS SUMMARY ---")
print(summary_df.to_string(index=False))