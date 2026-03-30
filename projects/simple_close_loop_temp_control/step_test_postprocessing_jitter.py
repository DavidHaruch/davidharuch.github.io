import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def analyze_jitter_and_rejection(filename):
    """ Loads CSV and computes jitter/rejection metrics between 120s and 300s. """
    df = pd.read_csv(filename)
    
    # 1. Normalize time based on step input
    sp_diff = df['setpoint'].diff().abs()
    step_idx = sp_diff.idxmax()
    t_step = df.loc[step_idx, 'time_s']
    df['t_rel'] = df['time_s'] - t_step
    
    # 2. Slice the steady-state window (120s to 300s)
    mask = (df['t_rel'] >= 140) & (df['t_rel'] <= 300)
    df_window = df[mask].copy()
    
    if df_window.empty:
        return None, None

    # 3. Compute Peak-to-Peak (Pk2Pk) values
    y_pk2pk = df_window['temp_sec'].max() - df_window['temp_sec'].min()
    amb_pk2pk = df_window['temp_amb'].max() - df_window['temp_amb'].min()
    
    # 4. Compute Reduction Ratio (Disturbance Rejection)
    # Ratio = Ambient Swing / Feedback Swing
    # Example: If Amb swings 100mC and Feedback only 2mC, Ratio = 50.
    reduction_ratio = amb_pk2pk / y_pk2pk if y_pk2pk > 0 else np.inf
    
    metrics = {
        'File': os.path.basename(filename),
        'Smith': "Yes" if df['is_smith'].iloc[0] == 1 else "No",
        'Kp': df['kp'].iloc[0],
        'Feedback Pk2Pk (mC)': round(y_pk2pk * 1000, 2),
        'Ambient Pk2Pk (mC)': round(amb_pk2pk * 1000, 2),
        'Reduction Ratio': round(reduction_ratio, 2)
    }
    
    return df_window, metrics

# --- DATA FILES ---
data_files = [
    r".\saved_data\step_response_1774814889_stdPI.csv", 
    r".\saved_data\step_response_1774815566_smithpred.csv"
]

results = []
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
fig.patch.set_facecolor('#F8F9FA')

for f in data_files:
    if not os.path.exists(f): continue
        
    df_win, m = analyze_jitter_and_rejection(f)
    if df_win is not None:
        results.append(m)
        
        # Plot 1: Feedback Jitter (Secondary Sensor)
        lbl = f"{m['Smith']} Smith | Pk2Pk: {m['Feedback Pk2Pk (mC)']}(mC) | Ratio: {m['Reduction Ratio']}x"
        ax1.plot(df_win['t_rel'], df_win['temp_sec'], label=lbl, linewidth=1.5)
        
        # Plot 2: Ambient Jitter
        ax2.plot(df_win['t_rel'], df_win['temp_amb'], label=f"Amb ({m['File']}) | | Pk2Pk: {m['Ambient Pk2Pk (mC)']}(mC)")

# Formatting
ax1.set_ylabel('Feedback Temp AIN2 (°C)', fontsize=12)
ax1.set_title('Steady-State Jitter & Disturbance Rejection (120s - 300s)', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right', fontsize=9)

ax2.set_ylabel('Ambient Temp AIN0 (°C)', fontsize=12)
ax2.set_xlabel('Time from Step (s)', fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right', fontsize=8)

plt.tight_layout()
plt.show()

# Display Results Table
summary_df = pd.DataFrame(results)
print("\n--- THERMAL ISOLATION SUMMARY ---")
cols = ['File', 'Smith', 'Feedback Pk2Pk (mC)', 'Ambient Pk2Pk (mC)', 'Reduction Ratio']
print(summary_df[cols].to_string(index=False))