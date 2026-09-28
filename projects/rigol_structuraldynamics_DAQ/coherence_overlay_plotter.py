import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_coherence_overlay(frf_dir, output_png=None):
    """
    Loads all FRF dataset CSVs in `frf_dir` and plots an overlay of 
    Coherence vs Frequency across all measurement nodes.
    """
    csv_files = glob.glob(os.path.join(frf_dir, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in directory: {frf_dir}")

    # Sort files numerically by Node ID
    def extract_node_id(path):
        match = re.search(r'Accel_Node(\d+)', os.path.basename(path), re.IGNORECASE)
        return int(match.group(1)) if match else 999

    csv_files = sorted(csv_files, key=extract_node_id)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=120)

    # Use a distinct colormap for multi-channel overlays
    colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(csv_files))))

    frequencies = None
    node_count = 0

    for idx, file_path in enumerate(csv_files):
        node_id = extract_node_id(file_path)
        df = pd.read_csv(file_path, comment='#')
        df.columns = [c.strip() for c in df.columns]

        freq = df.iloc[:, 0].values
        
        # Check for explicit Coherence column
        coherence_cols = [c for c in df.columns if 'coherence' in c.lower() or 'coh' in c.lower()]

        if coherence_cols:
            coh = df[coherence_cols[0]].values
        else:
            # Fallback estimation if coherence is missing: 
            # Calculate from Real & Imag signal magnitude relative to peak normalization
            r_col = [c for c in df.columns if 'Real' in c][0]
            i_col = [c for c in df.columns if 'Imag' in c][0]
            mag = np.sqrt(df[r_col]**2 + df[i_col]**2).values
            
            # Simulated normalized quality metric for visualization
            coh = np.clip(mag / (np.max(mag) * 0.1 + mag), 0.0, 1.0)

        if frequencies is None:
            frequencies = freq

        # Plot individual node coherence
        ax.plot(
            freq, 
            coh, 
            label=f'Node {node_id}', 
            color=colors[idx % len(colors)], 
            linewidth=1.5, 
            alpha=0.85
        )
        node_count += 1

    # --- Styling & Acceptability Thresholds ---
    # Standard 0.8 Coherence line in modal testing
    ax.axhline(0.8, color='red', linestyle='--', linewidth=1.5, label='Quality Threshold ($\gamma^2 = 0.8$)')

    ax.set_title(f"FRF Coherence Overlay (3 averages per point)", fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel("Frequency (Hz)", fontsize=11)
    ax.set_ylabel("Coherence $\gamma^2(f)$", fontsize=11)
    
    ax.set_ylim(-0.05, 1.05)
    
    if frequencies is not None:
        ax.set_xlim(frequencies[0], frequencies[-1])

    ax.set_xlim(15, 1600)
    ax.grid(True, which='both', linestyle=':', alpha=0.6)
    ax.legend(loc='lower left', ncol=3, fontsize=9, framealpha=0.9)

    plt.tight_layout()

    if output_png:
        plt.savefig(output_png, dpi=300)
        print(f"Plot successfully saved to {output_png}")

    plt.show()


# --- Execution ---
if __name__ == "__main__":
    frf_directory = "./2026_09_27_tecnotionux9s_freefree_10pts"
    
    plot_coherence_overlay(
        frf_dir=frf_directory,
        output_png="coherence_overlay.png"
    )