import os
import re
import glob
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D

ffmpeg_executable = r"C:\Users\david\AppData\Local\Microsoft\WinGet\Packages\Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe\ffmpeg-9.0.2-full_build\bin\ffmpeg.exe"
mpl.rcParams['animation.ffmpeg_path'] = ffmpeg_executable

class ODSVideoExporter:
    def __init__(self, frf_dir, nodes_csv, lines_csv=None):
        self.frf_dir = frf_dir
        self.nodes_csv = nodes_csv
        self.lines_csv = lines_csv

        self.nodes_df = None
        self.lines_list = []
        self.frf_data = {}
        self.frequencies = None

        self._load_data()

    def _load_data(self):
        # 1. Load Node Coordinates
        self.nodes_df = pd.read_csv(self.nodes_csv)
        self.nodes_df.columns = [c.strip().upper() for c in self.nodes_df.columns]

        # 2. Load Line Connectivity
        if self.lines_csv and os.path.exists(self.lines_csv):
            lines_df = pd.read_csv(self.lines_csv)
            lines_df.columns = [c.strip().upper() for c in lines_df.columns]
            self.lines_list = [(int(row.iloc[0]), int(row.iloc[1])) for _, row in lines_df.iterrows()]
        else:
            # Fallback: Consecutive lines along node order
            nodes = sorted(self.nodes_df['NODE'].tolist())
            self.lines_list = [(nodes[i], nodes[i+1]) for i in range(len(nodes)-1)]

        # 3. Load FRF Data
        csv_files = glob.glob(os.path.join(self.frf_dir, "*.csv"))
        for file_path in csv_files:
            fname = os.path.basename(file_path)
            match = re.search(r'Accel_Node(\d+)', fname, re.IGNORECASE)
            if match:
                node_id = int(match.group(1))
                df = pd.read_csv(file_path, comment='#')
                df.columns = [c.strip() for c in df.columns]
                self.frf_data[node_id] = df

        if self.frf_data:
            sample_node = list(self.frf_data.keys())[0]
            self.frequencies = self.frf_data[sample_node].iloc[:, 0].values

    def export_video(self, target_freq_hz, output_filename, anim_gain=15.0, fps=30, dpi=150, elev=30, azim=45):
        """
        Exports a seamless MP4 or WebM video loop for a specific frequency.
        """
        if self.frequencies is None:
            raise ValueError("No FRF data loaded.")

        # Find closest frequency index
        freq_idx = (np.abs(self.frequencies - target_freq_hz)).argmin()
        actual_freq = self.frequencies[freq_idx]
        print(f"Exporting ODS animation for {actual_freq:.2f} Hz...")

        # Extract node geometry
        node_coords = {int(row['NODE']): (row['X'], row['Y'], row['Z']) for _, row in self.nodes_df.iterrows()}
        node_ids = np.array(list(node_coords.keys()))
        x_base = np.array([node_coords[n][0] for n in node_ids])
        y_base = np.array([node_coords[n][1] for n in node_ids])
        z_base = np.array([node_coords[n][2] for n in node_ids])

        # Calculate peak magnitude across nodes for scaling normalization
        mags = []
        for nid in node_ids:
            if nid in self.frf_data:
                df = self.frf_data[nid]
                r_col = [c for c in df.columns if 'Real' in c][0]
                i_col = [c for c in df.columns if 'Imag' in c][0]
                r = df[r_col].iloc[freq_idx]
                i = df[i_col].iloc[freq_idx]
                mags.append(np.sqrt(r**2 + i**2))
        peak_mag = np.max(mags) if len(mags) > 0 and np.max(mags) > 1e-12 else 1.0

        # Create Matplotlib Figure
        fig = plt.figure(figsize=(8, 8), dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')

        # Compute bounding limits for fixed 1:1 spatial aspect ratio
        max_span = max(
            np.ptp(x_base),
            np.ptp(y_base),
            max(anim_gain * 2.0, np.ptp(z_base))
        )
        if max_span == 0:
            max_span = 1.0

        x_mid = (np.max(x_base) + np.min(x_base)) / 2.0
        y_mid = (np.max(y_base) + np.min(y_base)) / 2.0
        z_mid = (np.max(z_base) + np.min(z_base)) / 2.0
        half_span = max_span * 0.55

        # Render function for each frame (30 frames total = 12 deg phase step per frame)
        num_frames = 30
        phase_angles = np.linspace(0, 360, num_frames, endpoint=False)

        def update_frame(frame_idx):
            ax.clear()
            phase_deg = phase_angles[frame_idx]
            phase_rad = np.radians(phase_deg)

            # Compute current phase displacement
            z_disp = np.zeros(len(node_ids))
            for i, nid in enumerate(node_ids):
                if nid in self.frf_data:
                    df = self.frf_data[nid]
                    r_col = [c for c in df.columns if 'Real' in c][0]
                    i_col = [c for c in df.columns if 'Imag' in c][0]
                    r = df[r_col].iloc[freq_idx]
                    im = df[i_col].iloc[freq_idx]
                    complex_val = complex(r, im)
                    z_disp[i] = (complex_val * np.exp(1j * phase_rad)).real

            scaled_z_disp = (z_disp / peak_mag) * anim_gain
            z_deformed = z_base + scaled_z_disp
            deformed_coords = {nid: (x, y, zd) for nid, x, y, zd in zip(node_ids, x_base, y_base, z_deformed)}

            # Plot Undeformed Wireframe
            for na, nb in self.lines_list:
                if na in node_coords and nb in node_coords:
                    xa, ya, za = node_coords[na]
                    xb, yb, zb = node_coords[nb]
                    ax.plot([xa, xb], [ya, yb], [za, zb], 'k--', alpha=0.25)

            # Plot Deformed Wireframe
            for na, nb in self.lines_list:
                if na in deformed_coords and nb in deformed_coords:
                    xa, ya, za = deformed_coords[na]
                    xb, yb, zb = deformed_coords[nb]
                    ax.plot([xa, xb], [ya, yb], [za, zb], 'r-', linewidth=2.5)

            # Plot Scatter Nodes
            ax.scatter(x_base, y_base, z_base, color='gray', s=15, alpha=0.3)
            ax.scatter(x_base, y_base, z_deformed, color='red', s=25)

            for nid, x, y, zd in zip(node_ids, x_base, y_base, z_deformed):
                ax.text(x, y, zd, f' {nid}', color='blue', fontsize=8)

            # Axis limits & 1:1 scaling
            ax.set_xlim(x_mid - half_span, x_mid + half_span)
            ax.set_ylim(y_mid - half_span, y_mid + half_span)
            ax.set_zlim(z_mid - half_span, z_mid + half_span)
            ax.set_box_aspect((1, 1, 1))

            ax.set_title(f"ODS @ {actual_freq:.0f} Hz", fontsize=12, pad=10)
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.view_init(elev=elev, azim=azim)

        # Build Animation Object
        anim = FuncAnimation(fig, update_frame, frames=num_frames, interval=1000/fps)

        # Select FFmpeg Codec based on output extension
        ext = os.path.splitext(output_filename)[1].lower()
        if ext == '.webm':
            writer = FFMpegWriter(fps=fps, codec='libvpx-vp9', bitrate=2000)
        else:
            # Default MP4 (H.264)
            writer = FFMpegWriter(fps=fps, codec='libx264', extra_args=['-pix_fmt', 'yuv420p'])

        anim.save(output_filename, writer=writer)
        plt.close(fig)
        print(f"Successfully saved animation to {output_filename}")


if __name__ == "__main__":
    # Ensure FFmpeg is installed and accessible in your system PATH
    exporter = ODSVideoExporter(
        frf_dir="./2026_09_27_tecnotionux9s_freefree_10pts",
        nodes_csv="./tecnotion_nodes.csv",
        lines_csv="./tecnotion_lines.csv"
    )

    # Export MP4 Loop (H.264)
    # exporter.export_video(
    #     target_freq_hz=240.0,
    #     output_filename="ods_240hz.mp4",
    #     anim_gain=50.0,
    #     fps=30
    # )
    # exporter.export_video(
    #     target_freq_hz=314.0,
    #     output_filename="ods_314hz.mp4",
    #     anim_gain=50.0,
    #     fps=30
    # )
    # exporter.export_video(
    #     target_freq_hz=542.0,
    #     output_filename="ods_542hz.mp4",
    #     anim_gain=70.0,
    #     fps=30
    # )
    # exporter.export_video(
    #     target_freq_hz=686.0,
    #     output_filename="ods_686hz.mp4",
    #     anim_gain=35.0,
    #     fps=30
    # )
    # exporter.export_video(
    #     target_freq_hz=914.0,
    #     output_filename="ods_914hz.mp4",
    #     anim_gain=35.0,
    #     fps=30
    # )
    exporter.export_video(
        target_freq_hz=1232.0,
        output_filename="ods_1232hz.mp4",
        anim_gain=35.0,
        fps=30
    )

    # # Export WebM Loop (VP9)
    # exporter.export_video(
    #     target_freq_hz=245.0,
    #     output_filename="ods_245hz.webm",
    #     anim_gain=15.0,
    #     fps=30
    # )