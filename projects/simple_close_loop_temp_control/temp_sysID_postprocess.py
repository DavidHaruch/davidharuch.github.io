import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Configuration ---
FILENAME = "step_test_1774681949.csv"
FS = 40.0
TS = 1.0 / FS
KP_USER, KI_USER = 0.5, 0.02 

# --- 1. Second-Order Model & Fit ---
def second_order_step(t, K, t1, t2):
    """Step response of 2 staggered RC stages:
    y(t) = K * (1 - (t1*exp(-t/t1) - t2*exp(-t/t2))/(t1 - t2))
    """
    t = np.maximum(t, 1e-6)
    # Prevent divide by zero if t1 == t2
    if abs(t1 - t2) < 1e-4: t2 = t1 * 1.001
    
    term1 = (t1 * np.exp(-t/t1)) / (t1 - t2)
    term2 = (t2 * np.exp(-t/t2)) / (t1 - t2)
    return K * (1 - (term1 - term2))

df = pd.read_csv(FILENAME)
t_data = df['time'].values - df['time'].values[0]
y_rel = (df['y'].values - df['d'].values) - (df['y'].values[0] - df['d'].values[0])
u_step = df['u'].values[0]

# Initial Guesses
y_ss = np.median(y_rel[-50:])
K_init = y_ss / u_step
t1_init = t_data[np.where(y_rel >= 0.63 * y_ss)[0][0]] * 0.8

popt, _ = curve_fit(
    lambda t, K, t1, t2: second_order_step(t, K * u_step, t1, t2),
    t_data, y_rel,
    p0=[K_init, t1_init, t1_init*0.2],
    bounds=([0, 0.1, 0.01], [np.inf, np.inf, np.inf])
)
K_f, t1_f, t2_f = popt

# --- 2. Figure 1: Time Domain Fit ---
plt.figure(figsize=(10, 4))
plt.plot(t_data, y_rel, 'k.', alpha=0.15, label='Measured Data')
plt.plot(t_data, second_order_step(t_data, K_f * u_step, t1_f, t2_f), 'r-', lw=2, label='2nd Order RC Fit')
plt.title(f"Figure 1: 2nd-Order Fit (K={K_f:.3f}, t1={t1_f:.1f}s, t2={t2_f:.1f}s)")
plt.xlabel("Time (s)"); plt.ylabel("Delta T (C)")
plt.legend(); plt.grid(True, alpha=0.3)
plt.show()

# --- 3. Frequency Domain (2nd Order + ZOH) ---
f = np.logspace(-4, np.log10(FS/2), 1000)
s = 1j * 2 * np.pi * f
G_plant = K_f / ((1 + t1_f*s) * (1 + t2_f*s))
G_zoh = np.sinc(f / FS) * np.exp(-s * TS / 2.0)
G = G_plant * G_zoh

# --- 4. Tuning & Metrics ---
def analyze(kp, ki):
    C = kp + ki/s
    L = G * C
    S = 1 / (1 + L)
    T = L / (1 + L)
    return L, S, T

# 45deg PM Tuning
idx_pm = np.argmin(np.abs(np.angle(G) - np.radians(-120)))
wc_opt = 2 * np.pi * f[idx_pm]
KP_OPT = 1.0 / np.abs(G[idx_pm])
KI_OPT = KP_OPT * (wc_opt / 10.0)

L_user, S_user, T_user = analyze(KP_USER, KI_USER)
L_opt, S_opt, T_opt = analyze(KP_OPT, KI_OPT)

# --- 5. Eight-Panel Visualization ---
fig, axs = plt.subplots(4, 2, figsize=(15, 16), sharex=True)
plt.subplots_adjust(hspace=0.25, wspace=0.2)

# Row 1: Plant Magnitude & Phase
axs[0,0].semilogx(f, 20*np.log10(np.abs(G_plant)), 'k')
axs[0,0].set_title("Plant Magnitude (dB)")
axs[0,1].semilogx(f, np.degrees(np.angle(G_plant)), 'k')
axs[0,1].set_title("Plant Phase (Max -180 deg)")

# Row 2: Open Loop (L) Magnitude & Phase
for L, lab, c in [(L_user, 'User', 'blue'), (L_opt, '45deg PM', 'green')]:
    axs[1,0].semilogx(f, 20*np.log10(np.abs(L)), color=c, label=lab)
    axs[1,1].semilogx(f, np.degrees(np.angle(L)), color=c)
axs[1,0].axhline(0, color='k', ls='--'); axs[1,1].axhline(-180, color='r', ls=':')
axs[1,0].set_title("Open Loop Gain (L)"); axs[1,0].legend()
axs[1,1].set_title("Open Loop Phase (L)")

# Row 3: Sensitivity (S) Magnitude & Phase
for S, c in [(S_user, 'blue'), (S_opt, 'green')]:
    axs[2,0].semilogx(f, 20*np.log10(np.abs(S)), color=c)
    axs[2,1].semilogx(f, np.degrees(np.angle(S)), color=c)
axs[2,0].axhline(6, color='orange', ls='--', label='6dB (Ms=2) Peak')
axs[2,0].set_title("Sensitivity Magnitude (S)")
axs[2,1].set_title("Sensitivity Phase (S)")

# Row 4: Closed Loop (T) Magnitude & Phase
for T, c in [(T_user, 'blue'), (T_opt, 'green')]:
    axs[3,0].semilogx(f, 20*np.log10(np.abs(T)), color=c)
    axs[3,1].semilogx(f, np.degrees(np.angle(T)), color=c)
axs[3,0].axhline(-3, color='k', ls=':'); axs[3,0].set_title("Closed Loop Magnitude (T)")
axs[3,1].set_title("Closed Loop Phase (T)")

for ax in axs.flat:
    ax.grid(True, which="both", alpha=0.2); ax.set_xlabel("Frequency (Hz)")

plt.show()

print(f"2nd Order Fit: K={K_f:.4f}, t1={t1_f:.2f}s, t2={t2_f:.2f}s")
print(f"45deg PM PI: Kp={KP_OPT:.4f}, Ki={KI_OPT:.6f}")