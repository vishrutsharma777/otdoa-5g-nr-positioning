"""Draw the README figures from the saved Monte Carlo results (EE666_results_N50.mat).

Usage:  python make_figures.py      (needs numpy, scipy, matplotlib)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
m = sio.loadmat(os.path.join(HERE, "EE666_results_N50.mat"), squeeze_me=True)
out = os.path.join(HERE, "images")
os.makedirs(out, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
N = int(m["N_MC"])

# 1. ToA estimator accuracy versus SNR (single link, 100 MHz configuration)
fig, ax = plt.subplots(figsize=(7.5, 4.6))
for name, rms, mk in zip(m["methods"], m["rms_est"], ["o-", "s-", "^-"]):
    ax.semilogy(m["snr_est"], rms, mk, ms=4, lw=1.6, label=str(name))
ax.axhline(1 / np.sqrt(12), color="grey", ls="--", lw=1)
ax.text(m["snr_est"][-1], 1 / np.sqrt(12) * 1.15, "integer quantisation 1/√12", ha="right", fontsize=9, color="grey")
ax.set_xlabel("Per-sample SNR (dB)")
ax.set_ylabel("RMS ToA error (samples)")
ax.set_title("ToA estimator accuracy versus SNR (one link, 100 MHz)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(out, "toa_estimators.png"), dpi=150)

# 2. Median positioning error versus SNR, per bandwidth and channel
snr = m["SNR_dB_MC"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
for b, (ax, bw) in enumerate(zip(axes, ["20 MHz", "100 MHz"])):
    for c, (ch, mk) in enumerate(zip(["AWGN", "TDL-A", "TDL-C"], ["o-", "s-", "^-"])):
        med = [np.median(m["errors"][b, c, s]) for s in range(len(snr))]
        ax.semilogy(range(len(snr)), med, mk, ms=5, lw=1.6, label=ch)
    ax.set_xticks(range(len(snr)), [str(int(s)) for s in snr])
    ax.set_xlabel("SNR (dB)")
    ax.set_title(f"PRS bandwidth {bw}")
axes[0].set_ylabel("Median positioning error (m)")
axes[1].legend()
fig.suptitle(f"OTDOA positioning error, band-limited ToA estimator (N = {N} UEs per point)")
fig.tight_layout()
fig.savefig(os.path.join(out, "positioning_vs_snr.png"), dpi=150)

# 3. Doppler sensitivity (TDL-A)
fd = m["doppler_vals"]
fig, ax = plt.subplots(figsize=(7.5, 4.4))
for b, (bw, mk, col) in enumerate(zip(["20 MHz", "100 MHz"], ["o-", "s-"], ["tab:blue", "tab:red"])):
    p50 = [np.median(e) for e in m["errors_dop"][b]]
    p90 = [np.percentile(e, 90) for e in m["errors_dop"][b]]
    ax.plot(range(len(fd)), p50, mk, color=col, lw=1.8, label=f"{bw} median")
    ax.plot(range(len(fd)), p90, mk[0] + "--", color=col, lw=1, alpha=0.5, label=f"{bw} 90th percentile")
ax.set_xticks(range(len(fd)), [f"{int(f)} Hz" for f in fd])
ax.set_xlabel("Maximum Doppler shift")
ax.set_ylabel("Positioning error (m)")
ax.set_title(f"Effect of UE mobility, TDL-A channel (N = {N})")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(out, "doppler.png"), dpi=150)
print("figures written to", out)
