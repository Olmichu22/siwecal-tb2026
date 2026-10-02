"""MIP intercalibration of the fixed tables: channel MPV distribution and per-slab medians, three threshold sets."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi_fixed/mips"
F = {"th210": f"{D}/th210/MIP_pedestalsubmode1_TB2026CERN_run_000th210_highgain.txt",
     "th220": f"{D}/th220/MIP_pedestalsubmode1_TB2026CERN_run_000th220_highgain.txt",
     "th230": f"{D}/th230/MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt"}
COL = {"th210": "#2B6CB0", "th220": "#B5541A", "th230": "#5E8C4A"}
plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})
for th, f in F.items():
    d = np.loadtxt(f, comments="#")
    fit = d[:, 4] >= 0
    std = fit & ~np.isin(d[:, 0], [12, 14])
    q = np.percentile(d[std, 3], [16, 50, 84]); spread = (q[2] - q[0]) / 2 / q[1]
    a1.hist(d[std, 3] / q[1], bins=np.linspace(0.6, 1.4, 81), histtype="step", lw=1.8, color=COL[th],
            label=f"{th}: {std.sum()} channels, MPV {q[1]:.1f} ADC, spread {100*spread:.1f} %")
    med = [np.median(d[fit & (d[:, 0] == s), 3]) if (fit & (d[:, 0] == s)).sum() >= 10 else np.nan for s in range(15)]
    a2.plot(range(15), np.array(med) / q[1], "o-", color=COL[th], ms=6, lw=1.5, label=th)
a1.set_xlabel("channel MIP MPV / median of the set"); a1.set_ylabel("channels")
a1.set_title("Channel-to-channel (FEV10 slabs, 500 µm)", fontsize=13)
a1.legend(frameon=False, fontsize=10, loc="upper left"); a1.set_ylim(0, a1.get_ylim()[1] * 1.3)
a2.axhline(1, color="#999", lw=1, ls="--")
a2.set_xlabel("slab"); a2.set_ylabel("slab median MPV / median of the set"); a2.set_xticks(range(15))
a2.set_title("Slab to slab", fontsize=13)
a2.annotate("slab 12: own threshold (DAC 243)", (12, 1.85), (6.2, 1.75), fontsize=10, arrowprops=dict(arrowstyle="->", color="#555"))
a2.annotate("slab 14: 650 µm sensor", (14, 1.27), (8.0, 1.5), fontsize=10, arrowprops=dict(arrowstyle="->", color="#555"))
a2.legend(frameon=False, fontsize=11, loc="upper left")
fig.suptitle("MIP intercalibration, high gain (MuonCalib_gaudi_fixed)", fontsize=14)
fig.tight_layout(); fig.savefig(f"{__file__.rsplit('/',1)[0]}/intercalibration.png", dpi=160)
