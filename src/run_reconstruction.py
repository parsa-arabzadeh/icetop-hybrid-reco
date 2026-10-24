import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from icetop_array import SurfaceHybridArray

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 9.5,
    "figure.dpi": 300
})

class ShowerReconstructor:
    def __init__(self, array):
        self.array = array
        self.pos = array.coords
        self.c_ns = 0.299792458

    def fit_event(self, sig, times, thresh=0.3):
        mask = sig > thresh
        if np.sum(mask) < 5:
            return None

        # Center-of-gravity seed
        w = sig[mask] ** 0.5
        x_cog = np.sum(self.pos[mask, 0] * w) / np.sum(w)
        y_cog = np.sum(self.pos[mask, 1] * w) / np.sum(w)

        xh, yh = self.pos[mask, 0], self.pos[mask, 1]
        th, sh = times[mask], sig[mask]
        sig_t = 2.5 / np.sqrt(np.maximum(sh, 0.5))

        def chi2_timing(t0, nx, ny):
            dx = xh - x_cog
            dy = yh - y_cog
            proj = dx * nx + dy * ny
            r2 = np.maximum(dx**2 + dy**2 - proj**2, 1.0)
            t_exp = t0 + (proj / self.c_ns) + 1.5e-4 * r2
            pen = 1e6 * max(0.0, (nx**2 + ny**2) - 0.95)
            return np.sum(((th - t_exp) / sig_t) ** 2) + pen

        m_t = Minuit(chi2_timing, t0=np.median(th), nx=0.0, ny=0.0)
        m_t.limits["nx"] = (-0.85, 0.85)
        m_t.limits["ny"] = (-0.85, 0.85)
        m_t.migrad()

        nx_r, ny_r = m_t.values["nx"], m_t.values["ny"]
        n_mag = min(np.hypot(nx_r, ny_r), 0.98)
        theta_rec = np.degrees(np.arcsin(n_mag))
        phi_rec = (np.degrees(np.arctan2(-ny_r, -nx_r)) + 360.0) % 360.0

        # Combined LDF log-likelihood including non-triggered stations
        def nll_ldf(xc, yc, lg_s125):
            s125 = 10.0 ** lg_s125
            dx = self.pos[:, 0] - xc
            dy = self.pos[:, 1] - yc
            proj = dx * nx_r + dy * ny_r
            rp = np.sqrt(np.maximum(dx**2 + dy**2 - proj**2, 2.0))

            s_pred = self.array.ldf_nkg(rp, s125, age=1.25)
            var_h = np.maximum(s_pred[mask] * 0.1 + 0.05**2, 0.02)
            nll_hit = 0.5 * np.sum(((sig[mask] - s_pred[mask]) ** 2) / var_h + np.log(var_h))
            nll_nohit = np.sum(s_pred[~mask] * 2.0)
            return nll_hit + nll_nohit

        s_seed = np.clip(np.percentile(sh, 75), 1.0, 500.0)
        m_l = Minuit(nll_ldf, xc=x_cog, yc=y_cog, lg_s125=np.log10(s_seed))
        m_l.limits["xc"] = (-150.0, 150.0)
        m_l.limits["yc"] = (-150.0, 150.0)
        m_l.limits["lg_s125"] = (-0.5, 3.5)
        m_l.migrad()

        return {
            "xc": m_l.values["xc"],
            "yc": m_l.values["yc"],
            "s125": 10.0 ** m_l.values["lg_s125"],
            "theta": theta_rec,
            "phi": phi_rec,
            "ok": m_l.valid and m_t.valid
        }

    @staticmethod
    def space_angle(th1, ph1, th2, ph2):
        t1, p1 = np.radians(th1), np.radians(ph1)
        t2, p2 = np.radians(th2), np.radians(ph2)
        c = np.sin(t1) * np.sin(t2) * np.cos(p1 - p2) + np.cos(t1) * np.cos(t2)
        return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))

def run_benchmark(n_mc=250):
    os.makedirs("plots", exist_ok=True)
    arr = SurfaceHybridArray(nx=7, ny=7, spacing=50.0)
    reco = ShowerReconstructor(arr)
    rng = np.random.default_rng(42)

    dr_list, dpsi_list, s_ratio, s_true_list = [], [], [], []

    print(f"Reconstructing {n_mc} simulated showers with iminuit...")
    for _ in tqdm(range(n_mc)):
        xc_t = rng.uniform(-85.0, 85.0)
        yc_t = rng.uniform(-85.0, 85.0)
        s125_t = 10.0 ** rng.uniform(0.5, 2.3)
        th_t = rng.uniform(5.0, 35.0)
        ph_t = rng.uniform(0.0, 360.0)

        sig, times, _, _ = arr.simulate_event((xc_t, yc_t), s125_t, th_t, ph_t, rng=rng)
        fit = reco.fit_event(sig, times)
        if fit and fit["ok"]:
            dr_list.append(np.hypot(fit["xc"] - xc_t, fit["yc"] - yc_t))
            dpsi_list.append(reco.space_angle(th_t, ph_t, fit["theta"], fit["phi"]))
            s_ratio.append(fit["s125"] / s125_t)
            s_true_list.append(s125_t)

    dr_arr = np.array(dr_list)
    dpsi_arr = np.array(dpsi_list)
    s_rat = np.array(s_ratio)
    s_tr = np.array(s_true_list)

    r68 = np.percentile(dr_arr, 68)
    psi68 = np.percentile(dpsi_arr, 68)
    print(f"Core resolution (68%): {r68:.2f} m | Angular PSF (68%): {psi68:.2f} deg")

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.5, 4.5))

    ax1.hist(dr_arr, bins=24, range=(0, 25), color="#1f77b4", alpha=0.75, edgecolor="k")
    ax1.axvline(r68, color="crimson", ls="--", lw=1.8, label=f"68%: {r68:.1f} m")
    ax1.set_xlabel(r"$\Delta R_{\mathrm{core}}$ [m]")
    ax1.set_ylabel("Events")
    ax1.set_title("Core Position Resolution")
    ax1.grid(True, ls=":", alpha=0.5)
    ax1.legend()

    ax2.hist(dpsi_arr, bins=24, range=(0, 4.0), color="#2ca02c", alpha=0.75, edgecolor="k")
    ax2.axvline(psi68, color="crimson", ls="--", lw=1.8, label=f"$\\Psi_{{68}} = {psi68:.2f}^\\circ$")
    ax2.set_xlabel(r"Opening Angle $\Delta \Psi$ [deg]")
    ax2.set_ylabel("Events")
    ax2.set_title("Angular Resolution")
    ax2.grid(True, ls=":", alpha=0.5)
    ax2.legend()

    ax3.scatter(s_tr, s_rat, c="#ff7f0e", alpha=0.7, edgecolor="k", s=26)
    ax3.axhline(1.0, color="k", ls="--", lw=1.3)
    ax3.set_xscale("log")
    ax3.set_ylim(0.5, 1.5)
    ax3.set_xlabel(r"$S_{125}^{\mathrm{true}}$ [VEM]")
    ax3.set_ylabel(r"$S_{125}^{\mathrm{rec}} / S_{125}^{\mathrm{true}}$")
    ax3.set_title(r"$S_{125}$ Reconstruction Linearity")
    ax3.grid(True, which="both", ls=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig("plots/reconstruction_performance.png")
    plt.close()

if __name__ == "__main__":
    run_benchmark()
