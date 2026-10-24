import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 9.5,
    "figure.dpi": 300
})

class SurfaceHybridArray:
    """
    Toy geometry and response model for an extended surface array of
    scintillation detectors and radio antennas at the South Pole.
    """
    def __init__(self, nx=7, ny=7, spacing=50.0):
        self.spacing = spacing
        gx = np.linspace(-(nx // 2) * spacing, (nx // 2) * spacing, nx)
        gy = np.linspace(-(ny // 2) * spacing, (ny // 2) * spacing, ny)
        xx, yy = np.meshgrid(gx, gy)
        self.coords = np.column_stack([xx.ravel(), yy.ravel()])
        self.n_det = len(self.coords)

        # Unit vector in local coordinate system
        b_inc = np.radians(-72.0)
        self.b_hat = np.array([0.0, np.cos(b_inc), np.sin(b_inc)])

    @staticmethod
    def ldf_nkg(r_m, s125, age=1.25, r_mol=125.0):
        r = np.maximum(r_m, 2.0)
        x = r / r_mol
        x_ref = 125.0 / r_mol
        val = (x ** (age - 2.0)) * ((1.0 + x) ** (age - 4.5))
        ref = (x_ref ** (age - 2.0)) * ((1.0 + x_ref) ** (age - 4.5))
        return s125 * (val / ref)

    def simulate_event(self, core_xy, s125, theta_deg, phi_deg, age=1.25, rng=None):
        if rng is None:
            rng = np.random.default_rng()

        th = np.radians(theta_deg)
        ph = np.radians(phi_deg)
        
        # Shower propagation unit vector
        n_vec = np.array([-np.sin(th) * np.cos(ph), -np.sin(th) * np.sin(ph), -np.cos(th)])
        nx, ny = n_vec[0], n_vec[1]

        dx = self.coords[:, 0] - core_xy[0]
        dy = self.coords[:, 1] - core_xy[1]

        proj = dx * nx + dy * ny
        r_perp = np.sqrt(np.maximum(dx**2 + dy**2 - proj**2, 1.0))

        s_exp = self.ldf_nkg(r_perp, s125, age=age)

        # Scintillator response: photoelectron statistics + baseline noise
        pe = rng.poisson(np.maximum(s_exp * 10.0, 0.0)) / 10.0
        noise = rng.normal(0.0, 0.15, size=self.n_det)
        s_meas = np.where(pe > 0.3, np.maximum(pe + noise, 0.0), 0.0)

        # Shower front arrival time (ns)
        c_ns = 0.299792458
        t_curv = 1.5e-4 * (r_perp ** 2)
        sigma_t = 2.5 / np.sqrt(np.maximum(s_meas, 0.5))
        t_meas = 1000.0 + (proj / c_ns) + t_curv + rng.normal(0.0, sigma_t)

        # Geomagnetic radio amplitude proxy: E_radio ~ E_em * |v x B| * exp(-r / r0)
        sin_alpha = np.linalg.norm(np.cross(n_vec, self.b_hat))
        radio_amp = 12.0 * s125 * sin_alpha * np.exp(-r_perp / 110.0)

        return s_meas, t_meas, r_perp, radio_amp

if __name__ == "__main__":
    os.makedirs("plots", exist_ok=True)
    arr = SurfaceHybridArray(nx=7, ny=7, spacing=50.0)
    rng = np.random.default_rng(101)

    core_true = (18.0, -22.0)
    s125_true = 25.0
    th_true, ph_true = 22.0, 135.0

    sig, t_arr, r_p, _ = arr.simulate_event(core_true, s125_true, th_true, ph_true, rng=rng)
    hits = sig > 0.3

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))

    ax1.scatter(arr.coords[~hits, 0], arr.coords[~hits, 1], c="lightgray", s=30, marker="s", label="Sub-threshold (< 0.3 MIP)")
    sc = ax1.scatter(arr.coords[hits, 0], arr.coords[hits, 1], c=t_arr[hits],
                     s=np.clip(np.log10(sig[hits] + 1.0) * 130, 35, 420),
                     cmap="plasma", edgecolor="k", lw=0.7, label="Triggered Stations")
    ax1.scatter(core_true[0], core_true[1], c="cyan", edgecolor="k", marker="*", s=240, zorder=5, label="True Core")

    cb = plt.colorbar(sc, ax=ax1)
    cb.set_label("Arrival Time [ns]")
    ax1.set_xlabel("X [m]")
    ax1.set_ylabel("Y [m]")
    ax1.set_title(f"Simulated Event Footprint ($S_{{125}}={s125_true}$ VEM, $\\theta={th_true}^\\circ$)")
    ax1.grid(True, ls=":", alpha=0.5)
    ax1.legend(loc="upper right", fontsize=9)

    r_grid = np.logspace(0.5, 2.5, 200)
    ax2.loglog(r_grid, arr.ldf_nkg(r_grid, s125_true), "r-", lw=2.0, label="NKG Model")
    ax2.errorbar(r_p[hits], sig[hits], yerr=np.sqrt(sig[hits] * 0.1 + 0.04),
                 fmt="ko", ms=4.5, capsize=2.5, alpha=0.8, label="Station Signals")
    ax2.axvline(125.0, color="royalblue", ls="--", lw=1.4, label=r"$R_{\mathrm{ref}} = 125\ \mathrm{m}$")

    ax2.set_xlabel(r"Shower Axis Distance $r_\perp$ [m]")
    ax2.set_ylabel("Signal [VEM]")
    ax2.set_title("Lateral Distribution Profile")
    ax2.set_xlim(5, 300)
    ax2.set_ylim(0.1, 1e3)
    ax2.grid(True, which="both", ls=":", alpha=0.5)
    ax2.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig("plots/array_event_footprint.png")
    plt.close()
    print("Saved plots/array_event_footprint.png")
