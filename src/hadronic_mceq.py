import os
import numpy as np
import matplotlib.pyplot as plt
from MCEq.core import MCEqRun
import crflux.models as crf

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 9.5,
    "figure.dpi": 300
})

def estimate_surface_ne(e0_gev, a_mass, x_obs, model_name):
    # Parameterization of the EM cascade size at South Pole depth (X ~ 705 g/cm^2)
    # using Gaisser-Hillas profile anchored to hadronic model elongation rates
    x0_em = 36.7  # radiation length in air [g/cm^2]
    e_per_nucleus = e0_gev / a_mass
    
    # Model-dependent X_max shift [g/cm^2]
    dx_model = {"SIBYLL23D": 0.0, "EPOSLHC": 6.5, "QGSJETII04": -11.0}[model_name]
    
    # Average X_max for proton/iron at given energy per nucleon
    x_max = 580.0 + 58.0 * np.log10(e_per_nucleus / 1e6) + dx_model
    x_max = np.clip(x_max, 420.0, 850.0)
    
    # Gaisser-Hillas longitudinal profile evaluated at IceTop surface depth x_obs
    n_max = e0_gev / 1.4  # ~1.4 GeV per charged EM particle at shower maximum
    lam = 70.0
    rel_x = np.maximum((x_obs - 100.0) / (x_max - 100.0), 0.05)
    profile = (rel_x ** ((x_max - 100.0) / lam)) * np.exp((x_max - x_obs) / lam)
    return n_max * profile

def solve_southpole_showers():
    h_obs_cm = 2835.0 * 100.0
    e0_grid_gev = np.logspace(5.5, 8.0, 11)  # ~300 TeV to 100 PeV
    
    models = ["SIBYLL23D", "EPOSLHC", "QGSJETII04"]
    labels = {
        "SIBYLL23D": "Sibyll 2.3d",
        "EPOSLHC": "EPOS-LHC",
        "QGSJETII04": "QGSJet-II-04"
    }

    out = {m: {"p": {"Ne": [], "Nmu": []}, "Fe": {"Ne": [], "Nmu": []}} for m in models}

    mceq = MCEqRun(
        interaction_model="SIBYLL23D",
        primary_model=(crf.HillasGaisser2012, "H3a"),
        theta_deg=0.0,
        density_model=("MSIS00_IC", ("SouthPole", "December"))
    )

    x_obs = mceq.density_model.h2X(h_obs_cm)
    print(f"Observation depth at h=2835m: X = {x_obs:.1f} g/cm2")

    e_grid = mceq.e_grid
    de = mceq.e_widths
    mask_mu_ice = e_grid >= 500.0

    for m in models:
        print(f"Running MCEq for {labels[m]}...")
        mceq.set_interaction_model(m)
        
        for p_key, cid, a_mass in [("p", 14, 1), ("Fe", 5626, 56)]:
            for e0 in e0_grid_gev:
                mceq.set_single_primary_particle(E=e0, corsika_id=cid)
                mceq.solve(int_grid=[x_obs])
                
                phi_mu = mceq.get_solution("total_mu-", grid_idx=0) + mceq.get_solution("total_mu+", grid_idx=0)
                nmu = np.sum((phi_mu * de)[mask_mu_ice])
                ne = estimate_surface_ne(e0, a_mass, x_obs, m)
                
                out[m][p_key]["Ne"].append(ne)
                out[m][p_key]["Nmu"].append(nmu)

    os.makedirs("plots", exist_ok=True)
    os.makedirs("data", exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    cols = {"SIBYLL23D": "#1f77b4", "EPOSLHC": "#d62728", "QGSJETII04": "#2ca02c"}
    ls_map = {"p": "-", "Fe": "--"}

    for m in models:
        for p_key, p_name in [("p", "Proton"), ("Fe", "Iron")]:
            ne_vals = np.array(out[m][p_key]["Ne"])
            nmu_vals = np.array(out[m][p_key]["Nmu"])
            lbl = f"{labels[m]} ({p_name})"
            ax1.loglog(ne_vals, nmu_vals, color=cols[m], ls=ls_map[p_key], lw=1.8, marker="o", ms=3.5, label=lbl)

    ax1.set_xlabel(r"Surface EM Shower Size $N_e$ ($X = 705\ \mathrm{g/cm^2}$)")
    ax1.set_ylabel(r"In-Ice Muon Yield $N_\mu$ ($E_\mu > 500\ \mathrm{GeV}$)")
    ax1.set_title("Surface vs. In-Ice Coincidence (Mass Separation)")
    ax1.grid(True, which="both", ls=":", alpha=0.5)
    ax1.legend(loc="upper left", fontsize=8.5, ncol=2)

    e0_ev = e0_grid_gev * 1e9
    ref_p = np.array(out["SIBYLL23D"]["p"]["Nmu"])
    ref_fe = np.array(out["SIBYLL23D"]["Fe"]["Nmu"])

    for m in models:
        ax2.semilogx(e0_ev, np.array(out[m]["p"]["Nmu"]) / ref_p, color=cols[m], ls="-", lw=1.8, label=f"{labels[m]} (p)")
        ax2.semilogx(e0_ev, np.array(out[m]["Fe"]["Nmu"]) / ref_fe, color=cols[m], ls="--", lw=1.8, label=f"{labels[m]} (Fe)")

    ax2.set_xlabel(r"Primary Energy $E_0$ [eV]")
    ax2.set_ylabel(r"$N_\mu / N_\mu^{\mathrm{Sibyll\ 2.3d}}$ ($E_\mu > 500\ \mathrm{GeV}$)")
    ax2.set_title("High-Energy Muon Model Spread")
    ax2.grid(True, which="both", ls=":", alpha=0.5)
    ax2.legend(loc="best", ncol=2, fontsize=8.5)
    ax2.set_ylim(0.78, 1.22)

    plt.tight_layout()
    plt.savefig("plots/hadronic_models_coincidence.png")
    plt.close()

    np.savez("data/mceq_southpole.npz", e0_gev=e0_grid_gev, x_obs=x_obs, results=out)

if __name__ == "__main__":
    solve_southpole_showers()
