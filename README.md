# Surface Hybrid Array Simulation & Air-Shower Reconstruction for IceTop/IceCube

This repository provides a standalone simulation and event reconstruction framework for cosmic-ray extensive air showers (EAS) at the South Pole ($h = 2835\ \mathrm{m}$, atmospheric depth $X \approx 705\ \mathrm{g/cm^2}$). The study addresses two core physical and experimental challenges relevant to the IceTop surface enhancement at the IceCube Neutrino Observatory:

1. **Hadronic Interaction Model Uncertainties:** Solving one-dimensional matrix cascade equations (`MCEq`) in the South Pole atmosphere to quantify high-energy in-ice muon yields ($E_\mu > 500\ \mathrm{GeV}$) and surface electromagnetic shower sizes ($N_e$) across post-LHC hadronic interaction models (`Sibyll 2.3d`, `EPOS-LHC`, and `QGSJet-II-04`).
2. **Maximum-Likelihood Event Reconstruction:** Modeling the spatial and temporal response of a surface scintillation detector array and reconstructing the shower core position $(x_c, y_c)$, arrival direction $(\theta, \phi)$, and lateral distribution scaling parameter $S_{125}$ using `iminuit` (`MINUIT2`).

---

## Repository Structure

```text
├── src/
│   ├── hadronic_mceq.py        # South Pole atmospheric cascade solver (MCEq + MSIS-00)
│   ├── icetop_array.py         # Surface array geometry, NKG LDF, timing front, and radio proxy
│   └── run_reconstruction.py   # Two-step Chi2 / Log-Likelihood reconstruction via iminuit
├── plots/                      # Diagnostic and reconstruction benchmark figures
├── data/                       # Precomputed MCEq lookup tables (mceq_southpole.npz)
└── requirements.txt            # Python package dependencies
```

---

## Physics & Methodology

### 1. Atmospheric Cascade Equations & Surface / In-Ice Coincidence

The longitudinal evolution of the hadronic and leptonic cascade in the South Pole atmosphere (modeled via the December `MSIS00_IC` profile) is governed by the coupled one-dimensional matrix cascade equation solved by `MCEq`:

$$\frac{\mathrm{d}\mathbf{\Phi}}{\mathrm{d}X} = -\left( \mathbf{\Lambda}_{\mathrm{int}} + \frac{1}{\rho(X)} \mathbf{\Lambda}_{\mathrm{dec}} \right) \mathbf{\Phi} + \mathbf{S}_{\mathrm{prod}} \mathbf{\Phi}$$

Air showers initiated by primary protons ($A=1$) and iron nuclei ($A=56$) are simulated across $3 \times 10^{14}\ \mathrm{eV} \le E_0 \le 10^{17}\ \mathrm{eV}$. The high-energy muon bundle yield capable of penetrating the glacial ice down to the deep IceCube optical array ($E_\mu > 500\ \mathrm{GeV}$) and the corresponding Gaisser-Hillas surface electromagnetic cascade size $N_e$ evaluated at $X_{\mathrm{obs}} = 705.3\ \mathrm{g/cm^2}$ are given by:

$$N_\mu(E_\mu > 500\ \mathrm{GeV}) = \int_{500\ \mathrm{GeV}}^{E_0} \left[ \Phi_{\mu^-}(E, X_{\mathrm{obs}}) + \Phi_{\mu^+}(E, X_{\mathrm{obs}}) \right] \mathrm{d}E$$

$$N_e(X_{\mathrm{obs}}) = N_{\max} \left( \frac{X_{\mathrm{obs}} - X_0}{X_{\max} - X_0} \right)^{\frac{X_{\max} - X_0}{\lambda}} \exp\left( \frac{X_{\max} - X_{\mathrm{obs}}}{\lambda} \right)$$

![Hadronic Models Coincidence](plots/hadronic_models_coincidence.png)

* **Left Panel (Mass Separation):** Following the Heitler-Matthews superposition scaling ($N_\mu \propto A^{1-\beta} E_0^\beta$), heavier primaries ($A=56$) yield systematically higher TeV muon multiplicities for a fixed surface electromagnetic size $N_e$ compared to proton primaries ($A=1$).
* **Right Panel (Model Spread):** Highlights the **10% to 15%** systematic spread in high-energy muon predictions among `Sibyll 2.3d`, `EPOS-LHC`, and `QGSJet-II-04`.

---

### 2. Surface Array Response, LDF & Shower-Front Timing

Surface scintillation detector stations are arranged on a regular $7 \times 7$ grid with $50\ \mathrm{m}$ spacing. For an incoming shower axis $\hat{\mathbf{n}} = (-\sin\theta\cos\phi, -\sin\theta\sin\phi, -\cos\theta)$ and core location $\mathbf{r}_c = (x_c, y_c)$, station coordinates $\mathbf{r}_i$ are projected into the perpendicular shower plane:

$$d_{\parallel, i} = (\mathbf{r}_i - \mathbf{r}_c) \cdot \hat{\mathbf{n}}, \qquad r_{\perp, i} = \sqrt{\Vert{}\mathbf{r}_i - \mathbf{r}_c\Vert{}^2 - d_{\parallel, i}^2}$$

Local particle densities are parameterized using a Nishimura-Kamata-Greisen (NKG) lateral distribution function normalized to the standard IceTop reference distance $R_{\mathrm{ref}} = 125\ \mathrm{m}$ ($S_{125}$), where Molière radius $R_{\mathrm{M}} = 125\ \mathrm{m}$ and shower age $s = 1.25$:

$$S(r_\perp) = S_{125} \, \frac{f_{\mathrm{NKG}}(r_\perp)}{f_{\mathrm{NKG}}(125\ \mathrm{m})}, \qquad f_{\mathrm{NKG}}(r_\perp) = \left(\frac{r_\perp}{R_{\mathrm{M}}}\right)^{s - 2} \left(1 + \frac{r_\perp}{R_{\mathrm{M}}}\right)^{s - 4.5}$$

Station signals incorporate Poisson photoelectron statistics, Gaussian electronic baseline noise, and a hardware trigger threshold of $0.3\ \mathrm{MIP}$. Expected arrival times follow a planar propagation front modified by a parabolic shower-front curvature delay ($a_{\mathrm{curv}} = 1.5 \times 10^{-4}\ \mathrm{ns/m^2}$) and signal-dependent timing resolution $\sigma_{t, i} = 2.5\ \mathrm{ns} / \sqrt{S_i}$:

$$t_{\mathrm{exp}, i} = t_0 + \frac{d_{\parallel, i}}{c} + a_{\mathrm{curv}} \, r_{\perp, i}^2$$

A geomagnetic radio amplitude proxy scaling with the Lorentz force ($E_{\mathrm{radio}} \propto S_{125} \sin\alpha \, \exp(-r_\perp / r_0)$, where $\alpha$ is the geomagnetic angle for the South Pole field inclination $I \approx -72^\circ$) is also implemented in the array module.

![Simulated Event Footprint](plots/array_event_footprint.png)

---

### 3. Maximum-Likelihood Event Reconstruction (`iminuit`)

For triggered events ($\ge 5$ stations above $0.3\ \mathrm{MIP}$), parameter estimation is carried out in two stages using `iminuit` (`MINUIT2`), initialized from a signal-weighted center-of-gravity seed $\mathbf{r}_{\mathrm{CoG}} = \sum S_i^{1/2} \mathbf{r}_i / \sum S_i^{1/2}$:

1. **Shower Axis Direction $(\theta, \phi)$:** Minimizes the timing $\chi^2$ functional over triggered stations:

$$\chi^2_{\mathrm{time}}(t_0, n_x, n_y) = \sum_{i \in \mathrm{hits}} \left( \frac{t_{\mathrm{meas}, i} - t_{\mathrm{exp}, i}(t_0, n_x, n_y)}{\sigma_{t, i}} \right)^2$$

2. **Shower Core $(x_c, y_c)$ & Scaling Parameter $S_{125}$:** Minimizes a combined negative log-likelihood ($\mathrm{NLL}$) accounting for both the charge variance of triggered stations and the Poisson zero-hit probability of silent (sub-threshold) stations:

$$\mathrm{NLL}(x_c, y_c, \log_{10} S_{125}) = \frac{1}{2} \sum_{i \in \mathrm{hits}} \left[ \frac{\left(S_{\mathrm{meas}, i} - S(r_{\perp, i})\right)^2}{\sigma_{S, i}^2} + \ln\left(\sigma_{S, i}^2\right) \right] + \sum_{j \in \mathrm{silent}} \mu_{\mathrm{eff}}(r_{\perp, j})$$

![Reconstruction Performance](plots/reconstruction_performance.png)

#### Monte Carlo Benchmark Summary ($N = 250$ Showers, $5^\circ \le \theta \le 35^\circ$)

| Observable / Metric | Symbol | 68% Containment / Performance |
| :--- | :---: | :---: |
| **Shower Core Position Resolution** | $\Delta R_{68}$ | $0.52\ \mathrm{m}$ |
| **Angular Resolution (Point Spread Function)** | $\Psi_{68}$ | $0.20^\circ$ |
| **Energy Proxy Recovery Ratio** | $S_{125}^{\mathrm{rec}} / S_{125}^{\mathrm{true}}$ | $1.00 \pm 0.01$ (Unbiased) |

---

## Installation & Usage

Clone the repository and install the required scientific Python packages:

```bash
git clone [https://github.com/parsa-arabzadeh/icetop-hybrid-reco.git](https://github.com/parsa-arabzadeh/icetop-hybrid-reco.git)
cd icetop-hybrid-reco
pip install -r requirements.txt
```

Run the simulation and reconstruction modules sequentially from the repository root:

```bash
python3 src/hadronic_mceq.py
python3 src/icetop_array.py
python3 src/run_reconstruction.py
```

---

## References
* A. Fedynitch, R. Engel, T. K. Gaisser, F. Riehn, and T. Stanev, *Calculation of conventional and prompt lepton fluxes at very high energy*, EPJ Web Conf. **99**, 08001 (2015) [`MCEq`].
* F. Riehn, R. Engel, A. Fedynitch, T. K. Gaisser, and T. Stanev, *Hadronic interaction model Sibyll 2.3d and limits on prompt lepton fluxes*, Phys. Rev. D **102**, 063002 (2020).
* R. Abbasi *et al.* (IceCube Collaboration), *IceTop: The surface component of IceCube*, Nucl. Instrum. Meth. A **700**, 188–220 (2013).

---

## Author
**Parsa Arabzadeh**  
MSc in Physics
University of Tehran
