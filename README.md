# BBGT-GEM-Maxwell VLBI DSP Pipeline (`v1.0-alpha`)

**A Generalized Non-Invasive Quadrature Residual Filter and Planetary Kernel Demodulation Engine Anchored at GRAO (Kuntunse)**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24323557.svg)](https://doi.org/10.5281/zenodo.24323557)
[![Release](https://img.shields.io/badge/release-v1.0--alpha-blue.svg)](#)
[![Observatory Centroid](https://img.shields.io/badge/centroid-GRAO%20Kuntunse%2032m-green.svg)](#)
[![Networks](https://img.shields.io/badge/networks-AVN%20%7C%20SKA%20%7C%20EVN%20%7C%20VLBA-orange.svg)](#)

## Overview

`bbgt-gem-maxwell-vlbi-pipeline` (`v1.0-alpha`) is a non-invasive digital signal processing (DSP) and astrometric residual-modeling software backend developed at the **Darkon Research Laboratory (Ghana)** for single-dish radio telescopes and Very Long Baseline Interferometry (VLBI) arrays worldwide. Utilizing the near-equatorial **Ghana Radio Astronomy Observatory (GRAO)** 32-meter parabolic antenna at Kuntunse as its primary astrometric centroid and reference anchor, the pipeline operates strictly downstream of standard software correlators (e.g., DiFX, SFXL) on standard AIPS/CASA Random-Groups UVFITS visibility datasets without altering physical telescope feeds, cryogenic low-noise amplifiers (LNAs), or pointing hardware.

### Core DSP Transformations & Physics Modules

1. **Generalized Multi-Observatory ITRF2014 Registry (`OBSERVATORY_REGISTRY`):**
   Computes 3D Euclidean baselines, $(u,v,w)$ coordinates, and geometric delay residuals anchored at `GRAO` (Kuntunse) across partner global stations:
   * `HART` — Hartebeesthoek 26m (South Africa, AVN/EVN): $4,533.169\text{ km}$ ($\tau_{g,\max} = 15.121\text{ ms}$)
   * `SKAM` — SKA-Mid / MeerKAT Core (Karoo, South Africa, SKA/AVN): $4,556.402\text{ km}$ ($\tau_{g,\max} = 15.199\text{ ms}$)
   * `EFLS` — Effelsberg 100m (Germany, EVN): $4,871.411\text{ km}$ ($\tau_{g,\max} = 16.249\text{ ms}$)
   * `PTVL` — VLBA Pie Town 25m (New Mexico, USA, VLBA): $9,842.435\text{ km}$ ($\tau_{g,\max} = 32.831\text{ ms}$)
   * `SKAL` — SKA-Low Core (Murchison, Australia, SKA): $10,849.037\text{ km}$ ($\tau_{g,\max} = 36.188\text{ ms}$)
2. **4:1 Electroweak-to-Gravity Quadrature Bridge ($n=2 \to n=0$):**
   Maps free-space wave impedance ($Z_2 = 376.7303\,\Omega$) to the Gravitational Kernel ground state ($Z_0^{\text{grav}} = 94.1826\,\Omega$) via channelized per-IF $1/4$ amplitude scaling and $+90^\circ$ ($+i$) octonionic quadrature rotation, isolating deterministic geometric delay residuals ($\Delta\tau$) while preserving $100\%$ electromagnetic phase-closure stability ($|V_{\text{EM}}|$).
3. **8:1 Einston Substrate Filter ($n=2 \to n=-1$):**
   Maps free-space reception ($Z_2 = 376.7303\,\Omega$) to the Einston Substrate impedance ($Z_{-1} = 47.0913\,\Omega$), aligning within a $6.18\%$ VSWR offset of standard $50\,\Omega$ coaxial LNA receiver chains.
4. **C-Band / X-Band Planetary Core Density Demodulation ($\omega_k = \sqrt{G_{\text{BBGT}}\rho}$):**
   Uses a Hann-windowed coherent lock-in scan ($G_{\text{BBGT}} = 6.51408576 \times 10^{-11}\text{ m}^3\text{kg}^{-1}\text{s}^{-2}$) to isolate narrowband planetary kernel frequencies—anchored by **Neptune** ($\rho = 1,638\text{ kg/m}^3$, $f_k = 51.988\,\mu\text{Hz}$, $T_k = 5.343\text{ hr}$) at C-Band ($5.0\text{ GHz}$) for **BBGT-$\phi$-Molniya** phase-locking studies, alongside **Jupiter** ($f_k = 46.776\,\mu\text{Hz}$, $T_k = 5.939\text{ hr}$) and **Saturn** ($f_k = 33.669\,\mu\text{Hz}$, $T_k = 8.250\text{ hr}$).

---

## Repository Structure

* `gem_maxwell_pipeline.py` — Core physics constants, global `OBSERVATORY_REGISTRY`, 3D ITRF baseline modeler, 4:1/8:1 quadrature filters, and planetary kernel demodulator.
* `uvfits_benchmark_and_diagnostics.py` — Generalized Random-Groups UVFITS generator and 4-panel verification suite.
* `BBGT_GEM_Maxwell_Verification.ipynb` — Interactive Google Colab / Jupyter computational notebook.
* `PROPOSAL.md` — Technical proposal for archival and C/X-Band verification at GRAO (Kuntunse) and partner AVN/SKA stations.
* `requirements.txt` — Standard Python package requirements (`pip install -r requirements.txt`).
* `Dockerfile` — Reproducible container environment (`python:3.11-slim`).

---

## Quick Start

### Option A: Standard Python (`pip`)
```bash
pip install -r requirements.txt
python3 gem_maxwell_pipeline.py
python3 uvfits_benchmark_and_diagnostics.py
