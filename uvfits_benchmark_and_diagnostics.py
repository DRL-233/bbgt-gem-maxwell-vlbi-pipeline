#!/usr/bin/env python3
"""
uvfits_benchmark_and_diagnostics.py
Phase III Open-Science Blueprint: Generalized Random-Groups UVFITS Generator & Diagnostic Suite (v1.0-alpha)
Principal Investigator: Matthias Asare-Darko (Darkon Research Laboratory, Ghana)
Zenodo DOI: 10.5281/zenodo.24323557
Primary Reference Centroid: Ghana Radio Astronomy Observatory (GRAO), Kuntunse (32m Dish / AVN Node)
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from gem_maxwell_pipeline import (
    C_LIGHT,
    G_BBGT,
    BBGTGravitationalFilter,
    PlanetaryKernelDemodulator,
    VLBISimulationSuite,
)


def compute_uvw_and_delays(
    vlbi: VLBISimulationSuite, ha_rad: np.ndarray, dec_rad: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Computes interferometric (u, v, w) coordinates in meters and BBGT delay residuals (seconds)."""
    bx, by, bz = vlbi.baseline_vec
    sin_h, cos_h = np.sin(ha_rad), np.cos(ha_rad)
    sin_d, cos_d = np.sin(dec_rad), np.cos(dec_rad)

    u_m = bx * sin_h + by * cos_h
    v_m = -bx * sin_d * cos_h + by * sin_d * sin_h + bz * cos_d
    w_m = bx * cos_d * cos_h - by * cos_d * sin_h + bz * sin_d

    dtau_s = np.array([vlbi.bbgt_residual_delay(h, dec_rad) for h in ha_rad])
    return u_m, v_m, w_m, dtau_s


def generate_random_groups_uvfits(
    station_a_code: str = "GRAO",
    station_b_code: str = "HART",
    filename: str | None = None,
    source_name: str = "CEN_A",
    ra_deg: float = 201.36506,
    dec_deg: float = -43.01911,
    freq_hz: float = 8.4e9,
    bandwidth_hz: float = 32.0e6,
    n_channels: int = 16,
    n_times: int = 480,
    flux_jy: float = 1.25,
    noise_sigma_jy: float = 0.015,
    planet_modulation: str | None = None,
    mod_depth_jy: float = 0.045,
) -> dict:
    """
    Writes an AIPS/CASA-compatible Random-Groups UVFITS file for any pair of stations
    in OBSERVATORY_REGISTRY, encoded with BBGT delay residuals and optional planetary modulation.
    """
    rng = np.random.default_rng(42)
    vlbi = VLBISimulationSuite(station_a_code, station_b_code)
    st_a, st_b = vlbi.station_a, vlbi.station_b

    if filename is None:
        freq_ghz = freq_hz / 1e9
        filename = f"{st_a.code}_{st_b.code}_{source_name}_{freq_ghz:.1f}GHz_BBGT.uvfits"

    # 8-hour tracking arc: HA from -4.0h to +4.0h (60-second integrations)
    ha_hours = np.linspace(-4.0, 4.0, n_times)
    time_s = (ha_hours - ha_hours[0]) * 3600.0
    ha_rad = np.deg2rad(ha_hours * 15.0)
    dec_rad = np.deg2rad(dec_deg)

    u_m, v_m, w_m, dtau_s = compute_uvw_and_delays(vlbi, ha_rad, dec_rad)
    chan_freqs = freq_hz + (np.arange(n_channels) - n_channels // 2) * (bandwidth_hz / n_channels)

    # Optional planetary kernel carrier modulation
    planet_signal = np.zeros(n_times, dtype=np.complex128)
    if planet_modulation is not None:
        demod = PlanetaryKernelDemodulator()
        omega_k, _, _ = demod.kernel_frequency(planet_modulation)
        planet_signal = mod_depth_jy * np.exp(1j * omega_k * time_s)

    # Vectorized visibility grid: (n_times, n_channels)
    phases = 2.0 * np.pi * dtau_s[:, np.newaxis] * chan_freqs[np.newaxis, :]
    noise = (
        rng.normal(0.0, noise_sigma_jy, size=(n_times, n_channels))
        + 1j * rng.normal(0.0, noise_sigma_jy, size=(n_times, n_channels))
    )
    v_grid = flux_jy * np.exp(1j * phases) + planet_signal[:, np.newaxis] + noise

    # Random-Groups UVFITS data array: (n_times, 1, 1, n_channels, 2_pols, 3_[Re,Im,Wt])
    vis_data = np.zeros((n_times, 1, 1, n_channels, 2, 3), dtype=np.float32)
    weight = 1.0 / (noise_sigma_jy ** 2)
    for pol in range(2):
        vis_data[:, 0, 0, :, pol, 0] = np.real(v_grid)
        vis_data[:, 0, 0, :, pol, 1] = np.imag(v_grid)
        vis_data[:, 0, 0, :, pol, 2] = weight

    parnames = ["UU", "VV", "WW", "BASELINE", "DATE"]
    pardata = [
        (u_m / C_LIGHT).astype(np.float32),
        (v_m / C_LIGHT).astype(np.float32),
        (w_m / C_LIGHT).astype(np.float32),
        np.full(n_times, 258.0, dtype=np.float32),
        (2461362.0 + time_s / 86400.0).astype(np.float32),
    ]

    group_data = fits.GroupData(vis_data, parnames=parnames, pardata=pardata, bitpix=-32)
    primary_hdu = fits.GroupsHDU(group_data)

    # Standard AIPS/CASA FITS headers
    hdr = primary_hdu.header
    hdr["OBJECT"] = source_name
    hdr["TELESCOP"] = f"{st_a.code}-{st_b.code}"
    hdr["INSTRUME"] = "GEM-MAXWELL-DSP"
    hdr["OBSERVER"] = "M. Asare-Darko (DRL)"
    hdr["OBSRA"] = ra_deg
    hdr["OBSDEC"] = dec_deg
    hdr["CTYPE2"] = "COMPLEX"
    hdr["CRVAL2"] = 1.0
    hdr["CDELT2"] = 1.0
    hdr["CRPIX2"] = 1.0
    hdr["CTYPE3"] = "STOKES"
    hdr["CRVAL3"] = -1.0
    hdr["CDELT3"] = -1.0
    hdr["CRPIX3"] = 1.0
    hdr["CTYPE4"] = "FREQ"
    hdr["CRVAL4"] = freq_hz
    hdr["CDELT4"] = bandwidth_hz / n_channels
    hdr["CRPIX4"] = float(n_channels // 2 + 1)
    hdr["CTYPE5"] = "RA"
    hdr["CRVAL5"] = ra_deg
    hdr["CTYPE6"] = "DEC"
    hdr["CRVAL6"] = dec_deg
    hdr["HISTORY"] = f"Generated by Darkon Research Laboratory BBGT Engine (G_BBGT={G_BBGT:.8e})"

    # Dynamic AIPS AN (Antenna) Binary Table Extension
    an_cols = fits.ColDefs([
        fits.Column(name="ANNAME", format="8A", array=[f"{st_a.code}-ANT", f"{st_b.code}-ANT"]),
        fits.Column(
            name="STABXYZ",
            format="3D",
            unit="METERS",
            array=np.vstack([st_a.coords, st_b.coords]),
        ),
        fits.Column(name="NOSTA", format="1J", array=[1, 2]),
        fits.Column(name="MNTSTA", format="1J", array=[0, 0]),
        fits.Column(name="STAXOF", format="1E", unit="METERS", array=[0.0, 0.0]),
    ])
    an_hdu = fits.BinTableHDU.from_columns(an_cols, name="AIPS AN")
    an_hdu.header["ARRNAM"] = "GLOBAL-VLBI"
    an_hdu.header["FREQ"] = freq_hz

    hdul = fits.HDUList([primary_hdu, an_hdu])
    hdul.writeto(filename, overwrite=True)

    return {
        "filename": filename,
        "baseline_label": f"{st_a.code}–{st_b.code}",
        "ha_hours": ha_hours,
        "time_s": time_s,
        "dtau_s": dtau_s,
        "vis_grid": v_grid,
        "chan_freqs": chan_freqs,
        "freq_hz": freq_hz,
    }


def apply_channelized_bbgt_filter(obs: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Applies 4:1 and 8:1 BBGT filters per IF channel before coherent band averaging."""
    v_grid = obs["vis_grid"]
    dtau_s = obs["dtau_s"]
    chan_freqs = obs["chan_freqs"]

    v_em_ch = np.zeros_like(v_grid)
    g_quad_ch = np.zeros_like(v_grid)
    s_sub_ch = np.zeros_like(v_grid)

    for c_idx, nu in enumerate(chan_freqs):
        filt = BBGTGravitationalFilter(freq_hz=nu)
        v_em_ch[:, c_idx], g_quad_ch[:, c_idx] = filt.apply_4to1_quadrature(v_grid[:, c_idx], dtau_s)
        s_sub_ch[:, c_idx] = filt.apply_8to1_substrate_filter(v_grid[:, c_idx], dtau_s)

    return np.mean(v_em_ch, axis=1), np.mean(g_quad_ch, axis=1), np.mean(s_sub_ch, axis=1)


def demodulate_planetary_spectrum(
    obs: dict, scan_freqs_uhz: np.ndarray, continuum_flux_jy: float
) -> np.ndarray:
    """Removes continuum quadrature envelope and applies a Hann-windowed coherent lock-in scan."""
    _, g_quad, _ = apply_channelized_bbgt_filter(obs)
    phase_center = 2.0 * np.pi * obs["freq_hz"] * obs["dtau_s"]
    continuum_quad_model = 1j * 0.25 * continuum_flux_jy * (1.0 - np.exp(-1j * phase_center))
    residual_quad = g_quad - continuum_quad_model
    residual_quad -= np.mean(residual_quad)

    window = np.hanning(len(obs["time_s"]))
    window_norm = np.sum(window) / 2.0

    return np.array([
        np.abs(np.sum(residual_quad * window * np.exp(-1j * 2.0 * np.pi * (f * 1e-6) * obs["time_s"])))
        / window_norm
        for f in scan_freqs_uhz
    ])


def run_diagnostic_suite(
    station_a: str = "GRAO",
    station_b: str = "HART",
    output_png: str = "GRAO_BBGT_Diagnostic_Report.png",
) -> None:
    """Runs Phase I (Centaurus A) & Phase II (Neptune & Jupiter at C-Band 5.0 GHz) verification."""
    cen_a = generate_random_groups_uvfits(
        station_a_code=station_a,
        station_b_code=station_b,
        source_name="CEN_A",
        dec_deg=-43.02,
        freq_hz=8.4e9,
        flux_jy=1.25,
    )
    v_em, g_quad, s_sub = apply_channelized_bbgt_filter(cen_a)

    neptune_obs = generate_random_groups_uvfits(
        station_a_code=station_a,
        station_b_code=station_b,
        source_name="NEPTUNE",
        dec_deg=-2.5,
        freq_hz=5.0e9,
        flux_jy=1.80,
        planet_modulation="Neptune",
        mod_depth_jy=0.060,
    )
    jupiter_obs = generate_random_groups_uvfits(
        station_a_code=station_a,
        station_b_code=station_b,
        source_name="JUPITER",
        dec_deg=-20.0,
        freq_hz=5.0e9,
        flux_jy=2.50,
        planet_modulation="Jupiter",
        mod_depth_jy=0.060,
    )

    demod = PlanetaryKernelDemodulator()
    _, f_nep_hz, _ = demod.kernel_frequency("Neptune")
    _, f_jup_hz, _ = demod.kernel_frequency("Jupiter")

    scan_freqs_uhz = np.linspace(15.0, 90.0, 1200)
    nep_spec_jy = demodulate_planetary_spectrum(neptune_obs, scan_freqs_uhz, continuum_flux_jy=1.80)
    jup_spec_jy = demodulate_planetary_spectrum(jupiter_obs, scan_freqs_uhz, continuum_flux_jy=2.50)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=150)
    fig.suptitle(
        f"{cen_a['baseline_label']} BBGT-GEM-Maxwell Software Backend Verification Suite (v1.0-alpha)",
        fontsize=13,
        fontweight="bold",
    )

    ax1 = axes[0, 0]
    ax1.plot(cen_a["ha_hours"], cen_a["dtau_s"] * 1e9, color="#1f4e78", lw=2.0)
    ax1.set_title("A. Deterministic BBGT Delay Residual (Cen A, 8.4 GHz)")
    ax1.set_xlabel("Hour Angle HA (hours)")
    ax1.set_ylabel("Delay Residual Δτ (ns)")
    ax1.grid(True, alpha=0.3)

    ax2 = axes[0, 1]
    ax2.plot(cen_a["ha_hours"], np.abs(v_em), color="#2e75b6", lw=1.5, label="Channelized EM |V_EM|")
    ax2.axhline(1.25, color="black", ls="--", lw=1.0, label="Nominal Continuum (1.25 Jy)")
    ax2.set_title("B. Phase-Closure Stability on De-rotated EM Channel")
    ax2.set_xlabel("Hour Angle HA (hours)")
    ax2.set_ylabel("Visibility Amplitude (Jy)")
    ax2.set_ylim(1.20, 1.30)
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)

    ax3 = axes[1, 0]
    ax3.plot(cen_a["ha_hours"], np.abs(g_quad), color="#c55a11", lw=1.6, label="4:1 Gravity ||G|| (94.18 Ω)")
    ax3.plot(cen_a["ha_hours"], np.abs(s_sub), color="#548235", lw=1.6, label="8:1 Einston ||S|| (47.09 Ω)")
    ax3.set_title("C. Isolated Quadrature Channels (4:1 & 8:1 Bridges)")
    ax3.set_xlabel("Hour Angle HA (hours)")
    ax3.set_ylabel("Quadrature Flux Density (Jy)")
    ax3.legend(loc="upper right")
    ax3.grid(True, alpha=0.3)

    ax4 = axes[1, 1]
    ax4.plot(scan_freqs_uhz, nep_spec_jy * 1e3, color="#008080", lw=1.9, label="Neptune C-Band (5.0 GHz)")
    ax4.plot(scan_freqs_uhz, jup_spec_jy * 1e3, color="#7030a0", lw=1.7, label="Jupiter C-Band (5.0 GHz)")
    ax4.axvline(f_nep_hz * 1e6, color="#c00000", ls="--", lw=1.4, label=f"Neptune f_k = {f_nep_hz*1e6:.3f} μHz (5.34 h)")
    ax4.axvline(f_jup_hz * 1e6, color="#7030a0", ls=":", lw=1.3, label=f"Jupiter f_k = {f_jup_hz*1e6:.3f} μHz (5.94 h)")
    ax4.set_title("D. C-Band Planetary Kernel Demodulation (Neptune & Jupiter)")
    ax4.set_xlabel("Modulation Frequency (μHz)")
    ax4.set_ylabel("Demodulated Amplitude (mJy)")
    ax4.legend(loc="upper right", fontsize=8)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    plt.savefig(output_png)
    plt.show()
    print(f"Generated UVFITS benchmarks: {cen_a['filename']}, {neptune_obs['filename']}, {jupiter_obs['filename']}")


if __name__ == "__main__":
    run_diagnostic_suite(station_a="GRAO", station_b="HART")
