#!/usr/bin/env python3
"""
gem_maxwell_pipeline.py
GEM-Maxwell Astrometric & Global VLBI Simulation / Post-Processing Engine (v1.0-alpha)
Principal Investigator: Matthias Asare-Darko (Darkon Research Laboratory, Ghana)
Repository: https://github.com/DRL-233/bbgt-gem-maxwell-vlbi-pipeline
Zenodo DOI: 10.5281/zenodo.24323557
Primary Reference Centroid: Ghana Radio Astronomy Observatory (GRAO), Kuntunse (32m Dish / AVN Node)

Implements:
1. Generalized Multi-Observatory ITRF2014 Station Registry anchored at GRAO (Kuntunse),
   supporting AVN (HartRAO), SKA-Mid (Karoo/MeerKAT), SKA-Low (Murchison), EVN (Effelsberg), and VLBA.
2. BBGTGravitationalFilter:
   - 4:1 Electroweak-to-Gravity (n=2 -> n=0, 376.7303 Ohm -> 94.1826 Ohm) channelized quadrature bridge
   - 8:1 Electroweak-to-Einston Substrate (n=2 -> n=-1, 376.7303 Ohm -> 47.0913 Ohm / ~50 Ohm RF) bridge
3. Narrowband Planetary Core Density Kernel Demodulator (omega_k = sqrt(G_BBGT * rho)),
   including NASA standard bulk densities for Solar System planets plus dwarf planets Ceres, Pluto,
   and Eris predicted by the BBGT-Phi-Molniya planetary orbit model.
"""

from dataclasses import dataclass
from typing import Dict, Tuple
import numpy as np

# ==============================================================================
# 1. BBGT & CODATA FUNDAMENTAL CONSTANTS
# ==============================================================================
G_BBGT = 6.51408576e-11       # Darkon / BBGT Gravitational Constant (m^3 kg^-1 s^-2)
C_LIGHT = 299792458.0         # Speed of light in vacuum (m/s)
EPSILON_0 = 8.8541878128e-12  # Vacuum permittivity (F/m)
MU_0 = 1.25663706212e-6       # Vacuum permeability (H/m)
ALPHA_1 = 7.2973525693e-3     # Fine-structure constant (n=1 EM)
R_K = 25812.80745             # von Klitzing constant h/e^2 (Ohms)

# Binary Wave-Impedance Ladder: Z_n = alpha_n * R_K = Z_2 * 2^(n-2)
Z_2_EW = 2.0 * ALPHA_1 * R_K   # 376.73031367 Ohms (n=2 Weak Isospin / Free Space)
Z_1_EM = Z_2_EW / 2.0          # 188.36515683 Ohms (n=1 Pure Electromagnetism)
Z_0_GRAV = Z_2_EW / 4.0        #  94.18257842 Ohms (n=0 Gravitational Kernel)
Z_M1_SUB = Z_2_EW / 8.0        #  47.09128921 Ohms (n=-1 Einston Substrate / ~50 Ohm RF)

# NASA Standard Solar System Bulk Densities (kg/m^3) including Phi-Molniya Dwarf Planets
PLANETARY_DENSITIES: Dict[str, float] = {
    "Sun": 1408.0,
    "Mercury": 5427.0,
    "Venus": 5243.0,
    "Earth": 5514.0,
    "Mars": 3933.0,
    "Ceres": 2162.0,    # NASA Dawn mission standard mean density (2.162 g/cm^3)
    "Jupiter": 1326.0,
    "Saturn": 687.0,
    "Uranus": 1271.0,
    "Neptune": 1638.0,
    "Pluto": 1854.0,    # NASA New Horizons standard mean density (1.854 g/cm^3)
    "Eris": 2430.0,     # NASA / Sicardy et al. occultation standard mean density (2.43 g/cm^3)
}


# ==============================================================================
# 2. GENERALIZED GLOBAL OBSERVATORY REGISTRY (GRAO CENTROID)
# ==============================================================================
@dataclass(frozen=True)
class ITRFStation:
    code: str
    name: str
    network: str
    x_m: float
    y_m: float
    z_m: float

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x_m, self.y_m, self.z_m], dtype=np.float64)


# Primary Equatorial Centroid
GRAO_KUNTUNSE = ITRFStation(
    "GRAO", "GRAO Kuntunse 32m (Ghana)", "AVN / Centroid", 6333076.42, -48231.15, 639188.74
)

# Partner Global VLBI & SKA Reference Stations (ITRF2014 coordinates in meters)
HARTRAO_SA = ITRFStation(
    "HART", "Hartebeesthoek 26m (South Africa)", "AVN / EVN", 5085442.78, 2668263.48, -2768697.03
)
SKA_MID_KAROO = ITRFStation(
    "SKAM", "SKA-Mid / MeerKAT Core (Karoo, SA)", "SKA / AVN", 5109360.13, 2006852.59, -3238948.12
)
SKA_LOW_MURCHISON = ITRFStation(
    "SKAL", "SKA-Low Core (Murchison, Australia)", "SKA", -2559454.08, 5095372.14, -2849057.18
)
EFFELSBERG_EVN = ITRFStation(
    "EFLS", "Effelsberg 100m (Germany)", "EVN", 4033947.24, 486990.89, 4900430.80
)
VLBA_PIETOWN = ITRFStation(
    "PTVL", "VLBA Pie Town 25m (New Mexico, USA)", "VLBA", -1640953.70, -5014816.02, 3575411.88
)

OBSERVATORY_REGISTRY: Dict[str, ITRFStation] = {
    st.code: st
    for st in [
        GRAO_KUNTUNSE,
        HARTRAO_SA,
        SKA_MID_KAROO,
        SKA_LOW_MURCHISON,
        EFFELSBERG_EVN,
        VLBA_PIETOWN,
    ]
}


def get_station(code_or_station: str | ITRFStation) -> ITRFStation:
    """Resolves a station registry code (e.g., 'GRAO', 'SKAM') or returns an ITRFStation directly."""
    if isinstance(code_or_station, ITRFStation):
        return code_or_station
    key = code_or_station.upper().strip()
    if key not in OBSERVATORY_REGISTRY:
        raise KeyError(f"Unknown station '{code_or_station}'. Available: {list(OBSERVATORY_REGISTRY.keys())}")
    return OBSERVATORY_REGISTRY[key]


# ==============================================================================
# 3. 3D ITRF BASELINE & BBGT DELAY MODELER
# ==============================================================================
class VLBISimulationSuite:
    """Models 3D ITRF baselines, geometric delays, and BBGT metric scaling offsets."""

    def __init__(
        self,
        station_a: str | ITRFStation = GRAO_KUNTUNSE,
        station_b: str | ITRFStation = HARTRAO_SA,
    ):
        self.station_a = get_station(station_a)
        self.station_b = get_station(station_b)
        self.baseline_vec = self.station_b.coords - self.station_a.coords
        self.chord_length_m = float(np.linalg.norm(self.baseline_vec))

    def coarse_geometric_delay(self, ha_rad: float, dec_rad: float) -> float:
        """Computes standard geometric light-travel delay tau_g (seconds) for source (HA, Dec)."""
        s_hat = np.array([
            np.cos(dec_rad) * np.cos(ha_rad),
            -np.cos(dec_rad) * np.sin(ha_rad),
            np.sin(dec_rad),
        ])
        w_m = np.dot(self.baseline_vec, s_hat)
        return float(w_m / C_LIGHT)

    def bbgt_residual_delay(self, ha_rad: float, dec_rad: float, delta_b_m: float = 4.2656) -> float:
        """
        Computes BBGT metric scaling delay residual (seconds) along baseline projection.
        Default delta_b_m = +4.2656 m (14.2284 ns max residual on Kuntunse-Hartebeesthoek).
        """
        s_hat = np.array([
            np.cos(dec_rad) * np.cos(ha_rad),
            -np.cos(dec_rad) * np.sin(ha_rad),
            np.sin(dec_rad),
        ])
        proj_factor = np.abs(np.dot(self.baseline_vec, s_hat)) / self.chord_length_m
        return float((delta_b_m * proj_factor) / C_LIGHT)


# ==============================================================================
# 4. BBGT QUADRATURE FILTERS (4:1 GRAVITY & 8:1 EINSTON SUBSTRATE)
# ==============================================================================
class BBGTGravitationalFilter:
    """
    Non-invasive DSP quadrature filter implementing:
    - 4:1 Electroweak-to-Gravity (Z_2 -> Z_0: 376.73 -> 94.18 Ohm) bridge
    - 8:1 Electroweak-to-Einston Substrate (Z_2 -> Z_-1: 376.73 -> 47.09 Ohm) bridge
    """

    def __init__(self, freq_hz: float = 8.4e9):
        self.freq_hz = freq_hz
        self.ratio_4_to_1 = Z_0_GRAV / Z_2_EW   # Exact 0.25 (1/4)
        self.ratio_8_to_1 = Z_M1_SUB / Z_2_EW   # Exact 0.125 (1/8)

    def apply_4to1_quadrature(
        self, visibilities: np.ndarray, delta_tau_s: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Applies 4:1 impedance scaling (1/4) and +90 deg (+1j) octonionic quadrature rotation
        to separate standard EM visibilities from the isolated gravitational channel ||G||.
        """
        phase_shift = 2.0 * np.pi * self.freq_hz * delta_tau_s
        v_em_derotated = visibilities * np.exp(-1j * phase_shift)
        delta_v = visibilities - v_em_derotated
        g_quadrature = 1j * self.ratio_4_to_1 * delta_v
        return v_em_derotated, g_quadrature

    def apply_8to1_substrate_filter(
        self, visibilities: np.ndarray, delta_tau_s: np.ndarray
    ) -> np.ndarray:
        """Applies 8:1 scaling (1/8) matching Z_-1 = 47.0913 Ohm (~50 Ohm coaxial LNA)."""
        phase_shift = 2.0 * np.pi * self.freq_hz * delta_tau_s
        delta_v = visibilities * (1.0 - np.exp(-1j * phase_shift))
        return 1j * self.ratio_8_to_1 * delta_v


# ==============================================================================
# 5. PLANETARY BULK-DENSITY KERNEL DEMODULATOR
# ==============================================================================
class PlanetaryKernelDemodulator:
    """Computes characteristic planetary bulk-density kernel frequencies omega_k = sqrt(G_BBGT * rho)."""

    def __init__(self, g_bbgt: float = G_BBGT):
        self.g_bbgt = g_bbgt

    def kernel_frequency(self, body_name: str) -> Tuple[float, float, float]:
        if body_name not in PLANETARY_DENSITIES:
            raise KeyError(f"Unknown body '{body_name}'. Available: {list(PLANETARY_DENSITIES.keys())}")
        rho = PLANETARY_DENSITIES[body_name]
        omega_k = np.sqrt(self.g_bbgt * rho)
        f_k_hz = omega_k / (2.0 * np.pi)
        period_hr = (1.0 / f_k_hz) / 3600.0
        return float(omega_k), float(f_k_hz), float(period_hr)


# ==============================================================================
# 6. STANDALONE SELF-TEST & GLOBAL CENTROID VERIFICATION
# ==============================================================================
if __name__ == "__main__":
    print("=== 1. BBGT Impedance Bridge Verification ===")
    print(f"Z_2 (Weak Isospin / Free Space): {Z_2_EW:.4f} Ohms")
    print(f"Z_1 (Pure Electromagnetism):     {Z_1_EM:.4f} Ohms")
    print(f"Z_0 (Gravitational Kernel):      {Z_0_GRAV:.4f} Ohms (4:1 ratio = {Z_2_EW/Z_0_GRAV:.1f})")
    print(f"Z_-1 (Einston Substrate / ~50R): {Z_M1_SUB:.4f} Ohms (8:1 ratio = {Z_2_EW/Z_M1_SUB:.1f})\n")

    print("=== 2. Global Baselines Anchored at GRAO (Kuntunse) Centroid ===")
    for code, partner in OBSERVATORY_REGISTRY.items():
        if code == "GRAO":
            continue
        suite = VLBISimulationSuite("GRAO", code)
        max_tau_ms = (suite.chord_length_m / C_LIGHT) * 1e3
        print(f"GRAO <-> {code:4s} ({partner.name:34s}) | Chord = {suite.chord_length_m/1e3:8.3f} km | Max tau_g = {max_tau_ms:6.3f} ms")

    print("\n=== 3. Solar System Planetary Kernel Frequencies (omega_k = sqrt(G_BBGT * rho)) ===")
    demod = PlanetaryKernelDemodulator()
    for body in PLANETARY_DENSITIES:
        w_k, f_k, t_hr = demod.kernel_frequency(body)
        tag = " <-- Phase II Primary Anchor" if body == "Neptune" else (
              " <-- Phi-Molniya Dwarf Planet" if body in ("Ceres", "Pluto", "Eris") else "")
        print(f"{body:10s} | rho = {PLANETARY_DENSITIES[body]:6.1f} kg/m^3 | "
              f"omega_k = {w_k:.4e} rad/s | f_k = {f_k*1e6:6.3f} uHz | T_k = {t_hr:6.3f} hr{tag}")

    print("\n=== 4. Benchmark: Centaurus A (Dec = -43.02 deg, 8.4 GHz, S_0 = 1.25 Jy) on GRAO-HART ===")
    vlbi = VLBISimulationSuite("GRAO", "HART")
    bbgt_filter = BBGTGravitationalFilter(freq_hz=8.4e9)
    dec_cen_a = np.deg2rad(-43.02)
    s0_jy = 1.25

    for ha_hr in [-4.0, 0.0, 4.0]:
        ha_rad = np.deg2rad(ha_hr * 15.0)
        dtau = vlbi.bbgt_residual_delay(ha_rad, dec_cen_a)
        vis_in = np.array([s0_jy * np.exp(1j * 2.0 * np.pi * 8.4e9 * dtau)])
        _, g_out = bbgt_filter.apply_4to1_quadrature(vis_in, np.array([dtau]))
        print(f"HA = {ha_hr:+4.1f}h | dtau = {dtau*1e9:6.3f} ns | ||G|| = {np.abs(g_out[0]):.3f} Jy")
