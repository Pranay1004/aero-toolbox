"""
reentry_lib.py
==============
Reentry Aerothermodynamics Calculator Library
Implements every numbered equation from "Re-entry Missions" Ch.15,
Sections 15.4.4 (Reentry Heating) and 15.5 (Reentry Mission Management).

Equation map (docstring on every function tells you exactly which
equation(s) it implements, by number, so you can trace every result
straight back to the book):

  15.28  q_dot = dQ/dt
  15.29  q_dot = f * D * V = f * rho0 * sigma * CD * S * V^3 / 2
  15.30  q_dot ~ sigma * V^3
  15.31  sigma = rho_inf / rho0
  15.32  Q = integral(q_dot dt) = q_dot_ave * dt
  15.33  Q ~ Ve^2 * sqrt( (m/(CD S)) * (1/sin(gamma_e)) )        [ballistic]
  15.34  Q ~ Ve^2 * sqrt( BC * (1/sin(gamma_e)) )                [ballistic]
  15.35  q_dot_max ~ V^3 * sqrt( BC * sin(gamma_e) )             [ballistic]
  15.36  Q ~ Ve^2 * sqrt( (L/D) * BC * (1/sin(gamma_e)) )        [lifting]
  15.37  q_dot_max ~ V^3 * sqrt( BC * sin(gamma_e) / (L/D) )     [lifting]
  15.38  q_t = 18300*sqrt(rho_inf)/sqrt(RN) * (V_inf/1e4)^3.05   [Scott et al., stagnation, W/cm^2]
  15.39  q_t = 11030/sqrt(RN) * (rho_inf/rho0)^0.5 * (V_inf/V0)^3.15  [alt. stagnation model]
  15.40  q_r = 100*RN*(V_inf/1e4)^8.5*(rho_inf/rho0)^1.6         [Martin, radiative stagnation]
  15.41  Q = m * Cv * dT                                        [Heat-Sink TPS]
  15.42  Q = hv * dm                                             [Ablative TPS]
  15.43  q_out = q_in ~= eps * sigma_SB * Tw^4                    [Radiative-equilibrium TPS]
  15.44  s = -integral( V dV / (D/m) )                            [range, velocity form]
  15.45  E = g*h + V^2/2                                          [specific energy]
  15.46  s = -integral( dE / D )                                  [range, energy form]
  15.47  (L/m) cos(sigma) <= [g - V^2/r] cos(gamma)                [equilibrium-glide / load constraint]
  15.48  q_t = (11030/sqrt(RN)) (rho/rho0)^0.5 (V/Vcir)^3.15 <= q_t_max   [heat-flux path constraint]
  15.49  n_z = (D/m) [1+(L/D)^2]^0.5 <= n_max                     [normal load-factor constraint]
  15.50  q = 0.5 * rho * V^2 <= q_max                              [dynamic pressure constraint]

Author: built for Pranay's reentry / TPS simulation workflow.
"""

import numpy as np

# ----------------------------------------------------------------------
# Input validation helpers
# ----------------------------------------------------------------------
def _validate_positive(val, name):
    """Raise ValueError if val <= 0."""
    if np.any(np.asarray(val) <= 0):
        raise ValueError(f"{name} must be positive, got {val}")

def _validate_non_negative(val, name):
    """Raise ValueError if val < 0."""
    if np.any(np.asarray(val) < 0):
        raise ValueError(f"{name} must be non-negative, got {val}")

def _validate_range(val, lo, hi, name):
    """Raise ValueError if val not in [lo, hi]."""
    arr = np.asarray(val)
    if np.any(arr < lo) or np.any(arr > hi):
        raise ValueError(f"{name} must be in [{lo}, {hi}], got {val}")

# ----------------------------------------------------------------------
# Physical constants & Multi-Planetary Configuration Database
# ----------------------------------------------------------------------
RHO0_SEALEVEL = 1.225          # kg/m^3, Earth sea-level density (Eq. 15.31, 15.39 rho0)
V0_CIRC = 7905.0               # m/s, Earth circular orbital velocity at sea level (Eq. 15.39/15.48)
G0 = 9.80665                   # m/s^2, Earth sea-level gravity
R_EARTH = 6.378e6              # m, Earth radius
SIGMA_SB = 5.670374e-8         # W/m^2-K^4, Stefan-Boltzmann constant (Eq. 15.43)
H_SCALE = 7200.0               # m, atmospheric scale height
R_GAS = 287.05287              # J/(kg·K), specific gas constant for dry air
GAMMA_AIR = 1.4                # ratio of specific heats for air
M_AIR = 0.0289644              # kg/mol, molar mass of dry air
R_UNIVERSAL = 8.314462618      # J/(mol·K), universal gas constant

# Multi-Planetary Database
PLANET_CONFIGS = {
    'Earth': {
        'name': 'Earth',
        'R_planet': 6.378e6,      # m
        'g0': 9.80665,            # m/s^2
        'rho0': 1.225,            # kg/m^3
        'H_scale': 7200.0,        # m
        'V_circ': 7905.0,         # m/s
        'R_gas': 287.05287,       # J/(kg K)
        'gamma': 1.40,
        'M_gas': 0.0289644,       # kg/mol
        'k_conv': 1.7415e-4,      # W s^0.5 / (m^1.5 kg^0.5)
        'h_entry': 120.0e3,       # m
        'P0': 101325.0,           # Pa
        'T0': 288.15,             # K
        'composition': {'N2': 0.7808, 'O2': 0.2095, 'Ar': 0.0093, 'CO2': 0.0004}
    },
    'Mars': {
        'name': 'Mars',
        'R_planet': 3.3895e6,
        'g0': 3.72076,
        'rho0': 0.020,
        'H_scale': 11100.0,
        'V_circ': 3550.0,
        'R_gas': 191.80,
        'gamma': 1.289,
        'M_gas': 0.04334,
        'k_conv': 1.9027e-4,
        'h_entry': 125.0e3,
        'P0': 610.0,
        'T0': 210.0,
        'composition': {'CO2': 0.9532, 'N2': 0.0270, 'Ar': 0.0160}
    },
    'Venus': {
        'name': 'Venus',
        'R_planet': 6.0518e6,
        'g0': 8.870,
        'rho0': 65.0,
        'H_scale': 15900.0,
        'V_circ': 7327.0,
        'R_gas': 191.30,
        'gamma': 1.286,
        'M_gas': 0.04345,
        'k_conv': 1.8900e-4,
        'h_entry': 200.0e3,
        'P0': 9.2e6,
        'T0': 737.0,
        'composition': {'CO2': 0.9650, 'N2': 0.0350}
    }
}


def get_planet_config(planet='Earth'):
    """Return physical parameters dict for given target body ('Earth', 'Mars', 'Venus')."""
    return PLANET_CONFIGS.get(str(planet).strip().capitalize(), PLANET_CONFIGS['Earth'])


# ----------------------------------------------------------------------
# International Standard Atmosphere & Multi-Planet Atmosphere Models
# ----------------------------------------------------------------------
_ISA_LAYERS = [
    (0,       288.15,  -0.0065),    # Troposphere
    (11000,   216.65,   0.0),       # Tropopause (isothermal)
    (20000,   216.65,  0.001),      # Lower stratosphere
    (32000,   228.65,  0.0028),     # Upper stratosphere
    (47000,   270.65,  0.0),        # Stratopause (isothermal)
    (51000,   270.65, -0.0028),     # Upper mesosphere
    (71000,   214.65, -0.002),      # Mesopause
    (84852,   186.87,  0.0),        # Lower thermosphere (isothermal)
]


def isa_temperature(h_m, planet='Earth'):
    """
    Atmospheric temperature [K] at geometric altitude h_m [m] for target planet.
    """
    h_m = float(max(h_m, 0.0))
    cfg = get_planet_config(planet)
    
    if cfg['name'] == 'Earth':
        for i in range(len(_ISA_LAYERS) - 1):
            h_base, T_base, L = _ISA_LAYERS[i]
            h_next = _ISA_LAYERS[i + 1][0]
            if h_m < h_next:
                return T_base if abs(L) < 1e-10 else T_base + L * (h_m - h_base)
        h_base, T_base, L = _ISA_LAYERS[-1]
        return T_base if abs(L) < 1e-10 else T_base + L * (h_m - h_base)
    elif cfg['name'] == 'Mars':
        # NASA Glenn Mars Atmosphere Model approximation
        T = 210.0 - 0.00099 * min(h_m, 70000.0)
        return float(max(T, 140.0))
    else: # Venus
        # VIRA Model approximation
        T = max(737.0 - 0.008 * min(h_m, 65000.0), 180.0)
        return float(T)


def isa_pressure(h_m, planet='Earth'):
    """
    Atmospheric pressure [Pa] at geometric altitude h_m [m] for target planet.
    """
    h_m = float(max(h_m, 0.0))
    cfg = get_planet_config(planet)
    
    if cfg['name'] == 'Earth':
        P = 101325.0
        for i in range(len(_ISA_LAYERS) - 1):
            h_base, T_base, L = _ISA_LAYERS[i]
            h_next = _ISA_LAYERS[i + 1][0]
            if h_m <= h_next:
                if abs(L) < 1e-10:
                    return P * np.exp(-G0 * M_AIR * (h_m - h_base) / (R_UNIVERSAL * T_base))
                else:
                    exponent = -G0 * M_AIR / (R_UNIVERSAL * L)
                    return P * ((T_base + L * (h_m - h_base)) / T_base) ** exponent
            else:
                if abs(L) < 1e-10:
                    P = P * np.exp(-G0 * M_AIR * (h_next - h_base) / (R_UNIVERSAL * T_base))
                else:
                    exponent = -G0 * M_AIR / (R_UNIVERSAL * L)
                    P = P * ((T_base + L * (h_next - h_base)) / T_base) ** exponent
        h_base, T_base, L = _ISA_LAYERS[-1]
        if abs(L) < 1e-10:
            return P * np.exp(-G0 * M_AIR * (h_m - h_base) / (R_UNIVERSAL * T_base))
        else:
            exponent = -G0 * M_AIR / (R_UNIVERSAL * L)
            return P * ((T_base + L * (h_m - h_base)) / T_base) ** exponent
    else:
        # Hydrostatic hydrostatic profile for Mars / Venus using exponential density scale
        rho = isa_density(h_m, planet=planet)
        return float(rho * cfg['R_gas'] * isa_temperature(h_m, planet=planet))


def isa_density(h_m, planet='Earth'):
    """
    Atmospheric density [kg/m^3] at geometric altitude h_m [m] for target planet.
    """
    cfg = get_planet_config(planet)
    h_m = float(max(h_m, 0.0))
    if cfg['name'] == 'Earth':
        T = isa_temperature(h_m, planet='Earth')
        P = isa_pressure(h_m, planet='Earth')
        return float(P / (R_GAS * T))
    else:
        return float(cfg['rho0'] * np.exp(-h_m / cfg['H_scale']))


def isa_speed_of_sound(h_m, planet='Earth'):
    """
    Speed of sound [m/s] at altitude h_m [m] for target planet: a = sqrt(gamma * R * T).
    """
    cfg = get_planet_config(planet)
    T = isa_temperature(h_m, planet=planet)
    return float(np.sqrt(cfg['gamma'] * cfg['R_gas'] * T))


def isa_dynamic_viscosity(h_m, planet='Earth'):
    """
    Dynamic viscosity [Pa·s] using Sutherland's law for specified planetary gas.
    """
    T = isa_temperature(h_m, planet=planet)
    T_ref = 273.15
    mu_ref = 1.716e-5
    S = 110.4
    return mu_ref * (T / T_ref) ** 1.5 * (T_ref + S) / (T + S)


# ----------------------------------------------------------------------
# Wrapper: density from altitude (uses ISA)
# ----------------------------------------------------------------------
def _density_at_km(h_km):
    """Density [kg/m^3] at one altitude in km. ISA to 86 km, exponential above."""
    h_km = float(h_km)
    if h_km <= 0:
        return RHO0_SEALEVEL
    if h_km <= 86.0:
        return max(isa_density(h_km * 1000.0), 1e-14)
    # Extend above ISA with exponential (H ≈ 7.5 km near mesopause)
    rho_86 = max(isa_density(86000.0), 1e-14)
    rho = rho_86 * np.exp(-(h_km - 86.0) / 7.5)
    return max(rho, 1e-30)  # floor prevents underflow to 0.0 at extreme altitudes


def density_exponential(h_km):
    """
    Free-stream density rho_inf [kg/m^3] at altitude h_km [km].
    ISA 1976 below ~86 km; exponential extension above for entry altitudes.
    Name kept for backward compatibility.
    """
    h_km = np.atleast_1d(np.asarray(h_km, dtype=float))
    rho = np.array([_density_at_km(h) for h in h_km])
    return rho if rho.size > 1 else float(rho[0])


# Cached altitude–density tables for fast inversion (built once)
_H_TABLE_KM = None
_RHO_TABLE = None


def _ensure_atm_tables():
    global _H_TABLE_KM, _RHO_TABLE
    if _H_TABLE_KM is None:
        _H_TABLE_KM = np.linspace(0.0, 150.0, 1501)
        _RHO_TABLE = np.array([_density_at_km(h) for h in _H_TABLE_KM])


def altitude_from_density(rho_target, h_grid_km=None):
    """
    Inverse of density function: altitude [km] for target density [kg/m^3].
    Returns altitude in kilometres (never meters). Vectorized over rho_target.
    """
    _ensure_atm_tables()
    if h_grid_km is None:
        h_grid = _H_TABLE_KM
        rho_grid = _RHO_TABLE
    else:
        h_grid = np.asarray(h_grid_km, dtype=float)
        rho_grid = density_exponential(h_grid)
    # rho decreases with h → reverse for np.interp
    out = np.interp(np.asarray(rho_target, dtype=float), rho_grid[::-1], h_grid[::-1])
    return float(out) if np.ndim(out) == 0 or out.size == 1 else out


# ----------------------------------------------------------------------
# Isentropic Relations
# ----------------------------------------------------------------------
def isentropic_T_ratio(M, gamma=GAMMA_AIR):
    """T_0/T = 1 + (gamma-1)/2 * M^2  (stagnation to static temperature ratio)."""
    return 1.0 + (gamma - 1.0) / 2.0 * np.asarray(M) ** 2


def isentropic_P_ratio(M, gamma=GAMMA_AIR):
    """P_0/P = (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1))  (stagnation to static pressure)."""
    return isentropic_T_ratio(M, gamma) ** (gamma / (gamma - 1.0))


def isentropic_rho_ratio(M, gamma=GAMMA_AIR):
    """rho_0/rho = (1 + (gamma-1)/2 * M^2)^(1/(gamma-1))  (stagnation to static density)."""
    return isentropic_T_ratio(M, gamma) ** (1.0 / (gamma - 1.0))


def isentropic_mach(P0_over_P, gamma=GAMMA_AIR):
    """Mach number from stagnation-to-static pressure ratio."""
    return np.sqrt(2.0 / (gamma - 1.0) * (np.asarray(P0_over_P) ** ((gamma - 1.0) / gamma) - 1.0))


def isentropic_T0(T_static, M, gamma=GAMMA_AIR):
    """Stagnation temperature T0 = T * (1 + (gamma-1)/2 * M^2)."""
    return np.asarray(T_static) * isentropic_T_ratio(M, gamma)


def isentropic_P0(P_static, M, gamma=GAMMA_AIR):
    """Stagnation pressure P0 = P * (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1))."""
    return np.asarray(P_static) * isentropic_P_ratio(M, gamma)


def isentropic_rho0(rho_static, M, gamma=GAMMA_AIR):
    """Stagnation density rho0 = rho * (1 + (gamma-1)/2 * M^2)^(1/(gamma-1))."""
    return np.asarray(rho_static) * isentropic_rho_ratio(M, gamma)


# ----------------------------------------------------------------------
# Gravity model
# ----------------------------------------------------------------------
def gravity(h_m, planet='Earth'):
    """g(h) using inverse-square law (used inside Eq. 15.45, 15.47)."""
    cfg = get_planet_config(planet)
    return cfg['g0'] * (cfg['R_planet'] / (cfg['R_planet'] + h_m)) ** 2


# ----------------------------------------------------------------------
# Eq. 15.28 - 15.32 : basic heat-rate / heat-load definitions
# ----------------------------------------------------------------------
def density_ratio(rho_inf, rho0=RHO0_SEALEVEL):
    """Eq. 15.31: sigma = rho_inf / rho0."""
    return np.asarray(rho_inf) / rho0


def heat_rate_total(f, rho0, sigma, CD, S, V):
    """
    Eq. 15.29/15.30: q_dot = f * rho0 * sigma * CD * S * V^3 / 2
    Total (not per-area) instantaneous heat-transfer rate to the
    vehicle [W], f = energy-conversion factor (Fig. 15.26, user input),
    rho0 sea-level density, sigma density ratio, CD drag coeff,
    S reference area [m^2], V velocity [m/s].
    """
    _validate_positive(rho0, "rho0")
    _validate_non_negative(CD, "CD")
    _validate_positive(S, "S")
    _validate_non_negative(V, "V")
    return f * rho0 * sigma * CD * S * np.asarray(V) ** 3 / 2.0


def total_heat_load_from_rate(q_dot_ave, dt):
    """Eq. 15.32: Q = q_dot_ave * delta_t  (or integral over time)."""
    return q_dot_ave * dt


def _trapz(y, x):
    """NumPy 1.x/2.x compatible trapezoidal integral."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def total_heat_load_integral(t, q_dot):
    """
    Eq. 15.28/15.32 exact form: Q = integral(q_dot dt).
    t [s], q_dot same units as desired Q/time (e.g. W/cm² → Q in J/cm²).
    """
    return _trapz(q_dot, t)


# ----------------------------------------------------------------------
# Eq. 15.33 - 15.37 scaling laws  →  dimensional via Allen–Eggers + Eq.15.39
# ----------------------------------------------------------------------
# The book prints 15.33–15.37 as proportionalities (no constant). Absolute
# values below integrate the dimensional stagnation model (Eq. 15.39) along
# an Allen–Eggers ρ–V path that preserves the same scaling in BC, γe, L/D.
# ----------------------------------------------------------------------
def allen_eggers_path(Ve, BC, gamma_e_deg, L_D=0.0, V_end=400.0, H=H_SCALE, n=600):
    """
    Simplified constant-γ reentry path (Allen–Eggers).
    ρ = −(2 BC sinγ / H) ln(V/Ve), with path stretched for L/D > 0.
    Returns V [m/s], rho [kg/m³], dt-compatible arrays (descending V).
    """
    gamma = abs(float(np.radians(gamma_e_deg)))
    sγ = max(np.sin(gamma), 1e-4)
    # Lifting: more time in atmosphere (matches 15.36 ∝ √(L/D) trend)
    BC_eff = BC / max(1.0 + 0.5 * max(float(L_D), 0.0), 1e-6)
    V = np.linspace(float(Ve), float(V_end), int(n))
    # avoid log(0)
    ratio = np.clip(V / float(Ve), 1e-12, 1.0)
    rho = -(2.0 * BC_eff * sγ / H) * np.log(ratio)
    rho = np.maximum(rho, 1e-14)
    return V, rho


def stagnation_heat_load_path(RN, Ve, BC, gamma_e_deg, L_D=0.0, V_end=400.0, H=H_SCALE, n=600):
    """
    Total stagnation-point heat load [J/cm²] = ∫ q_t dt along Allen–Eggers path,
    with q_t from Eq. 15.39 (W/cm²). Uses dV/dt ≈ −½ ρ V² / BC.
    Captures book trends 15.33–15.36 with physical units.
    """
    V, rho = allen_eggers_path(Ve, BC, gamma_e_deg, L_D=L_D, V_end=V_end, H=H, n=n)
    q = stagnation_heat_alt(rho, V, RN)  # W/cm²
    # dt = −dV / (½ ρ V² / BC)
    dV = np.diff(V)
    V_m = 0.5 * (V[:-1] + V[1:])
    rho_m = 0.5 * (rho[:-1] + rho[1:])
    q_m = 0.5 * (q[:-1] + q[1:])
    accel = 0.5 * rho_m * V_m ** 2 / max(BC, 1e-9)
    dt = np.where(accel > 1e-12, -dV / accel, 0.0)
    dt = np.maximum(dt, 0.0)
    # Lifting multiplies exposure (√(L/D) family of 15.36)
    LD = max(float(L_D), 0.0)
    stretch = np.sqrt(1.0 + LD) if LD > 0 else 1.0
    return float(np.sum(q_m * dt) * stretch)  # J/cm²


def stagnation_peak_rate_path(RN, Ve, BC, gamma_e_deg, L_D=0.0, V_end=400.0, H=H_SCALE, n=600):
    """Peak stagnation heat rate [W/cm²] on Allen–Eggers path (Eq. 15.39)."""
    V, rho = allen_eggers_path(Ve, BC, gamma_e_deg, L_D=L_D, V_end=V_end, H=H, n=n)
    q = stagnation_heat_alt(rho, V, RN)
    return float(np.max(q))


def ballistic_heat_load(Ve, BC, gamma_e_deg, RN=1.0, **_kw):
    """
    Dimensional ballistic total heat load [J/cm²].
    Scaling of Eq. 15.34 with absolute values from Eq. 15.39 path integral.
    Optional RN (m); default 1 m blunt body.
    """
    return stagnation_heat_load_path(RN, Ve, BC, gamma_e_deg, L_D=0.0)


def ballistic_heat_load_massform(Ve, m, CD, S, gamma_e_deg, RN=1.0, **_kw):
    """Eq. 15.33 form with BC = m/(CD S); returns J/cm²."""
    return ballistic_heat_load(Ve, m / (CD * S), gamma_e_deg, RN=RN)


def ballistic_max_heat_rate(V, BC, gamma_e_deg, RN=1.0, **_kw):
    """
    Peak heat rate [W/cm²] for ballistic entry (Eq. 15.35 trend).
    V is entry velocity Ve [m/s].
    """
    return stagnation_peak_rate_path(RN, V, BC, gamma_e_deg, L_D=0.0)


def lifting_heat_load(Ve, L_D, BC, gamma_e_deg, RN=1.0, **_kw):
    """Dimensional lifting total heat load [J/cm²] (Eq. 15.36 trend)."""
    if np.any(np.asarray(L_D) <= 0):
        raise ValueError("L_D must be positive for lifting heat load")
    return stagnation_heat_load_path(RN, Ve, BC, gamma_e_deg, L_D=L_D)


def lifting_max_heat_rate(V, BC, gamma_e_deg, L_D, RN=1.0, **_kw):
    """Peak heat rate [W/cm²] for lifting entry (Eq. 15.37 trend)."""
    if np.any(np.asarray(L_D) <= 0):
        raise ValueError("L_D must be positive for lifting max heat rate")
    return stagnation_peak_rate_path(RN, V, BC, gamma_e_deg, L_D=L_D)


# ----------------------------------------------------------------------
# Eq. 15.38 - 15.40 : stagnation-point convective & radiative heating
# ----------------------------------------------------------------------
def stagnation_heat_fay_riddell(rho_inf, V_inf, RN):
    """
    Eq. 15.38 (Scott et al. / Fay-Riddell correlation): q_t = 18300 * sqrt(rho_inf) / sqrt(RN) *
    (V_inf/1e4)^3.05
    rho_inf [kg/m^3], V_inf [m/s], RN [m]  ->  q_t in W/cm^2.
    Note (from text): this model consistently reads HIGH vs. measured data.
    """
    _validate_positive(rho_inf, "rho_inf")
    _validate_non_negative(V_inf, "V_inf")
    _validate_positive(RN, "RN")
    return 18300.0 * np.sqrt(rho_inf) / np.sqrt(RN) * (np.asarray(V_inf) / 1.0e4) ** 3.05


def stagnation_heat_allen_eggers(rho_inf, V_inf, RN, rho0=RHO0_SEALEVEL, V0=V0_CIRC):
    """
    Eq. 15.39: q_t = 11030/sqrt(RN) * (rho_inf/rho0)^0.5 * (V_inf/V0)^3.15
    rho_inf [kg/m^3], V_inf [m/s], RN [m], rho0 sea-level density
    [kg/m^3], V0 circular orbital velocity at sea level (=7950 m/s).
    Returns q_t in W/cm^2.
    """
    _validate_positive(rho_inf, "rho_inf")
    _validate_non_negative(V_inf, "V_inf")
    _validate_positive(RN, "RN")
    return (11030.0 / np.sqrt(RN)) * (np.asarray(rho_inf) / rho0) ** 0.5 * \
           (np.asarray(V_inf) / V0) ** 3.15


# Backward-compatibility aliases
stagnation_heat_scott = stagnation_heat_fay_riddell
stagnation_heat_alt = stagnation_heat_allen_eggers


def sutherland_viscosity(T):
    """
    Sutherland's law for the dynamic viscosity of air [Pa*s].
    T in Kelvin.
    """
    mu_0 = 1.716e-5
    T_0 = 273.15
    S = 110.4
    T = np.maximum(T, 1e-6)
    return mu_0 * ((T / T_0) ** 1.5) * ((T_0 + S) / (T + S))


def air_specific_heat(T):
    """
    Specific heat capacity of air Cp [J/(kg*K)].
    Uses an average specific heat for high temperatures (air up to 2500K).
    """
    return 1100.0


def stagnation_heat_fay_riddell_full(rho_inf, V_inf, T_inf, p_inf, RN, eps, max_iter=20):
    """
    High-fidelity Full Fay-Riddell convective heating.
    Returns qt (W/cm^2), Tw (K), hw (J/kg), He (J/kg)
    It iteratively solves for Tw and qw.
    """
    Cp = air_specific_heat(T_inf)
    h_inf = Cp * T_inf
    He = h_inf + 0.5 * V_inf**2
    
    gamma = 1.4
    M_inf_sq = (V_inf**2) / (gamma * R_GAS * max(T_inf, 1e-6))
    if M_inf_sq < 1.0:
        return 0.0, 300.0, Cp*300.0, He
    
    p_e = p_inf * (1.0 + 2.0 * gamma / (gamma + 1.0) * (M_inf_sq - 1.0))
    rho_e = rho_inf * ((gamma + 1.0) * M_inf_sq) / (2.0 + (gamma - 1.0) * M_inf_sq)
    T_e = p_e / (rho_e * R_GAS)
    mu_e = sutherland_viscosity(T_e)
    
    du_e_dx = (1.0 / RN) * np.sqrt(2.0 * max(p_e - p_inf, 0.0) / max(rho_e, 1e-6))
    Pr_term = 0.763 * (0.71 ** -0.6)
    
    Tw = 300.0
    qw_W_m2 = 0.0
    eps_safe = max(eps, 1e-6)
    
    for _ in range(max_iter):
        hw = air_specific_heat(Tw) * Tw
        rho_w = p_e / (R_GAS * Tw)
        mu_w = sutherland_viscosity(Tw)
        
        dh = max(He - hw, 0.0)
        qw_W_m2 = Pr_term * ((rho_e * mu_e)**0.4) * ((rho_w * mu_w)**0.1) * dh * np.sqrt(du_e_dx)
        
        Tw_new = (qw_W_m2 / (eps_safe * SIGMA_SB)) ** 0.25
        Tw_new = max(Tw_new, T_inf)
        
        if abs(Tw_new - Tw) < 1.0:
            Tw = Tw_new
            break
        Tw = 0.5 * Tw + 0.5 * Tw_new
        
    return qw_W_m2 / 10000.0, Tw, hw, He


def radiative_heat_martin(RN, V_inf, rho_inf, rho0=RHO0_SEALEVEL):
    """
    Eq. 15.40 (Martin): q_r = 100 * RN * (V_inf/1e4)^8.5 * (rho_inf/rho0)^1.6
    Gas-to-surface radiative heat rate at very high entry velocity.
    Returns q_r in W/cm^2 (RN in m).
    """
    _validate_positive(RN, "RN")
    _validate_non_negative(V_inf, "V_inf")
    _validate_positive(rho_inf, "rho_inf")
    return 100.0 * RN * (np.asarray(V_inf) / 1.0e4) ** 8.5 * (np.asarray(rho_inf) / rho0) ** 1.6


# ----------------------------------------------------------------------
# Lecture #2: Stagnation Point and Distributed Heating & Analytical Entry Equations
# ----------------------------------------------------------------------
def stagnation_heat_sutton_graves(rho_inf, V_inf, RN, k=1.7415e-4):
    """
    Sutton-Graves convective stagnation heating model (Lecture #2 page 45/54):
    q_s = k * sqrt(rho_inf / RN) * V_inf^3
    With k in SI units, returns W/m^2. We convert to W/cm^2 by multiplying by 1e-4.
    k defaults to 1.7415e-4 (standard Earth value). For Mars, use k = 1.9027e-4.
    """
    _validate_non_negative(rho_inf, "rho_inf")
    _validate_non_negative(V_inf, "V_inf")
    _validate_positive(RN, "RN")
    return 1e-4 * k * np.sqrt(np.asarray(rho_inf) / RN) * np.asarray(V_inf) ** 3


def local_heating_sphere(q_stag, theta_deg):
    """
    Distributed Heating along a spherical body (Lecture #2 page 33):
    q / q_stag = cos(theta)
    Returns local heat flux q [same units as q_stag] for angle theta_deg (degrees).
    """
    theta_rad = np.radians(np.asarray(theta_deg))
    return np.asarray(q_stag) * np.cos(theta_rad)


def total_heating_rate_sphere(q_stag, RN):
    """
    Total integrated heat rate into a hemisphere (Lecture #2 page 33):
    Q_dot = q_stag * (projected area of hemisphere) = q_stag * pi * RN^2
    Assumes q_stag is in W/cm^2 and RN is in meters.
    Returns total heat rate in Watts [W].
    """
    _validate_positive(RN, "RN")
    # projected area in cm^2: pi * (RN * 100)^2 = 10000 * pi * RN^2
    area_cm2 = np.pi * (np.asarray(RN) * 100.0) ** 2
    return np.asarray(q_stag) * area_cm2


def ballistic_peak_heating_density_analytical(BC, gamma_e_deg, H=H_SCALE):
    """
    Analytical density at peak convective heating for ballistic entry (Lecture #2 page 43):
    rho* = BC * sin(gamma_e) / (3 * H)
    BC: ballistic coefficient [kg/m^2], gamma_e_deg: entry angle [deg]
    """
    gamma_e = abs(np.radians(np.asarray(gamma_e_deg, dtype=float)))
    s_gamma = np.maximum(np.sin(gamma_e), 1e-4)
    return np.asarray(BC, dtype=float) * s_gamma / (3.0 * H)


def ballistic_peak_heating_altitude_analytical(BC, gamma_e_deg, H=H_SCALE, rho0=RHO0_SEALEVEL):
    """
    Analytical altitude at peak convective heating for ballistic entry (Lecture #2 page 44):
    h* = -H * ln( BC * sin(gamma_e) / (3 * H * rho0) )
    BC: ballistic coefficient [kg/m^2], gamma_e_deg: entry angle [deg]
    """
    rho_star = ballistic_peak_heating_density_analytical(BC, gamma_e_deg, H=H)
    ratio = rho_star / rho0
    # Handle array or scalar
    if np.ndim(ratio) == 0:
        return -H * np.log(ratio) / 1000.0 if ratio > 0 else np.nan
    else:
        out = np.full(ratio.shape, np.nan)
        ok = ratio > 0
        out[ok] = -H * np.log(ratio[ok]) / 1000.0
        return out


def ballistic_peak_heating_velocity_analytical(Vatm):
    """
    Analytical velocity at peak convective heating for ballistic entry (Lecture #2 page 44):
    V* = Vatm * exp(-1/6) = 0.8465 * Vatm
    """
    return np.asarray(Vatm, dtype=float) * np.exp(-1.0 / 6.0)


def ballistic_max_heat_rate_analytical(Vatm, BC, gamma_e_deg, RN, H=H_SCALE, k=1.7415e-4):
    """
    Analytical peak stagnation convective heat rate for ballistic entry (Lecture #2 page 45/48):
    q_s,max = k * sqrt(1/RN) * sqrt(BC * sin(gamma_e) / (3 * H)) * (0.846 * Vatm)^3 * 1e-4
    Returns W/cm^2.
    """
    rho_star = ballistic_peak_heating_density_analytical(BC, gamma_e_deg, H=H)
    V_star = ballistic_peak_heating_velocity_analytical(Vatm)
    return stagnation_heat_sutton_graves(rho_star, V_star, RN, k=k)


def ballistic_heat_load_analytical(Vatm, BC, gamma_e_deg, RN, H=H_SCALE, k=1.7415e-4):
    """
    Analytical total stagnation heat load for ballistic entry (Lecture #2 page 48):
    Q_s = k * Vatm^2 * sqrt(pi * H * BC / (RN * sin(gamma_e))) * 1e-4
    Returns J/cm^2.
    """
    gamma_e = abs(np.radians(np.asarray(gamma_e_deg, dtype=float)))
    s_gamma = np.maximum(np.sin(gamma_e), 1e-4)
    Q_m2 = k * (np.asarray(Vatm, dtype=float) ** 2) * np.sqrt(np.pi * H * np.asarray(BC, dtype=float) / (np.asarray(RN, dtype=float) * s_gamma))
    return Q_m2 * 1e-4


def lifting_max_heat_rate_analytical(Ve, BC, L_D, RN):
    """
    Analytical peak convective heat rate for lifting entry / equilibrium glide (Lecture #2 page 53):
    q_s,max = 1.94 * sqrt( BC / (RN * |L/D|) )  [W/cm^2] (after converting 1.94e4 W/m^2 -> 1.94 W/cm^2)
    """
    _validate_positive(RN, "RN")
    ld_abs = np.abs(np.asarray(L_D, dtype=float))
    _validate_positive(ld_abs, "abs(L_D)")
    return 1.94 * np.sqrt(np.asarray(BC, dtype=float) / (np.asarray(RN, dtype=float) * ld_abs))


def lifting_heat_load_analytical(Vatm, BC, L_D, RN, Vc=V0_CIRC):
    """
    Analytical total heat load for lifting entry / equilibrium glide (Lecture #2 page 53):
    Q_s = 2050 * sqrt( BC / (RN * |L/D|) ) * [ arcsin(Vatm/Vc) - (Vatm/Vc)*sqrt(1 - (Vatm/Vc)^2) ]  [J/cm^2]
    (after converting 2.05e7 J/m^2 -> 2050 J/cm^2)
    """
    _validate_positive(RN, "RN")
    ld_abs = np.abs(np.asarray(L_D, dtype=float))
    _validate_positive(ld_abs, "abs(L_D)")
    ratio = np.asarray(Vatm, dtype=float) / Vc
    ratio = np.clip(ratio, 0.0, 1.0)
    bracket = np.arcsin(ratio) - ratio * np.sqrt(1.0 - ratio ** 2)
    return 2050.0 * np.sqrt(np.asarray(BC, dtype=float) / (np.asarray(RN, dtype=float) * ld_abs)) * bracket


# ----------------------------------------------------------------------
# Eq. 15.41 - 15.43 : TPS sizing (heat-sink / ablative / radiative-eq)
# ----------------------------------------------------------------------
def heat_sink_mass(Q, Cv, dT):
    """Eq. 15.41: Q = m*Cv*dT  ->  m = Q / (Cv * dT).  [Heat-Sink TPS]"""
    return Q / (Cv * dT)


def ablative_mass_loss(Q, hv):
    """Eq. 15.42: Q = hv*dm  ->  dm = Q / hv.  [Ablative TPS]"""
    return Q / hv


def radiative_equilibrium_temp(q_in, eps, sigma_sb=SIGMA_SB):
    """
    Eq. 15.43: q_out = q_in ~= eps*sigma*Tw^4  ->  Tw = (q_in/(eps*sigma))^0.25
    q_in in W/m^2 (convert W/cm^2 * 1e4 if using Eq. 15.38-15.40 outputs).
    Returns Tw in Kelvin.
    """
    return (np.asarray(q_in) / (eps * sigma_sb)) ** 0.25


# ----------------------------------------------------------------------
# Eq. 15.44 - 15.46 : range / specific-energy relations
# ----------------------------------------------------------------------
def specific_energy(h, V, planet='Earth'):
    """Eq. 15.45: E = g*h + V^2/2  (uses altitude-varying gravity)."""
    h = np.asarray(h, dtype=float)
    V = np.asarray(V, dtype=float)
    return gravity(h, planet=planet) * h + V ** 2 / 2.0


def range_from_velocity(V, D_over_m):
    """
    Eq. 15.44: s = −∫ V dV / (D/m)  [m]
    V [m/s] descending; D_over_m [m/s²] drag acceleration.
    """
    V = np.asarray(V, dtype=float)
    D_over_m = np.asarray(D_over_m, dtype=float)
    return -_trapz(V / np.maximum(D_over_m, 1e-12), V)


def range_from_energy(E, D_over_m):
    """
    Eq. 15.46: s = −∫ dE / (D/m)  [m]
    E [J/kg] specific energy; D here is drag acceleration [m/s²].
    """
    E = np.asarray(E, dtype=float)
    D_over_m = np.asarray(D_over_m, dtype=float)
    return -_trapz(1.0 / np.maximum(D_over_m, 1e-12), E)


# ----------------------------------------------------------------------
# Eq. 15.47 - 15.50 : reentry-corridor path constraints
# ----------------------------------------------------------------------
def equilibrium_glide_density(V, L_D, BC, h_m=None, planet='Earth'):
    """
    Eq. 15.47 with bank σ=0 (full lift-up), cos γ ≈ 1:
        L/m = g − V²/r
        L/m = (L/D) · (½ ρ V² / BC)
        → ρ_eq = 2 BC (g − V²/r) / ((L/D) V²)
    Upper corridor boundary (min density / max altitude).
    If h_m is None, iterate once so g and r match the altitude of ρ_eq.
    """
    V = np.asarray(V, dtype=float)
    LD = max(float(L_D), 1e-6)
    cfg = get_planet_config(planet)
    R_p = cfg['R_planet']

    def _rho_at(h_guess_m):
        g = gravity(h_guess_m, planet=planet)
        rr = R_p + h_guess_m
        num = 2.0 * BC * (g - V ** 2 / rr)
        den = LD * V ** 2
        return np.where(num > 0, num / den, np.nan)

    if h_m is not None:
        return _rho_at(float(h_m))

    # Fixed-point: guess h, update from ρ, recompute once (vectorized)
    rho = _rho_at(40000.0)
    rho_arr = np.atleast_1d(np.asarray(rho, dtype=float))
    h_est = np.full(rho_arr.shape, 40000.0)
    ok = np.isfinite(rho_arr) & (rho_arr > 0)
    if np.any(ok):
        h_est[ok] = altitude_from_density(rho_arr[ok]) * 1000.0
    g = gravity(h_est, planet=planet)
    rr = R_p + h_est
    V1 = np.atleast_1d(V)
    num = 2.0 * BC * (g - V1 ** 2 / rr)
    den = LD * V1 ** 2
    rho2 = np.where(num > 0, num / den, np.nan)
    return float(rho2[0]) if np.ndim(V) == 0 else rho2


def heatflux_boundary_density(V, RN, q_t_max, rho0=RHO0_SEALEVEL, Vcir=V0_CIRC):
    """
    Eq. 15.48 inverted at q_t = q_t_max [W/cm²]:
        ρ = ρ0 ( q_t_max √R_N / 11030 / (V/V_cir)^3.15 )²
    Max density (min altitude) allowed by heat flux.
    """
    _validate_positive(V, "V")
    _validate_positive(RN, "RN")
    _validate_positive(q_t_max, "q_t_max")
    V = np.asarray(V, dtype=float)
    ratio = (q_t_max * np.sqrt(RN) / 11030.0) / (V / Vcir) ** 3.15
    return rho0 * ratio ** 2


def loadfactor_boundary_density(V, BC, L_D, n_max, g=G0):
    """
    Eq. 15.49 inverted. n_max is in g's (as in textbook: n_max = 2.5 g).

        n_z = (D/m) √(1+(L/D)²) / g  ≤  n_max
        →  ρ = 2 BC (n_max · g) / (V² √(1+(L/D)²))

    Returns max density allowed by normal load.
    """
    _validate_positive(V, "V")
    _validate_positive(BC, "BC")
    _validate_positive(n_max, "n_max")
    V = np.asarray(V, dtype=float)
    return 2.0 * BC * (n_max * g) / (V ** 2 * np.sqrt(1.0 + float(L_D) ** 2))


def dynpress_boundary_density(V, q_max):
    """Eq. 15.50 inverted: ρ = 2 q_max / V²  (q_max in Pa)."""
    _validate_positive(V, "V")
    _validate_positive(q_max, "q_max")
    V = np.asarray(V, dtype=float)
    return 2.0 * q_max / V ** 2


def reentry_corridor(V, RN, BC, L_D, q_t_max, n_max, q_max,
                      rho0=RHO0_SEALEVEL, Vcir=V0_CIRC, planet='Earth'):
    """
    Full reentry corridor (Fig. 15.34 / 15.35):
      upper altitude = equilibrium glide (Eq. 15.47, σ=0)
      lower altitude = most restrictive of heat / load / dynamic pressure
    Returns densities and altitudes in **km** (not meters).
    n_max in g's; q_t_max in W/cm²; q_max in Pa.
    """
    V = np.asarray(V, dtype=float)
    rho_upper = equilibrium_glide_density(V, L_D, BC, planet=planet)
    rho_heat = heatflux_boundary_density(V, RN, q_t_max, rho0=rho0, Vcir=Vcir)
    rho_load = loadfactor_boundary_density(V, BC, L_D, n_max)
    rho_q = dynpress_boundary_density(V, q_max)
    # most restrictive lower bound = smallest allowed max-density
    rho_lower = np.minimum.reduce([rho_heat, rho_load, rho_q])

    def _h(rho_arr):
        rho_arr = np.asarray(rho_arr, dtype=float)
        out = np.full(rho_arr.shape, np.nan, dtype=float)
        ok = np.isfinite(rho_arr) & (rho_arr > 0)
        if np.any(ok):
            out[ok] = altitude_from_density(rho_arr[ok])
        return out

    h_upper = _h(rho_upper)
    h_lower = _h(rho_lower)
    h_heat = _h(rho_heat)
    h_load = _h(rho_load)
    h_q = _h(rho_q)

    corridor_width = np.where(
        np.isfinite(h_upper) & np.isfinite(h_lower),
        np.maximum(h_upper - h_lower, 0.0),
        np.nan,
    )

    return dict(
        V=V,
        rho_upper=rho_upper, rho_lower=rho_lower,
        rho_heat=rho_heat, rho_load=rho_load, rho_q=rho_q,
        h_upper=h_upper, h_lower=h_lower,
        h_heat=h_heat, h_load=h_load, h_q=h_q,
        corridor_width=corridor_width,
    )


def ballistic_coefficient(m, CD, S):
    """BC = m / (Cd · S)  [kg/m²]."""
    return m / (CD * S)


def normal_load_factor(D_over_m, L_D, g=G0):
    """Eq. 15.49 as load factor in g's: n_z = (D/m)√(1+(L/D)²) / g."""
    return (np.asarray(D_over_m) * np.sqrt(1.0 + np.asarray(L_D) ** 2)) / g


def drag_from_density(rho, V, BC):
    """Drag acceleration D/m [m/s²] = ½ ρ V² / BC."""
    return 0.5 * np.asarray(rho) * np.asarray(V) ** 2 / BC


def dynamic_pressure(rho, V):
    """Eq. 15.50 (value): q = 0.5 * rho * V^2."""
    return 0.5 * np.asarray(rho) * np.asarray(V) ** 2


# ----------------------------------------------------------------------
# Representative vehicle-class presets, used to drive the stagnation-
# point-heating and heat-load comparisons across BALLISTIC, LIFTING and
# WINGED reentry classes (Sec. 15.4.4, Eq. 15.33-15.40 applied per class).
# Values are illustrative (typical order-of-magnitude), not a specific
# named vehicle's certified data -- edit freely for your own study.
# ----------------------------------------------------------------------
VEHICLE_PRESETS = {
    "ballistic": dict(
        label="Ballistic capsule",
        RN=1.0, BC=400.0, L_D=0.0, gamma_e=20.0, Ve=11000.0,
    ),
    "lifting": dict(
        label="Lifting body (low L/D)",
        RN=1.0, BC=150.0, L_D=0.3, gamma_e=1.5, Ve=7900.0,
    ),
    "winged": dict(
        label="Winged body (high L/D)",
        RN=0.3, BC=100.0, L_D=1.2, gamma_e=1.0, Ve=7800.0,
    ),
}


# ----------------------------------------------------------------------
# Qualitative-figure generators (Fig. 15.26, 15.29, 15.32, 15.36).
# The book gives these as empirical/illustrative charts with NO closed-
# form governing equation (f, the energy-conversion factor, is only
# ever used as an input to Eq. 15.29; the alpha- and bank-angle
# profiles and the ranging strategy are guidance-law *strategies*, not
# equations). These helpers synthesize representative curves of the
# correct qualitative shape purely so the figures can be regenerated
# and swept parametrically -- they are explicitly NOT derived from a
# governing equation in the text.
# ----------------------------------------------------------------------
def energy_conversion_factor_profile(h_km):
    """
    Qualitative reproduction of Fig. 15.26 f(altitude):
    free-molecular heating (f -> ~1) at high altitude, decreasing
    through laminar convection to a minimum (f ~ 0.01) around
    20-40 km, then a turbulent-convection uptick, then continuing
    down at low altitude (hot-gas radiation regime).
    Returns f (dimensionless, 0-1) for use in Eq. 15.29.
    """
    h = np.atleast_1d(np.asarray(h_km, dtype=float))
    f_free_molecule = 1.0 / (1.0 + np.exp(-(h - 95.0) / 8.0))
    f_laminar_min = 0.012 + 0.02 * np.exp(-((h - 30.0) / 12.0) ** 2)
    turbulent_bump = 0.05 * np.exp(-((h - 15.0) / 6.0) ** 2)
    f = np.maximum(f_free_molecule, f_laminar_min + turbulent_bump)
    f = np.clip(f, 1e-3, 1.0)
    return f if f.size > 1 else float(f[0])


def alpha_profile(mach, alpha_hi=40.0, alpha_lo=5.0, mach_break_hi=20.0, mach_break_lo=8.0):
    """
    Qualitative reproduction of Fig. 15.29: constant high angle of
    attack at high Mach, linear ramp-down through the break region,
    constant low angle of attack at low Mach. Returns alpha [deg].
    """
    M = np.atleast_1d(np.asarray(mach, dtype=float))
    alpha = np.where(
        M >= mach_break_hi, alpha_hi,
        np.where(M <= mach_break_lo, alpha_lo,
                 alpha_lo + (alpha_hi - alpha_lo) * (M - mach_break_lo) / (mach_break_hi - mach_break_lo))
    )
    return alpha if alpha.size > 1 else float(alpha[0])


def bank_angle_sawtooth(mach, period_mach=6.0, sigma_max=70.0):
    """
    Qualitative reproduction of Fig. 15.32: a sawtooth bank-angle
    profile that reverses sign periodically in Mach to null cross
    range while retaining down-range control authority.
    Returns sigma [deg].
    """
    M = np.atleast_1d(np.asarray(mach, dtype=float))
    saw = sigma_max * (2.0 * (M / period_mach - np.floor(M / period_mach + 0.5)))
    sign_flip = np.where((np.floor(M / period_mach).astype(int) % 2) == 0, 1.0, -1.0)
    sigma = saw * sign_flip
    return sigma if sigma.size > 1 else float(sigma[0])


def ranging_drag_profile(V, shift_frac=0.0, scale_frac=0.0):
    """
    Qualitative reproduction of Fig. 15.36 ranging strategy: shifts and
    scales a nominal trapezoidal drag-acceleration-vs-velocity profile
    to represent short-range (compressed/raised) or long-range
    (stretched/lower) adjustments, per the text's description of
    onboard magnitude adjustment of the stored reference profile.
    """
    Vn = np.linspace(7500, 0, 400)
    nominal = 3.0 + 17.0 / (1.0 + np.exp(-(7000 - Vn) / 700)) * (1.0 / (1.0 + np.exp((500 - Vn) / 700)))
    nominal = np.clip(nominal, 0, None)
    shifted_V = Vn * (1.0 + shift_frac)
    scaled = nominal * (1.0 + scale_frac)
    V = np.asarray(V, dtype=float)
    return np.interp(V, shifted_V[::-1], scaled[::-1])


def corridor_drag_segments(V, D_max=None, V1=7500.0, V2=5500.0, V3=4500.0,
                           V4=1500.0, D_min_quad=2.0, BC=100.0, L_D=1.2,
                           q_t_max=60.0, RN=0.3):
    """
    Generate the segmented drag-acceleration profile shown in the textbook
    Fig. 15.35 style corridor plot:

      Segment 1: 1st Quadratic (V1..V2) — rising from V1 to peak
      Segment 2: 2nd Quadratic (V2..V3) — flat-top constant drag plateau
      Segment 3: Constant Drag (V3..V4) — linear descent
      Segment 4: Linear Energy (V4..0) — tail-off

    If D_max is None, it is auto-computed as the minimum of:
      - load-factor limit at V_mid (2.5 g typical)
      - heat-flux limit at V_mid (Eq. 15.48)
    using the provided BC, L_D, q_t_max, RN.

    Parameters
    ----------
    V : array-like  — velocity array [m/s] (descending from V1 to 0)
    D_max : float or None — peak/constant drag acceleration [m/s²]
    V1 : float      — entry interface velocity [m/s]
    V2 : float      — start of constant-drag plateau [m/s]
    V3 : float      — end of constant-drag plateau [m/s]
    V4 : float      — start of linear-energy segment [m/s]
    D_min_quad : float — drag value at V=V1 (start of 1st quadratic)
    BC : float      — ballistic coefficient [kg/m²] (used if D_max is None)
    L_D : float     — lift-to-drag ratio (used if D_max is None)
    q_t_max : float — max stagnation heat flux [W/cm²] (used if D_max is None)
    RN : float      — nose radius [m] (used if D_max is None)

    Returns
    -------
    D : ndarray — drag acceleration at each V [m/s²]
    """
    V = np.asarray(V, dtype=float)

    if D_max is None:
        V_mid = (V1 + V3) / 2.0
        n_max = 2.5  # g, typical structural limit
        rho_load = loadfactor_boundary_density(V_mid, BC, L_D, n_max)
        D_from_load = 0.5 * rho_load * V_mid ** 2 / BC
        rho_heat = heatflux_boundary_density(V_mid, RN, q_t_max)
        D_from_heat = 0.5 * rho_heat * V_mid ** 2 / BC
        D_max = min(float(D_from_load), float(D_from_heat))
        D_max = max(D_max, D_min_quad + 1.0)

    D = np.zeros_like(V)

    # Segment 1: 1st Quadratic (V1 -> V2): parabolic rise from D_min_quad to D_max
    mask1 = (V <= V1) & (V > V2)
    t1 = (V1 - V[mask1]) / (V1 - V2)  # t1 goes 0..1 as V goes V1..V2
    D[mask1] = D_min_quad + (D_max - D_min_quad) * t1 ** 2

    # Segment 2: 2nd Quadratic / Constant Drag plateau (V2 -> V3)
    mask2 = (V <= V2) & (V > V3)
    D[mask2] = D_max

    # Segment 3: Linear descent (V3 -> V4)
    mask3 = (V <= V3) & (V > V4)
    t3 = (V3 - V[mask3]) / (V3 - V4)  # t3 goes 0..1 as V goes V3..V4
    D[mask3] = D_max * (1.0 - t3)

    # Segment 4: Linear Energy tail (V4 -> 0)
    mask4 = (V <= V4) & (V >= 0)
    t4 = V[mask4] / V4
    D[mask4] = D_max * 0.15 * t4  # small residual drag

    return D


# ----------------------------------------------------------------------
# 3-DOF Trajectory Integration (point-mass, spherical target body)
# ----------------------------------------------------------------------
def trajectory_integrate(V0, gamma0_deg, h0_km, RN, BC, L_D,
                         bank_angle_deg=0.0,
                         t_max=7200.0, dt=0.5, planet='Earth'):
    """
    3-DOF point-mass trajectory for a reentry vehicle.

    Equations of motion (non-rotating spherical planet):
        dh/dt     = V * sin(gamma)
        dV/dt     = -D/m - g*sin(gamma)
        dgamma/dt = (L*cos(sigma))/(m*V) - (g/V - V/r)*cos(gamma)

    Parameters
    ----------
    V0 : float         Entry velocity [m/s]
    gamma0_deg : float Entry flight-path angle [deg] (negative = descending)
    h0_km : float      Entry altitude [km]
    RN : float         Nose radius [m]
    BC : float         Ballistic coefficient m/(CD*S) [kg/m^2]
    L_D : float        Lift-to-drag ratio
    bank_angle_deg : float Bank angle [deg], 0 = full lift-up, 90 = no vertical lift
    t_max : float      Max simulation time [s] (0 or negative for auto-estimation)
    dt : float         Time step [s]
    planet : str       Target planet ('Earth', 'Mars', 'Venus')
    """
    _validate_positive(V0, "V0")
    _validate_positive(RN, "RN")
    _validate_positive(BC, "BC")
    if dt <= 0:
        raise ValueError(f"dt must be positive, got {dt}")

    cfg = get_planet_config(planet)
    R_p = cfg['R_planet']
    g0 = cfg['g0']
    k_conv = cfg['k_conv']
    R_gas = cfg['R_gas']
    gamma_gas = cfg['gamma']

    sigma = np.radians(bank_angle_deg)
    gamma0 = np.radians(gamma0_deg)
    H_FLOOR = 0.0
    V_MIN = 1.0

    # Auto-calculate expected reentry duration if t_max is 0 or negative
    if t_max is None or t_max <= 0:
        gamma_mag = abs(gamma0)
        sγ = max(np.sin(gamma_mag), 1e-4)
        t_est = (h0_km * 1000.0) / max(V0 * sγ, 1.0)
        if L_D > 0.1:
            t_max = max(min(t_est * 8.0, 4800.0), 1200.0)
        else:
            t_max = max(min(t_est * 3.0, 1500.0), 600.0)

    def derivs(h, V, gamma):
        h = max(h, H_FLOOR + 10.0)
        r = R_p + h
        g = g0 * (R_p / r) ** 2
        rho = isa_density(h, planet=planet)
        q = 0.5 * rho * V ** 2
        D_over_m = q / BC
        L_over_m = L_D * D_over_m

        dhdt = V * np.sin(gamma)
        dVdt = -D_over_m - g * np.sin(gamma)
        lift_vert = L_over_m * np.cos(sigma)
        grav_term = (g / V - V / r) * np.cos(gamma)
        dgdt = (lift_vert / V) - grav_term

        return [dhdt, dVdt, dgdt]

    def _rk4_step(h, V, gamma, dt_):
        k1 = derivs(h, V, gamma)
        k2 = derivs(h + 0.5*dt_*k1[0], V + 0.5*dt_*k1[1], gamma + 0.5*dt_*k1[2])
        k3 = derivs(h + 0.5*dt_*k2[0], V + 0.5*dt_*k2[1], gamma + 0.5*dt_*k2[2])
        k4 = derivs(h + dt_*k3[0], V + dt_*k3[1], gamma + dt_*k3[2])
        h_new = h + dt_/6.0*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        V_new = V + dt_/6.0*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        g_new = gamma + dt_/6.0*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
        return h_new, V_new, g_new

    # --- Accumulate trajectory ---
    t_all, h_all, V_all, gamma_all = [], [], [], []
    h_m = h0_km * 1000.0
    V = V0
    gamma = gamma0
    t_elapsed = 0.0
    max_steps = int(t_max / dt)
    stop_reason = "t_max reached"

    for step_i in range(max_steps):
        h_m, V, gamma = _rk4_step(h_m, V, gamma, dt)
        V = max(V, 0.0)
        t_elapsed += dt

        t_all.append(t_elapsed)
        h_all.append(h_m)
        V_all.append(V)
        gamma_all.append(gamma)

        if h_m <= H_FLOOR:
            stop_reason = "ground impact"
            break
        if V < V_MIN:
            stop_reason = "velocity depleted"
            break

        if h_m > 1e6 and t_elapsed > 200.0:
            stop_reason = "escaped orbit"
            break

        glide_window_s = 200.0
        glide_window_steps = max(int(glide_window_s / dt), 100)
        if t_elapsed > glide_window_s + 50.0 and len(h_all) > glide_window_steps:
            h_recent = np.array(h_all[-glide_window_steps:])
            dh_recent = np.abs(h_recent - h_recent[0])
            if np.max(dh_recent) < 50.0 and h_m > 30000.0:
                stop_reason = "sustained high-altitude equilibrium glide"
                break

    t_arr = np.array(t_all)
    h_arr = np.array(h_all)
    V_arr = np.array(V_all)
    gamma_arr = np.array(gamma_all)

    # Derived quantities
    h_arr = np.clip(h_arr, 0, None)
    V_arr = np.clip(V_arr, 0, None)

    rho_arr = np.array([isa_density(h, planet=planet) for h in h_arr])
    q_dyn = 0.5 * rho_arr * V_arr ** 2
    D_over_m_arr = q_dyn / BC
    g_arr = np.array([g0 * (R_p / (R_p + h)) ** 2 for h in h_arr])
    T_arr = np.array([isa_temperature(h, planet=planet) for h in h_arr])
    p_arr = np.array([isa_pressure(h, planet=planet) for h in h_arr])
    a_arr = np.sqrt(gamma_gas * R_gas * T_arr)
    mach_arr = V_arr / np.maximum(a_arr, 1e-6)

    L_over_m_arr = L_D * D_over_m_arr
    n_z = np.sqrt(D_over_m_arr ** 2 + L_over_m_arr ** 2) / g_arr

    q_dot_conv = np.array([stagnation_heat_sutton_graves(rho_arr[i], V_arr[i], RN, k=k_conv)
                           for i in range(len(t_arr))])
    q_dot_rad = np.array([radiative_heat_martin(RN, V_arr[i], rho_arr[i])
                          for i in range(len(t_arr))])
    q_total = q_dot_conv + q_dot_rad

    Q_cumul = np.cumsum(q_total * dt)
    range_km = np.cumsum(V_arr * np.cos(gamma_arr) * dt) / 1000.0

    return dict(
        t=t_arr,
        h_km=h_arr / 1000.0,
        V=V_arr,
        Mach=mach_arr,
        gamma_deg=np.degrees(gamma_arr),
        q_dot_conv=q_dot_conv,
        q_dot_rad=q_dot_rad,
        q_dot_total=q_total,
        n_z=n_z,
        q_dyn=q_dyn,
        p_atm=p_arr,
        Q_cumul=Q_cumul,
        D_over_m=D_over_m_arr,
        range_km=range_km,
        rho=rho_arr,
        planet=planet,
        stop_reason=stop_reason,
    )


# ----------------------------------------------------------------------
# CSV export utilities
# ----------------------------------------------------------------------
import csv
import os

def _csv_path(outdir, prefix, params_str):
    """Build CSV filepath: outdir/prefix_params.csv"""
    os.makedirs(outdir, exist_ok=True)
    safe = params_str.replace(" ", "_").replace("/", "-")
    return os.path.join(outdir, f"{prefix}_{safe}.csv")


def export_trajectory_csv(traj, outdir="outputs", Ve=0, BC=0, LD=0, gamma=0):
    """Save trajectory data to CSV."""
    fname = _csv_path(outdir, "trajectory",
                       f"Ve{Ve:.0f}_BC{BC:.0f}_LD{LD:.1f}_g{gamma:.1f}")
    keys = ["t", "h_km", "V", "Mach", "gamma_deg", "q_dot_conv", "q_dot_rad",
            "q_dot_total", "n_z", "q_dyn", "Q_cumul", "D_over_m", "range_km", "rho"]
    with open(fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(len(traj["t"])):
            w.writerow([traj[k][i] for k in keys])
    print(f"CSV saved: {fname}")
    return fname


def export_corridor_csv(corr, outdir="outputs", RN=0, BC=0, LD=0,
                         qtmax=0, nmax=0, qmax=0):
    """Save corridor data to CSV."""
    fname = _csv_path(outdir, "corridor",
                       f"RN{RN:.1f}_BC{BC:.0f}_LD{LD:.1f}_qt{qtmax:.0f}")
    keys = ["V", "h_upper", "h_lower", "h_heat", "h_load", "h_q",
            "corridor_width", "rho_upper", "rho_lower"]
    with open(fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(len(corr["V"])):
            row = []
            for k in keys:
                v = corr[k][i]
                row.append(f"{v:.6e}" if abs(v) < 1e-4 else f"{v:.4f}")
            w.writerow(row)
    print(f"CSV saved: {fname}")
    return fname


def export_stagnation_csv(outdir="outputs", alt=55, V=7800, RN=0.3):
    """Save stagnation point results to CSV."""
    rho = density_exponential(alt)
    fname = _csv_path(outdir, "stagnation",
                       f"alt{alt:.0f}_V{V:.0f}_RN{RN:.1f}")
    with open(fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "value", "unit"])
        w.writerow(["altitude_km", f"{alt:.1f}", "km"])
        w.writerow(["velocity_ms", f"{V:.1f}", "m/s"])
        w.writerow(["nose_radius_m", f"{RN:.3f}", "m"])
        w.writerow(["density_kgm3", f"{rho:.6e}", "kg/m³"])
        q_scott = stagnation_heat_scott(rho, V, RN)
        q_alt = stagnation_heat_alt(rho, V, RN)
        q_rad = radiative_heat_martin(RN, V, rho)
        w.writerow(["q_scott_Wcm2", f"{q_scott:.4f}", "W/cm²"])
        w.writerow(["q_alt_Wcm2", f"{q_alt:.4f}", "W/cm²"])
        w.writerow(["q_rad_Wcm2", f"{q_rad:.6f}", "W/cm²"])
    print(f"CSV saved: {fname}")
    return fname


def gliding_entry_analytical(Ve, L_D, bank_angle_deg=45.0, Vc=V0_CIRC):
    """
    Closed-form analytical equations for Gliding Entry & Footprint (TU Delft AE4870B Ch 4.3).
    
    Equations:
        Rf = -0.5 * Re * (L/D) * ln(1 - (Ve/Vc)²)                             [Eq. 4.29]
        t_flight = 0.5 * (Vc / g0) * (L/D) * ln((1 + Ve/Vc)/(1 - Ve/Vc))     [Eq. 4.32]
        Rfy_max = (pi²/24) * Re * (L/D)² * cos(sigma) * sin(sigma)           [Eq. 4.55]
    """
    x_E = min(float(Ve) / Vc, 0.999)
    LD = max(abs(float(L_D)), 1e-4)

    # Flight range (m)
    R_f = -0.5 * R_EARTH * LD * np.log(max(1.0 - x_E**2, 1e-10))

    # Analytical flight time (s)
    t_flight = 0.5 * (Vc / G0) * LD * np.log((1.0 + x_E) / max(1.0 - x_E, 1e-10))

    # Lateral cross-range footprint (m)
    sigma = np.radians(bank_angle_deg)
    R_fy = (np.pi**2 / 24.0) * R_EARTH * (LD**2) * np.cos(sigma) * np.sin(sigma)

    # Max cross-range occurs at bank angle sigma = 45 deg
    R_fy_max = (np.pi**2 / 24.0) * R_EARTH * (LD**2) * 0.5

    return {
        "Ve": Ve,
        "L_D": L_D,
        "bank_angle_deg": bank_angle_deg,
        "R_f_km": float(R_f / 1000.0),
        "t_flight_s": float(t_flight),
        "R_fy_km": float(R_fy / 1000.0),
        "R_fy_max_km": float(R_fy_max / 1000.0),
    }


def generate_openfoam_pato_boundary_conditions(traj, RN, eps=0.85, title="Earth Reentry Trajectory", dt_output=10.0):
    """
    Generate OpenFOAM / PATO boundaryConditions.txt file content with detailed engineering header comments.
    
    Parameters
    ----------
    dt_output : float
        Time increment between output rows (seconds). Default 10 s.
        The trajectory is downsampled to this interval.
    
    Columns generated:
      t(s)  p_total_w(Pa)  rhoeUeCH(kg/m2/s)  h_r(J/kg)  qRad(W/m2)  chemistryOn  lambda  T_inf(K)  p_dyn_w(Pa)
    """
    t = np.asarray(traj["t"])
    h_km = np.asarray(traj["h_km"])
    V = np.asarray(traj["V"])
    rho = np.asarray(traj["rho"])
    Mach = np.asarray(traj["Mach"])
    q_tot = np.asarray(traj["q_dot_total"])     # W/cm²
    q_conv = np.asarray(traj["q_dot_conv"])   # W/cm²
    q_rad = np.asarray(traj["q_dot_rad"])     # W/cm²
    Q_cumul = np.asarray(traj["Q_cumul"])     # J/cm²
    n_z = np.asarray(traj["n_z"])
    q_dyn = np.asarray(traj["q_dyn"])         # Pa

    # Peak events
    iq = int(np.argmax(q_tot))
    inz = int(np.argmax(n_z))
    iqdyn = int(np.argmax(q_dyn))

    # Header comments
    # --- Downsample trajectory to dt_output intervals ---
    t_full = t
    dt_out = max(float(dt_output), 0.1)
    if len(t_full) > 2:
        t_out_times = np.arange(0.0, t_full[-1], dt_out)
        # Always include the last point
        if t_out_times[-1] < t_full[-1]:
            t_out_times = np.append(t_out_times, t_full[-1])
        # Find nearest indices for each output time
        indices = []
        for tout in t_out_times:
            idx = int(np.argmin(np.abs(t_full - tout)))
            indices.append(idx)
        # Deduplicate while preserving order
        seen = set()
        unique_indices = []
        for idx in indices:
            if idx not in seen:
                seen.add(idx)
                unique_indices.append(idx)
        indices = unique_indices
    else:
        indices = list(range(len(t_full)))

    header = [
        "/*---------------------------------------------------------------------------*\\",
        "boundaryConditions — OpenFOAM / PATO Aerothermal Surface Boundary Conditions",
        f"Title: {title}",
        "-----------------------------------------------------------------------------",
        " 1. TRAJECTORY SUMMARY:",
        f"    - Altitude Corridor:  {h_km[0]:.1f} km  ->  {h_km[-1]:.1f} km",
        f"    - Velocity Range:     {V[0]:.1f} m/s ({Mach[0]:.1f} Mach)  ->  {V[-1]:.1f} m/s",
        f"    - Total Flight Time:  {t[-1]:.1f} s",
        f"    - Integrated Range:   {traj['range_km'][-1]:.1f} km",
        f"    - Output Time Step:   {dt_out:.1f} s ({len(indices)} rows)",
        "",
        " 2. VEHICLE SPECIFICATIONS:",
        f"    - Nose Radius (Rn):   {RN:.3f} m",
        f"    - Surface Emissivity: {eps:.2f}",
        "",
        " 3. PEAK AEROTHERMAL & FLIGHT EVENTS:",
        f"    - Peak Heat Flux:     {q_tot[iq]:.2f} W/cm² ({q_tot[iq]*1e4:.2e} W/m²) at t = {t[iq]:.1f} s (h = {h_km[iq]:.1f} km, V = {V[iq]:.0f} m/s, Mach = {Mach[iq]:.1f})",
        f"    - Peak Deceleration:  {n_z[inz]:.2f} g at t = {t[inz]:.1f} s (h = {h_km[inz]:.1f} km)",
        f"    - Peak Dyn Pressure:  {q_dyn[iqdyn]/1000:.2f} kPa at t = {t[iqdyn]:.1f} s",
        f"    - Total Heat Load Q:  {Q_cumul[-1]:.2f} J/cm² ({Q_cumul[-1]/100:.2f} MJ/m²)",
        "",
        " 4. ATMOSPHERIC REGIMES & CHEMISTRY:",
        "    - t = 0 to peak heat: Hypersonic shock dissociation regime (N2/O2 chemical reaction active, chemistryOn = 1)",
        "    - Post peak heat:     Supersonic / continuum deceleration phase",
        "    - Low altitude:       Subsonic terminal descent (chemistryOn = 0)",
        "\\*---------------------------------------------------------------------------*/",
        "/*",
        "t(s)       p_total_w(Pa)  rhoeUeCH(kg/m2/s)  h_r(J/kg)        qRad(W/m2)       chemistryOn  lambda  T_inf(K)   p_dyn_w(Pa)",
        "*/"
    ]

    rows_str = []
    rows_data = []

    for i in indices:
        ti = float(t[i])
        Vi = float(V[i])
        hi_m = float(h_km[i] * 1000.0)
        rho_i = float(rho[i])
        q_conv_Wm2 = float(q_conv[i] * 1e4)
        q_rad_Wm2 = float(q_rad[i] * 1e4)

        T_inf = isa_temperature(hi_m)
        P_inf = isa_pressure(hi_m)

        # Total stagnation pressure behind shock (Pa)
        P_total_w = P_inf + rho_i * Vi**2

        # Recovery enthalpy (J/kg): h_r = Cp * T_inf + 0.5 * V^2
        Cp = 1005.0
        h_r = Cp * T_inf + 0.5 * Vi**2

        # Wall enthalpy (J/kg) based on radiative equilibrium wall temp
        q_tot_Wm2 = max(q_tot[i] * 1e4, 1.0)
        Tw = (q_tot_Wm2 / (eps * SIGMA_SB)) ** 0.25
        hw = Cp * Tw

        # rhoeUeCH = q_conv / (h_r - h_w)  [kg/m²·s]
        dh = max(h_r - hw, 1e4)
        rhoeUeCH = q_conv_Wm2 / dh

        chem_flag = 1 if Vi > 1000.0 else 0
        lambda_val = 0.5
        p_dyn_w = 0.5 * rho_i * Vi**2

        line = f"{ti:10.4f}  {P_total_w:14.4f}  {rhoeUeCH:17.8e}  {h_r:16.1f}  {q_rad_Wm2:16.8f}     {chem_flag:1d}            {lambda_val:4.2f}  {T_inf:8.2f}  {p_dyn_w:12.4f}"
        rows_str.append(line)

        rows_data.append({
            "t": ti, "p_total_w": P_total_w, "rhoeUeCH": rhoeUeCH,
            "h_r": h_r, "qRad": q_rad_Wm2, "chemistryOn": chem_flag,
            "lambda": lambda_val, "T_inf": T_inf, "p_dyn_w": p_dyn_w
        })

    full_text = "\n".join(header + rows_str)
    return {
        "full_text": full_text,
        "header_comments": "\n".join(header),
        "rows_str": rows_str,
        "rows_data": rows_data,
        "n_rows": len(t),
    }


def solve_target_duration(target_t_s, V0, gamma0_deg, h0_km, RN, BC, L_D,
                          bank_angle_deg=0.0, vary="gamma", method="bounded"):
    """
    Find entry parameter value to match target reentry duration with feasibility bounds checking.

    Parameters
    ----------
    vary : str   Which parameter to optimise: 'gamma', 'BC', 'L_D', 'bank', 'V0', 'h0'
    """
    from scipy.optimize import minimize_scalar

    if V0 > 25000.0:
        raise ValueError(f"Entry speed V0 = {V0:.0f} m/s exceeds physical entry velocity limit (25,000 m/s)")

    target_t = float(target_t_s)

    def _sim_time(val):
        kw = dict(V0=V0, gamma0_deg=gamma0_deg, h0_km=h0_km,
                  RN=RN, BC=BC, L_D=L_D, bank_angle_deg=bank_angle_deg,
                  t_max=max(target_t * 3, 3600.0), dt=1.0)
        if vary == "gamma":
            kw["gamma0_deg"] = val
        elif vary == "BC":
            kw["BC"] = max(val, 5.0)
        elif vary == "L_D":
            kw["L_D"] = max(val, 0.0)
        elif vary == "bank":
            kw["bank_angle_deg"] = val
        elif vary == "V0":
            kw["V0"] = max(val, 100.0)
        elif vary == "h0":
            kw["h0_km"] = max(val, 40.0)
        else:
            raise ValueError(f"Unknown vary parameter: {vary}")

        try:
            res = trajectory_integrate(**kw)
            return res["t"][-1]
        except Exception:
            return 0.0

    # Physical search bounds for each parameter
    bounds = {
        "gamma": (-0.5, -40.0),
        "BC":    (10.0, 2000.0),
        "L_D":   (0.0, 3.0),
        "bank":  (0.0, 85.0),
        "V0":    (3000.0, 12000.0),
        "h0":    (60.0, 200.0),
    }
    lo, hi = bounds[vary]
    if lo > hi:
        lo, hi = hi, lo

    # Evaluate grid over bounds to determine physical min/max achievable times
    n_sample = 30
    sample_vals = np.linspace(lo, hi, n_sample)
    sample_times = np.array([_sim_time(v) for v in sample_vals])

    achievable_min = float(np.min(sample_times[sample_times > 0])) if np.any(sample_times > 0) else 0.0
    achievable_max = float(np.max(sample_times))

    is_attainable = (achievable_min <= target_t <= achievable_max)

    # Use minimize_scalar bounded or fine grid search
    def _cost(v):
        return ( _sim_time(v) - target_t ) ** 2

    if method == "bounded":
        res_opt = minimize_scalar(_cost, bounds=(lo, hi), method=method, options={'xatol': 1e-3})
    else:
        # Golden or brent
        res_opt = minimize_scalar(_cost, bracket=(lo, hi), method=method)
    best_val = float(res_opt.x)

    # Run final fine trajectory at best value
    final_kw = dict(V0=V0, gamma0_deg=gamma0_deg, h0_km=h0_km,
                    RN=RN, BC=BC, L_D=L_D, bank_angle_deg=bank_angle_deg,
                    t_max=max(target_t * 3, 3600.0), dt=0.5)
    if vary == "gamma":
        final_kw["gamma0_deg"] = best_val
    elif vary == "BC":
        final_kw["BC"] = max(best_val, 5.0)
    elif vary == "L_D":
        final_kw["L_D"] = max(best_val, 0.0)
    elif vary == "bank":
        final_kw["bank_angle_deg"] = best_val
    elif vary == "V0":
        final_kw["V0"] = max(best_val, 100.0)
    elif vary == "h0":
        final_kw["h0_km"] = max(best_val, 40.0)

    best_traj = trajectory_integrate(**final_kw)
    achieved_t = float(best_traj["t"][-1])
    error_s = float(abs(achieved_t - target_t))

    return {
        "target_t_s": target_t,
        "vary": vary,
        "matched_value": best_val,
        "achieved_t_s": achieved_t,
        "error_s": error_s,
        "is_attainable": is_attainable,
        "achievable_min_s": achievable_min,
        "achievable_max_s": achievable_max,
        "trajectory": best_traj,
        "final_params": final_kw,
    }


def solve_target_heating(target_value, target_type, V0, gamma0_deg, h0_km, RN, BC, L_D,
                           bank_angle_deg=0.0, vary="gamma", method="bounded"):
    """
    Find entry parameter value to match a target heating condition with feasibility bounds checking.

    Parameters
    ----------
    target_value : float Target value (W/cm^2 for flux, J/cm^2 for load)
    target_type : str    'q_max' for peak heat flux, 'Q_total' for total heat load
    vary : str           Which parameter to optimise: 'gamma', 'BC', 'L_D', 'bank', 'V0', 'h0'
    """
    from scipy.optimize import minimize_scalar

    if V0 > 25000.0:
        raise ValueError(f"Entry speed V0 = {V0:.0f} m/s exceeds physical entry velocity limit (25,000 m/s)")

    target_val = float(target_value)

    def _sim_heat(val):
        kw = dict(V0=V0, gamma0_deg=gamma0_deg, h0_km=h0_km,
                  RN=RN, BC=BC, L_D=L_D, bank_angle_deg=bank_angle_deg,
                  t_max=3600.0, dt=1.0)
        if vary == "gamma":
            kw["gamma0_deg"] = val
        elif vary == "BC":
            kw["BC"] = max(val, 5.0)
        elif vary == "L_D":
            kw["L_D"] = max(val, 0.0)
        elif vary == "bank":
            kw["bank_angle_deg"] = val
        elif vary == "V0":
            kw["V0"] = max(val, 100.0)
        elif vary == "h0":
            kw["h0_km"] = max(val, 40.0)
        else:
            raise ValueError(f"Unknown vary parameter: {vary}")

        try:
            res = trajectory_integrate(**kw)
            if target_type == 'q_max':
                return float(np.max(res["q_dot_total"]))
            elif target_type == 'Q_total':
                return float(res["Q_cumul"][-1])
            else:
                return 0.0
        except Exception:
            return 0.0

    # Physical search bounds for each parameter
    bounds = {
        "gamma": (-0.5, -40.0),
        "BC":    (10.0, 2000.0),
        "L_D":   (0.0, 3.0),
        "bank":  (0.0, 85.0),
        "V0":    (3000.0, 12000.0),
        "h0":    (60.0, 200.0),
    }
    lo, hi = bounds[vary]
    if lo > hi:
        lo, hi = hi, lo

    # Evaluate grid over bounds to determine physical min/max achievable
    n_sample = 30
    sample_vals = np.linspace(lo, hi, n_sample)
    sample_heats = np.array([_sim_heat(v) for v in sample_vals])

    valid = sample_heats > 0
    if np.any(valid):
        achievable_min = float(np.min(sample_heats[valid]))
        achievable_max = float(np.max(sample_heats[valid]))
    else:
        achievable_min = 0.0
        achievable_max = 0.0

    is_attainable = (achievable_min <= target_val <= achievable_max)

    def _cost(v):
        h_val = _sim_heat(v)
        if h_val <= 0:
            return 1e9 # Penalty for failed runs
        return (h_val - target_val) ** 2

    if method == "bounded":
        res_opt = minimize_scalar(_cost, bounds=(lo, hi), method=method, options={'xatol': 1e-3})
    else:
        res_opt = minimize_scalar(_cost, bracket=(lo, hi), method=method)
    best_val = float(res_opt.x)

    # Run final fine trajectory at best value
    final_kw = dict(V0=V0, gamma0_deg=gamma0_deg, h0_km=h0_km,
                    RN=RN, BC=BC, L_D=L_D, bank_angle_deg=bank_angle_deg,
                    t_max=3600.0, dt=0.5)
    if vary == "gamma":
        final_kw["gamma0_deg"] = best_val
    elif vary == "BC":
        final_kw["BC"] = max(best_val, 5.0)
    elif vary == "L_D":
        final_kw["L_D"] = max(best_val, 0.0)
    elif vary == "bank":
        final_kw["bank_angle_deg"] = best_val
    elif vary == "V0":
        final_kw["V0"] = max(best_val, 100.0)
    elif vary == "h0":
        final_kw["h0_km"] = max(best_val, 40.0)

    best_traj = trajectory_integrate(**final_kw)
    if target_type == 'q_max':
        achieved_val = float(np.max(best_traj["q_dot_total"]))
    else:
        achieved_val = float(best_traj["Q_cumul"][-1])
        
    error_val = float(abs(achieved_val - target_val))

    return {
        "target_val": target_val,
        "target_type": target_type,
        "vary": vary,
        "matched_value": best_val,
        "achieved_val": achieved_val,
        "error_val": error_val,
        "is_attainable": is_attainable,
        "achievable_min": achievable_min,
        "achievable_max": achievable_max,
        "trajectory": best_traj,
        "final_params": final_kw,
    }


def monte_carlo_trajectory(V0, gamma0_deg, h0_km, RN, BC, L_D,
                            n_runs=50, gamma_std=0.3, BC_std=10.0, LD_std=0.02, distribution="normal"):
    """
    Run Monte Carlo stochastic trajectory batch and return 3-sigma dispersional statistics.
    """
    t_finals = []
    q_maxs = []
    Q_totals = []
    n_maxs = []
    ranges = []

    np.random.seed(42)
    for _ in range(int(n_runs)):
        if distribution == "uniform":
            g_i = np.random.uniform(gamma0_deg - gamma_std*1.732, gamma0_deg + gamma_std*1.732)
            BC_i = max(np.random.uniform(BC - BC_std*1.732, BC + BC_std*1.732), 10.0)
            LD_i = max(np.random.uniform(L_D - LD_std*1.732, L_D + LD_std*1.732), 0.0)
        else:
            g_i = np.random.normal(gamma0_deg, gamma_std)
            BC_i = max(np.random.normal(BC, BC_std), 10.0)
            LD_i = max(np.random.normal(L_D, LD_std), 0.0)

        try:
            res = trajectory_integrate(V0=V0, gamma0_deg=g_i, h0_km=h0_km,
                                       RN=RN, BC=BC_i, L_D=LD_i, t_max=4800.0, dt=1.5)
            t_finals.append(res["t"][-1])
            q_maxs.append(np.max(res["q_dot_total"]))
            Q_totals.append(res["Q_cumul"][-1])
            n_maxs.append(np.max(res["n_z"]))
            ranges.append(res["range_km"][-1])
        except Exception:
            pass

    t_finals = np.array(t_finals)
    q_maxs = np.array(q_maxs)
    Q_totals = np.array(Q_totals)
    n_maxs = np.array(n_maxs)
    ranges = np.array(ranges)

    return {
        "n_runs": len(t_finals),
        "flight_time": {"mean": float(np.mean(t_finals)), "std": float(np.std(t_finals)), "min": float(np.min(t_finals)), "max": float(np.max(t_finals))},
        "peak_q": {"mean": float(np.mean(q_maxs)), "std": float(np.std(q_maxs)), "min": float(np.min(q_maxs)), "max": float(np.max(q_maxs))},
        "total_Q": {"mean": float(np.mean(Q_totals)), "std": float(np.std(Q_totals)), "min": float(np.min(Q_totals)), "max": float(np.max(Q_totals))},
        "peak_g": {"mean": float(np.mean(n_maxs)), "std": float(np.std(n_maxs)), "min": float(np.min(n_maxs)), "max": float(np.max(n_maxs))},
        "range_km": {"mean": float(np.mean(ranges)), "std": float(np.std(ranges)), "min": float(np.min(ranges)), "max": float(np.max(ranges))},
    }


def solve_1d_transient_tps(t_arr, q_total_Wm2, d_m=0.04, rho_tps=300.0, cp_tps=1200.0, k_tps=0.15, T_initial=300.0, eps=0.85):
    """
    Solves 1D transient heat conduction PDE across TPS thickness d_m:
        rho * cp * dT/dt = k * d^2T/dx^2

    Boundary Conditions:
        - Surface (x=0): Newton-Raphson non-linear re-radiation energy balance:
          q_inc = eps * sigma * T_surf^4 + k * (T_surf - T_1)/dx + (rho * cp * dx / (2*dt)) * (T_surf - T_surf_old)
        - Backwall (x=d): adiabatic (dT/dx = 0)
    """
    sigma_sb = 5.670374e-8
    N = 25  # spatial nodes across thickness
    dx = d_m / max(N - 1, 1)
    
    T = np.full(N, float(T_initial), dtype=float)
    t_len = len(t_arr)
    
    T_surface_hist = np.zeros(t_len)
    T_backwall_hist = np.zeros(t_len)
    
    alpha = k_tps / max(rho_tps * cp_tps, 1e-6) # thermal diffusivity [m^2/s]
    
    for i in range(t_len):
        dt_macro = float(t_arr[i] - t_arr[i-1]) if i > 0 else 0.5
        dt_macro = max(dt_macro, 0.01)
        q_inc = max(float(q_total_Wm2[i]), 0.0)
        
        # Sub-stepping for numerical stability: Fo = alpha * dt / dx^2 <= 0.4
        dt_sub_max = 0.4 * (dx**2) / max(alpha, 1e-10)
        n_sub = max(int(np.ceil(dt_macro / dt_sub_max)), 1)
        dt_step = dt_macro / n_sub
        
        c_cap = (rho_tps * cp_tps * dx) / (2.0 * dt_step)
        k_dx = k_tps / dx
        
        for _ in range(n_sub):
            T_old = np.copy(T)
            T_new = np.copy(T)
            
            # Newton-Raphson solve for surface node T_new[0]
            T0_guess = T_old[0]
            for _iter in range(8):
                f_val = eps * sigma_sb * (T0_guess**4) + k_dx * (T0_guess - T_old[1]) + c_cap * (T0_guess - T_old[0]) - q_inc
                df_val = 4.0 * eps * sigma_sb * (T0_guess**3) + k_dx + c_cap
                T0_guess = T0_guess - f_val / max(df_val, 1e-6)
                T0_guess = max(T0_guess, 100.0)
                
            T_new[0] = T0_guess
            
            # Interior nodes
            fo = alpha * dt_step / (dx**2)
            for j in range(1, N - 1):
                T_new[j] = T_old[j] + fo * (T_old[j+1] - 2.0 * T_old[j] + T_old[j-1])
                
            # Adiabatic backwall node (x=d)
            T_new[N-1] = T_old[N-1] + 2.0 * fo * (T_old[N-2] - T_old[N-1])
            T = T_new
            
        T_surface_hist[i] = T[0]
        T_backwall_hist[i] = T[N-1]
        
    i_max_bw = int(np.argmax(T_backwall_hist))
    return {
        "t": t_arr,
        "T_surface": T_surface_hist,
        "T_backwall": T_backwall_hist,
        "T_backwall_max": float(np.max(T_backwall_hist)),
        "t_soak_max_s": float(t_arr[i_max_bw]),
        "T_final_profile": T
    }


# ----------------------------------------------------------------------
# 3-Layer Composite TPS Stack (Glass Coating + 90% Porous Tile + Al-7xxx Metal)
# ----------------------------------------------------------------------
T_K_TABLE = np.array([310.0, 469.0, 675.0, 897.0, 1082.0, 1273.0])
K_TILE_TABLE = np.array([0.064, 0.086, 0.112, 0.148, 0.198, 0.308]) # W/m.K

T_CP_TABLE = np.array([323.0, 353.0, 373.0, 473.0, 675.0, 859.0])
CP_TILE_TABLE = np.array([930.0, 940.0, 960.0, 1010.0, 1150.0, 1270.0]) # J/kg.K

T_EPS_TABLE = np.array([303.0, 673.0, 1173.0])
EPS_TILE_TABLE = np.array([0.80, 0.83, 0.87])

def get_tile_k(T):
    return float(np.interp(T, T_K_TABLE, K_TILE_TABLE))

def get_tile_cp(T):
    return float(np.interp(T, T_CP_TABLE, CP_TILE_TABLE))

def get_tile_eps(T):
    return float(np.interp(T, T_EPS_TABLE, EPS_TILE_TABLE))
SUBSTRUCTURE_MATERIALS = {
    "Al-7xxx (Default)": {"rho": 2810.0, "cp": 960.0,  "T_limit": 450.0},
    "Al-2024":           {"rho": 2780.0, "cp": 875.0,  "T_limit": 450.0},
    "Al-6061":           {"rho": 2700.0, "cp": 896.0,  "T_limit": 450.0},
    "Inconel 718":       {"rho": 8190.0, "cp": 435.0,  "T_limit": 950.0},
    "SS-304":            {"rho": 8000.0, "cp": 500.0,  "T_limit": 1100.0},
    "Titanium (Ti-6Al-4V)":{"rho": 4430.0,"cp": 526.0, "T_limit": 650.0},
}


def solve_1d_transient_tps_multilayer(
    t_arr, q_total_Wm2,
    d_tile_m=0.04,
    d_glass_m=0.0003, # 300 um Reaction Glass Coating (SiO2-B2O3)
    d_metal_m=0.002,  # 2.0 mm Al-7xxx Substructure
    rho_glass=2200.0, k_glass=1.4, cp_glass=900.0,
    rho_tile=250.0, k_tile=0.15, cp_tile=1200.0,   # 0.25 g/cc (90% Porous Silica Tile)
    rho_metal=2810.0, k_metal=150.0, cp_metal=875.0,
    T_initial=300.0,
    use_tdependent=True
):
    """
    3-Layer Composite Transient Heat Conduction Solver:
    - Outer Coating: Reaction Glass Coating (SiO2-B2O3, 300 um)
    - Core Tile: 90% Porous Silica Tile (rho = 250 kg/m3, T-dependent k, cp, eps)
    - Substructure: Aluminum Alloy Al-7xxx (2.0 mm)
    """
    sigma_sb = 5.670374e-8
    N = 25 # Spatial nodes across Porous Silica Tile
    dx = d_tile_m / max(N - 1, 1)

    T = np.full(N, float(T_initial), dtype=float)
    T_metal = float(T_initial)
    t_len = len(t_arr)

    T_surface_hist = np.zeros(t_len)
    T_center_hist = np.zeros(t_len)
    T_backwall_hist = np.zeros(t_len)
    k_surf_hist = np.zeros(t_len)
    k_center_hist = np.zeros(t_len)
    k_back_hist = np.zeros(t_len)
    q_rerad_hist = np.zeros(t_len)
    q_cond_hist = np.zeros(t_len)
    T_all_hist = np.zeros((t_len, N))
    x_nodes = np.linspace(0, d_tile_m, N)

    R_glass = d_glass_m / k_glass  # Glass thermal resistance
    C_metal = rho_metal * cp_metal * d_metal_m # Metal heat capacity per m2

    for i in range(t_len):
        dt_macro = float(t_arr[i] - t_arr[i-1]) if i > 0 else 0.5
        dt_macro = max(dt_macro, 0.01)
        q_inc = max(float(q_total_Wm2[i]), 0.0)

        k_n = np.zeros(N)
        cp_n = np.zeros(N)
        for n in range(N):
            if use_tdependent:
                k_n[n] = get_tile_k(T[n])
                cp_n[n] = get_tile_cp(T[n])
            else:
                k_n[n] = float(k_tile)
                cp_n[n] = float(cp_tile)

        alpha_max = np.max(k_n / (rho_tile * cp_n))
        dt_sub_max = 0.4 * (dx**2) / max(alpha_max, 1e-10)
        n_sub = max(int(np.ceil(dt_macro / dt_sub_max)), 1)
        dt_step = dt_macro / n_sub

        for _ in range(n_sub):
            for n in range(N):
                if use_tdependent:
                    k_n[n] = get_tile_k(T[n])
                    cp_n[n] = get_tile_cp(T[n])

            T_old = np.copy(T)
            T_new = np.copy(T)
            T_m_old = T_metal

            eps_surf = get_tile_eps(T_old[0]) if use_tdependent else 0.85

            U_surf = 1.0 / (R_glass + dx / (2.0 * k_n[0]))
            c_cap0 = (rho_tile * cp_n[0] * dx) / (2.0 * dt_step)

            # Newton-Raphson for outer glass surface temperature T0
            T0_guess = T_old[0]
            for _iter in range(8):
                f_val = eps_surf * sigma_sb * (T0_guess**4) + U_surf * (T0_guess - T_old[0]) - q_inc
                df_val = 4.0 * eps_surf * sigma_sb * (T0_guess**3) + U_surf
                T0_guess = T0_guess - f_val / max(df_val, 1e-6)
                T0_guess = max(T0_guess, 100.0)

            T_surf_actual = T0_guess
            q_cond_in = U_surf * (T_surf_actual - T_old[0])

            # Node 0 update
            T_new[0] = T_old[0] + (dt_step / c_cap0) * (q_cond_in - k_n[0] * (T_old[0] - T_old[1]) / dx)

            # Interior nodes
            for j in range(1, N - 1):
                alpha_j = k_n[j] / (rho_tile * cp_n[j])
                fo = alpha_j * dt_step / (dx**2)
                T_new[j] = T_old[j] + fo * (T_old[j+1] - 2.0 * T_old[j] + T_old[j-1])

            # Bottom tile node coupled to Al-7xxx substructure
            U_bot = k_n[-1] / dx
            q_to_metal = U_bot * (T_old[-1] - T_m_old)
            T_new[-1] = T_old[-1] + (dt_step / c_cap0) * (k_n[-1] * (T_old[-2] - T_old[-1]) / dx - q_to_metal)

            # Al-7xxx Metal Substructure update
            T_metal = T_m_old + (dt_step / C_metal) * q_to_metal
            T = T_new

        T_surface_hist[i] = T_surf_actual
        T_center_hist[i] = T[N // 2]
        T_backwall_hist[i] = T_metal
        T_all_hist[i, :] = T
        k_surf_hist[i] = k_n[0]
        k_center_hist[i] = k_n[N // 2]
        k_back_hist[i] = k_n[-1]
        
        q_rerad_hist[i] = (eps_surf * sigma_sb * (T_surf_actual**4)) / 1e4 # W/cm2
        q_cond_hist[i] = q_cond_in / 1e4 # W/cm2

    i_max_bw = int(np.argmax(T_backwall_hist))
    x_nodes_full = np.concatenate([
        np.array([-d_glass_m]),
        x_nodes,
        np.array([d_tile_m + d_metal_m])
    ])
    T_soak_full = np.concatenate([
        np.array([T_surface_hist[i_max_bw]]),
        T_all_hist[i_max_bw],
        np.array([T_backwall_hist[i_max_bw]])
    ])

    return {
        "t": t_arr,
        "T_surface": T_surface_hist,
        "T_center": T_center_hist,
        "T_backwall": T_backwall_hist,
        "T_backwall_max": float(np.max(T_backwall_hist)),
        "t_soak_max_s": float(t_arr[i_max_bw]),
        "T_final_profile": T,
        "x_nodes": x_nodes,
        "x_nodes_full": x_nodes_full,
        "T_all_hist": T_all_hist,
        "T_soak_profile": T_all_hist[i_max_bw],
        "T_soak_profile_full": T_soak_full,
        "k_surf_hist": k_surf_hist,
        "k_center_hist": k_center_hist,
        "k_back_hist": k_back_hist,
        "q_rerad": q_rerad_hist,
        "q_cond": q_cond_hist
    }


def solve_required_tps_thickness(t_arr, q_total_Wm2, T_backwall_limit=450.0,
                                   rho_tps=250.0, cp_tps=1200.0, k_tps=0.15,
                                   d_metal=0.002, rho_metal=2810.0, cp_metal=960.0,
                                   T_initial=300.0, eps=0.85,
                                   d_min=0.005, d_max=0.25, use_multilayer=True, 
                                   use_tdependent=True, method="brentq"):
    """
    Sizes the minimum required TPS tile thickness d_m (meters) using Brent's method root-finding
    such that peak transient backwall temperature T_backwall_max <= T_backwall_limit.
    """
    def _eval(d_val):
        if use_multilayer:
            return solve_1d_transient_tps_multilayer(t_arr, q_total_Wm2, d_tile_m=d_val, rho_tile=rho_tps, 
                                                     k_tile=k_tps, cp_tile=cp_tps,
                                                     d_metal_m=d_metal, rho_metal=rho_metal, cp_metal=cp_metal,
                                                     T_initial=T_initial, use_tdependent=use_tdependent)
        else:
            return solve_1d_transient_tps(t_arr, q_total_Wm2, d_m=d_val, rho_tps=rho_tps, cp_tps=cp_tps, k_tps=k_tps, T_initial=T_initial, eps=eps)

    from scipy.optimize import root_scalar

    def _obj(d_val):
        sol = _eval(d_val)
        return sol["T_backwall_max"] - T_backwall_limit

    sol_min = _eval(d_min)
    sol_max = _eval(d_max)

    if sol_min["T_backwall_max"] <= T_backwall_limit:
        d_req = d_min
        sol_opt = sol_min
        is_attainable = True
    elif sol_max["T_backwall_max"] > T_backwall_limit:
        d_req = d_max
        sol_opt = sol_max
        is_attainable = False
    else:
        is_attainable = True
        if method == "secant":
            res = root_scalar(_obj, method=method, x0=d_min, x1=d_max, xtol=1e-3)
        elif method == "bisect":
            res = root_scalar(_obj, method=method, bracket=[d_min, d_max], xtol=1e-3)
        else:
            res = root_scalar(_obj, method="brentq", bracket=[d_min, d_max], xtol=1e-3)
        
        d_req = float(res.root)
        sol_opt = _eval(d_req)

    mass_areal_g_cm2 = (rho_tps * d_req) * 0.1
    return {
        "is_attainable": is_attainable,
        "d_req_m": float(d_req),
        "d_req_cm": float(d_req * 100.0),
        "d_req_mm": float(d_req * 1000.0),
        "d_max_cm": float(d_max * 100.0),
        "d_max_mm": float(d_max * 1000.0),
        "mass_areal_g_cm2": float(mass_areal_g_cm2),
        "T_backwall_max": sol_opt["T_backwall_max"],
        "t_soak_max_s": sol_opt["t_soak_max_s"],
        "T_backwall_limit": float(T_backwall_limit),
        "sol_tps": sol_opt
    }


def solve_target_surface_temp(V0, h0_km, RN, BC, L_D, T_surf_target_K=1533.15,
                              eps=0.85, bank_angle_deg=0.0, planet='Earth', method='bisect',
                              vary='gamma', gamma_min=-30.0, gamma_max=-0.5,
                              RN_min=0.1, RN_max=10.0, BC_min=10.0, BC_max=1000.0):
    """
    Solves for entry parameter (gamma0, RN, or BC) such that the peak 
    equilibrium surface temperature T_surf_max equals T_surf_target_K.
    """
    sigma_sb = 5.670374e-8

    def _eval(p_val, current_vary):
        if current_vary == 'gamma':
            g_val, r_val, b_val = p_val, RN, BC
        elif current_vary == 'RN':
            g_val = gamma_min if gamma_min > -10.0 else -6.0
            r_val, b_val = p_val, BC
        elif current_vary == 'BC':
            g_val = gamma_min if gamma_min > -10.0 else -6.0
            r_val, b_val = RN, p_val
        else:
            g_val, r_val, b_val = p_val, RN, BC

        traj = trajectory_integrate(V0=V0, gamma0_deg=g_val, h0_km=h0_km, RN=r_val, BC=b_val, L_D=L_D,
                                    bank_angle_deg=bank_angle_deg, planet=planet)
        q_peak_Wm2 = np.max(traj["q_dot_total"]) * 1e4
        T_surf_peak = (max(q_peak_Wm2, 0.0) / (max(eps, 0.01) * sigma_sb)) ** 0.25
        return traj, float(T_surf_peak), float(g_val), float(r_val), float(b_val)

    from scipy.optimize import root_scalar

    cur_vary = vary
    if cur_vary == 'RN':
        p_min, p_max = RN_min, RN_max
    elif cur_vary == 'BC':
        p_min, p_max = BC_min, BC_max
    else:
        p_min, p_max = gamma_min, gamma_max

    def _obj(p_val):
        _, T_surf_peak, _, _, _ = _eval(p_val, cur_vary)
        return T_surf_peak - T_surf_target_K

    traj_min, T_min, g_min, r_min, b_min = _eval(p_min, cur_vary)
    traj_max, T_max, g_max, r_max, b_max = _eval(p_max, cur_vary)

    if T_min > T_max:
        is_attainable = (T_max <= T_surf_target_K <= T_min)
    else:
        is_attainable = (T_min <= T_surf_target_K <= T_max)

    if not is_attainable and cur_vary == 'gamma':
        # Auto-fallback: search across Nose Radius RN
        cur_vary = 'RN'
        p_min, p_max = RN_min, RN_max
        traj_min, T_min, g_min, r_min, b_min = _eval(p_min, cur_vary)
        traj_max, T_max, g_max, r_max, b_max = _eval(p_max, cur_vary)
        if T_min > T_max:
            is_attainable = (T_max <= T_surf_target_K <= T_min)
        else:
            is_attainable = (T_min <= T_surf_target_K <= T_max)

    if not is_attainable:
        p_opt = p_min if abs(T_min - T_surf_target_K) < abs(T_max - T_surf_target_K) else p_max
        traj_opt, T_opt, g_opt, r_opt, b_opt = _eval(p_opt, cur_vary)
    else:
        if method == "secant":
            res = root_scalar(_obj, method=method, x0=p_min, x1=p_max, xtol=1e-3)
        else:
            bracket = [p_min, p_max] if T_min <= T_max else [p_max, p_min]
            res = root_scalar(_obj, method=method, bracket=bracket, xtol=1e-3)
        p_opt = float(res.root)
        traj_opt, T_opt, g_opt, r_opt, b_opt = _eval(p_opt, cur_vary)

    return {
        "is_attainable": is_attainable,
        "vary_param": cur_vary,
        "matched_value": float(p_opt),
        "gamma_opt_deg": float(g_opt),
        "RN_opt_m": float(r_opt),
        "BC_opt_kgm2": float(b_opt),
        "T_surf_peak_K": float(T_opt),
        "T_surf_peak_C": float(T_opt - 273.15),
        "T_surf_target_K": float(T_surf_target_K),
        "T_surf_target_C": float(T_surf_target_K - 273.15),
        "trajectory": traj_opt
    }
