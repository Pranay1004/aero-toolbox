"""
FastAPI backend for Reentry Aerothermodynamics Calculator Web Edition.
Mirrors all 10 GUI tabs as JSON API endpoints.
"""
from __future__ import annotations
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import reentry_lib as R

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Optional, List
import json

app = FastAPI(title="Reentry Aerothermodynamics Calculator")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

def _ndarray_to_list(arr):
    if isinstance(arr, np.ndarray):
        return arr.tolist()
    if isinstance(arr, (np.floating, np.integer)):
        return float(arr)
    return arr

def _traj_to_json(traj):
    return {k: _ndarray_to_list(v) for k, v in traj.items()}

# ============================================================
# Tab 1: Stagnation Heating (full trajectory)
# ============================================================
class StagnationRequest(BaseModel):
    planet: str = "Earth"
    h0_km: float = 120.0
    V0: float = 11000.0
    gamma0: float = -20.0
    RN: float = 1.0
    BC: float = 400.0
    LD: float = 0.0
    eps: float = 0.85
    t_max: float = 0.0

@app.post("/api/stagnation")
def calc_stagnation(req: StagnationRequest):
    try:
        traj = R.trajectory_integrate(
            V0=req.V0, gamma0_deg=req.gamma0, h0_km=req.h0_km,
            RN=req.RN, BC=req.BC, L_D=req.LD, t_max=req.t_max, dt=0.5,
            planet=req.planet
        )
        h_km = np.array(traj["h_km"])
        V = np.array(traj["V"])
        rho = np.array(traj["rho"])
        t = np.array(traj["t"])

        q_fr = np.array([R.stagnation_heat_fay_riddell(rho[i], V[i], req.RN) for i in range(len(t))])
        q_ae = np.array([R.stagnation_heat_allen_eggers(rho[i], V[i], req.RN) for i in range(len(t))])
        T_inf = np.array([R.isa_temperature(h*1000.0, planet=req.planet) for h in h_km])
        p_inf = np.array([R.isa_pressure(h*1000.0, planet=req.planet) for h in h_km])

        q_full = np.zeros(len(t))
        Tw_full = np.zeros(len(t))
        hw_full = np.zeros(len(t))
        He_full = np.zeros(len(t))
        for i in range(len(t)):
            q_full[i], Tw_full[i], hw_full[i], He_full[i] = R.stagnation_heat_fay_riddell_full(
                rho[i], V[i], T_inf[i], p_inf[i], req.RN, req.eps
            )

        eps_safe = max(req.eps, 1e-6)
        Tw_fr = R.radiative_equilibrium_temp(q_fr * 1e4, eps_safe)
        Tw_ae = R.radiative_equilibrium_temp(q_ae * 1e4, eps_safe)

        # Velocity sweep bar chart
        V_bar = np.array([2000, 6000, 8000, 10000, 12000])
        q_fr_bar = []
        q_ae_bar = []
        for vb in V_bar:
            V_path, rho_path = R.allen_eggers_path(vb, req.BC, req.gamma0, L_D=req.LD)
            q_fr_path = np.array([R.stagnation_heat_fay_riddell(rho_path[j], V_path[j], req.RN) for j in range(len(V_path))])
            q_ae_path = np.array([R.stagnation_heat_allen_eggers(rho_path[j], V_path[j], req.RN) for j in range(len(V_path))])
            q_fr_bar.append(float(np.max(q_fr_path)))
            q_ae_bar.append(float(np.max(q_ae_path)))

        iq_fr = int(np.argmax(q_fr))
        iq_ae = int(np.argmax(q_ae))
        iq_full = int(np.argmax(q_full))
        iTw_full = int(np.argmax(Tw_full))

        report = {
            "stop_reason": traj["stop_reason"],
            "flight_time": float(t[-1]),
            "end_h": float(h_km[-1]),
            "end_V": float(V[-1]),
            "end_Mach": float(traj["Mach"][-1]),
            "max_mach": float(np.max(traj["Mach"])),
            "fr_peak_q": float(q_fr[iq_fr]),
            "fr_peak_h": float(h_km[iq_fr]),
            "fr_peak_V": float(V[iq_fr]),
            "fr_peak_t": float(t[iq_fr]),
            "ae_peak_q": float(q_ae[iq_ae]),
            "full_peak_q": float(q_full[iq_full]),
            "full_peak_Tw": float(Tw_full[iTw_full]),
        }

        return {
            "report": report,
            "t": _ndarray_to_list(t),
            "h_km": _ndarray_to_list(h_km),
            "V": _ndarray_to_list(V),
            "Mach": _ndarray_to_list(traj["Mach"]),
            "q_fr": _ndarray_to_list(q_fr),
            "q_ae": _ndarray_to_list(q_ae),
            "q_full": _ndarray_to_list(q_full),
            "Tw_fr": _ndarray_to_list(Tw_fr),
            "Tw_ae": _ndarray_to_list(Tw_ae),
            "Tw_full": _ndarray_to_list(Tw_full),
            "n_z": _ndarray_to_list(traj["n_z"]),
            "V_bar": V_bar.tolist(),
            "q_fr_bar": q_fr_bar,
            "q_ae_bar": q_ae_bar,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 2: BC Calculator
# ============================================================
class BCCalcRequest(BaseModel):
    mass: Optional[float] = None
    cd: Optional[float] = None
    area: Optional[float] = None
    bc: Optional[float] = None

@app.post("/api/bc")
def calc_bc(req: BCCalcRequest):
    try:
        if req.mass is not None and req.cd is not None and req.area is not None:
            if req.mass <= 0 or req.cd <= 0 or req.area <= 0:
                raise ValueError("All values must be positive")
            bc = req.mass / (req.cd * req.area)
            return {"bc": bc, "unit": "kg/m²"}
        elif req.bc is not None:
            results = {}
            if req.cd is not None and req.area is not None:
                results["mass"] = req.bc * req.cd * req.area
                results["mode"] = "mass"
            elif req.mass is not None and req.area is not None:
                results["cd"] = req.mass / (req.bc * req.area)
                results["mode"] = "cd"
            elif req.mass is not None and req.cd is not None:
                results["area"] = req.mass / (req.bc * req.cd)
                results["mode"] = "area"
            return results
        else:
            raise ValueError("Provide either (mass, cd, area) or (bc + two others)")
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 3: Ballistic Reentry
# ============================================================
class BallisticRequest(BaseModel):
    Ve: float = 11000.0
    BC: float = 400.0
    gamma: float = 20.0
    RN: float = 1.0

@app.post("/api/ballistic")
def calc_ballistic(req: BallisticRequest):
    try:
        Q = R.ballistic_heat_load(req.Ve, req.BC, req.gamma, RN=req.RN)
        qmax = R.ballistic_max_heat_rate(req.Ve, req.BC, req.gamma, RN=req.RN)

        gs = np.linspace(0.5, 30, 80)
        Qs = [R.ballistic_heat_load(req.Ve, req.BC, g, RN=req.RN) for g in gs]
        qms = [R.ballistic_max_heat_rate(req.Ve, req.BC, g, RN=req.RN) for g in gs]

        return {
            "Q": Q,
            "qmax": qmax,
            "gs": _ndarray_to_list(gs),
            "Qs": [float(x) for x in Qs],
            "qms": [float(x) for x in qms],
            "gamma": req.gamma,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 4: Lifting Reentry
# ============================================================
class LiftingRequest(BaseModel):
    planet: str = "Earth"
    Ve: float = 7900.0
    BC: float = 150.0
    LD: float = 0.3
    gamma0: float = 1.5
    RN: float = 1.0
    gamma_angles: Optional[List[float]] = None

@app.post("/api/lifting")
def calc_lifting(req: LiftingRequest):
    try:
        if req.gamma_angles and len(req.gamma_angles) >= 3:
            gamma_angles = req.gamma_angles[:5]
        else:
            gamma_angles = [req.gamma0, req.gamma0 + 2.0, req.gamma0 + 4.0, req.gamma0 - 2.0, req.gamma0 - 4.0]

        Q = R.lifting_heat_load(req.Ve, req.LD, req.BC, req.gamma0, RN=req.RN)
        qmax = R.lifting_max_heat_rate(req.Ve, req.BC, req.gamma0, req.LD, RN=req.RN)

        cfg = R.PLANET_CONFIGS.get(req.planet, R.PLANET_CONFIGS['Earth'])
        h0_km = cfg['h_entry'] / 1000.0

        trajectories = []
        for g_val in gamma_angles:
            traj = R.trajectory_integrate(V0=req.Ve, gamma0_deg=-abs(g_val), h0_km=h0_km,
                                          RN=req.RN, BC=req.BC, L_D=req.LD, planet=req.planet)
            q_peak = float(np.max(traj["q_dot_conv"]))
            trajectories.append({
                "gamma": g_val,
                "q_peak": q_peak,
                "duration": float(traj["t"][-1]),
                "t": _ndarray_to_list(traj["t"]),
                "q_dot_conv": _ndarray_to_list(traj["q_dot_conv"]),
            })

        return {
            "Q": Q,
            "qmax": qmax,
            "trajectories": trajectories,
            "gamma_angles": gamma_angles,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 5: TPS Sizing
# ============================================================
class TPSRequest(BaseModel):
    planet: str = "Earth"
    eps: float = 0.85
    d_mm: float = 50.0
    rho_tile: float = 250.0
    k_tile: float = 0.15
    cp_tile: float = 1200.0
    Cv: float = 900.0
    hv: float = 2.5e6
    glass_d_um: float = 300.0
    glass_k: float = 1.4
    sub_d_mm: float = 5.0
    sub_rho: float = 2810.0
    sub_cp: float = 960.0
    sub_mat_name: str = "Al-7xxx (Default)"
    target_T: float = 450.0
    use_tdep: bool = True
    heat_source: str = "trajectory"
    custom_q_peak: float = 50.0
    custom_duration: float = 600.0

@app.post("/api/tps")
def calc_tps(req: TPSRequest):
    try:
        d_m = max(req.d_mm / 1000.0, 0.005)
        d_cm = req.d_mm / 10.0
        glass_d_m = req.glass_d_um / 1e6
        sub_d_m = req.sub_d_mm / 1000.0

        if "custom" in req.heat_source.lower():
            t_arr = np.linspace(0.0, max(req.custom_duration, 10.0), 500)
            q_pulse = req.custom_q_peak * np.sin(np.pi * t_arr / max(req.custom_duration, 10.0))
            q_pulse = np.maximum(q_pulse, 0.0)
            Q_cumul = np.cumsum(q_pulse) * (t_arr[1] - t_arr[0])
            traj_bal = {"t": t_arr, "q_dot_total": q_pulse, "Q_cumul": Q_cumul}
            traj_glide = traj_bal
        else:
            cfg = R.PLANET_CONFIGS.get(req.planet, R.PLANET_CONFIGS['Earth'])
            h0_km = cfg['h_entry'] / 1000.0
            traj_bal = R.trajectory_integrate(V0=0.9*cfg['V_circ'], gamma0_deg=-15.0, h0_km=h0_km, RN=1.0, BC=400.0, L_D=0.0, planet=req.planet)
            traj_glide = R.trajectory_integrate(V0=0.9*cfg['V_circ'], gamma0_deg=-6.0, h0_km=h0_km, RN=1.0, BC=150.0, L_D=0.3, planet=req.planet)

        q_bal_Wm2 = traj_bal["q_dot_total"] * 1e4
        q_glide_Wm2 = traj_glide["q_dot_total"] * 1e4
        delta_T = max(req.target_T - 300.0, 50.0)

        m_sink = R.heat_sink_mass(float(traj_glide["Q_cumul"][-1] * 1e4), req.Cv, delta_T)
        dm_abl = R.ablative_mass_loss(float(traj_glide["Q_cumul"][-1] * 1e4), req.hv)
        Tw_eq = R.radiative_equilibrium_temp(float(np.max(q_glide_Wm2)), req.eps)

        tps_bal = R.solve_1d_transient_tps_multilayer(
            traj_bal["t"], q_bal_Wm2, d_tile_m=d_m, d_glass_m=glass_d_m, k_glass=req.glass_k,
            rho_tile=req.rho_tile, k_tile=req.k_tile, cp_tile=req.cp_tile,
            d_metal_m=sub_d_m, rho_metal=req.sub_rho, cp_metal=req.sub_cp, use_tdependent=req.use_tdep
        )
        tps_glide = R.solve_1d_transient_tps_multilayer(
            traj_glide["t"], q_glide_Wm2, d_tile_m=d_m, d_glass_m=glass_d_m, k_glass=req.glass_k,
            rho_tile=req.rho_tile, k_tile=req.k_tile, cp_tile=req.cp_tile,
            d_metal_m=sub_d_m, rho_metal=req.sub_rho, cp_metal=req.sub_cp, use_tdependent=req.use_tdep
        )

        T_surf_max_C = float(np.max(tps_glide['T_surface']) - 273.15)
        glass_status = "CRITICAL" if T_surf_max_C >= 1600.0 else "SAFE"

        return {
            "m_sink": m_sink * 0.1,
            "dm_abl": dm_abl * 0.1,
            "Tw_eq": Tw_eq,
            "glass_status": glass_status,
            "T_surf_max_C": T_surf_max_C,
            "T_backwall_max": float(tps_glide['T_backwall_max']),
            "t_soak": float(tps_glide['t_soak_max_s']),
            "t_bal": _ndarray_to_list(traj_bal["t"]),
            "q_bal": _ndarray_to_list(traj_bal["q_dot_total"]),
            "t_glide": _ndarray_to_list(traj_glide["t"]),
            "q_glide": _ndarray_to_list(traj_glide["q_dot_total"]),
            "tps_t": _ndarray_to_list(tps_glide["t"]),
            "T_surface": _ndarray_to_list(tps_glide["T_surface"] - 273.15),
            "T_tile_surface": _ndarray_to_list(tps_glide["T_all_hist"][:, 0] - 273.15),
            "T_center": _ndarray_to_list(tps_glide["T_center"] - 273.15),
            "T_backwall": _ndarray_to_list(tps_glide["T_backwall"] - 273.15),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 5b: TPS Tile Thickness Solver
# ============================================================
class TPSSolveRequest(BaseModel):
    planet: str = "Earth"
    eps: float = 0.85
    rho_tile: float = 250.0
    k_tile: float = 0.15
    cp_tile: float = 1200.0
    glass_d_um: float = 300.0
    glass_k: float = 1.4
    sub_d_mm: float = 5.0
    sub_rho: float = 2810.0
    sub_cp: float = 960.0
    sub_mat_name: str = "Al-7xxx (Default)"
    target_T: float = 450.0
    use_tdep: bool = True
    method: str = "brentq"

@app.post("/api/tps/solve")
def solve_tps(req: TPSSolveRequest):
    try:
        cfg = R.PLANET_CONFIGS.get(req.planet, R.PLANET_CONFIGS['Earth'])
        h0_km = cfg['h_entry'] / 1000.0
        traj_glide = R.trajectory_integrate(V0=0.9*cfg['V_circ'], gamma0_deg=-6.0, h0_km=h0_km,
                                            RN=1.0, BC=150.0, L_D=0.3, planet=req.planet)
        q_glide_Wm2 = traj_glide["q_dot_total"] * 1e4
        sub_d_m = req.sub_d_mm / 1000.0

        res = R.solve_required_tps_thickness(
            traj_glide["t"], q_glide_Wm2, T_backwall_limit=req.target_T,
            rho_tps=req.rho_tile, cp_tps=req.cp_tile, k_tps=req.k_tile, eps=req.eps,
            d_metal=sub_d_m, rho_metal=req.sub_rho, cp_metal=req.sub_cp,
            use_multilayer=True, use_tdependent=req.use_tdep, method=req.method
        )
        sol = res["sol_tps"]
        return {
            "is_attainable": res["is_attainable"],
            "d_req_mm": res["d_req_mm"],
            "T_backwall_max": res["T_backwall_max"],
            "mass_areal_g_cm2": res["mass_areal_g_cm2"],
            "tps_t": _ndarray_to_list(sol["t"]),
            "T_surface": _ndarray_to_list(sol["T_surface"] - 273.15),
            "T_backwall": _ndarray_to_list(sol["T_backwall"] - 273.15),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 6: Reentry Corridor
# ============================================================
class CorridorRequest(BaseModel):
    RN: float = 0.3
    BC: float = 100.0
    LD: float = 1.2
    qtmax: float = 60.0
    nmax: float = 2.5
    qmax: float = 10000.0

@app.post("/api/corridor")
def calc_corridor(req: CorridorRequest):
    try:
        V = np.linspace(200, 7800, 400)
        corr = R.reentry_corridor(V, req.RN, req.BC, req.LD, req.qtmax, req.nmax, req.qmax)

        D_heat = R.drag_from_density(corr["rho_heat"], V, req.BC)
        D_load = R.drag_from_density(corr["rho_load"], V, req.BC)
        D_q = R.drag_from_density(corr["rho_q"], V, req.BC)
        D_u = R.drag_from_density(corr["rho_upper"], V, req.BC)
        D_des = R.corridor_drag_segments(V, BC=req.BC, L_D=req.LD, q_t_max=req.qtmax, RN=req.RN)

        return {
            "V": _ndarray_to_list(V),
            "h_upper": _ndarray_to_list(corr["h_upper"]),
            "h_lower": _ndarray_to_list(corr["h_lower"]),
            "h_heat": _ndarray_to_list(corr["h_heat"]),
            "h_load": _ndarray_to_list(corr["h_load"]),
            "h_q": _ndarray_to_list(corr["h_q"]),
            "corridor_width": _ndarray_to_list(corr["corridor_width"]),
            "D_heat": _ndarray_to_list(D_heat),
            "D_load": _ndarray_to_list(D_load),
            "D_q": _ndarray_to_list(D_q),
            "D_u": _ndarray_to_list(D_u),
            "D_des": _ndarray_to_list(D_des),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 7: Trajectory
# ============================================================
class TrajectoryRequest(BaseModel):
    planet: str = "Earth"
    V0: float = 7800.0
    gamma0: float = -6.0
    h0_km: float = 120.0
    RN: float = 0.3
    BC: float = 150.0
    LD: float = 0.3
    bank: float = 0.0
    dt: float = 0.5
    t_max: float = 0.0

@app.post("/api/trajectory")
def calc_trajectory(req: TrajectoryRequest):
    try:
        traj = R.trajectory_integrate(
            V0=req.V0, gamma0_deg=req.gamma0, h0_km=req.h0_km,
            RN=req.RN, BC=req.BC, L_D=req.LD,
            bank_angle_deg=req.bank, t_max=req.t_max, dt=req.dt,
            planet=req.planet
        )
        iq = int(np.argmax(traj["q_dot_total"]))
        inz = int(np.argmax(traj["n_z"]))
        iqdyn = int(np.argmax(traj["q_dyn"]))
        q_peak_Wm2 = float(traj['q_dot_total'][iq]) * 1e4
        T_surf_max_K = (max(q_peak_Wm2, 0.0) / (0.85 * 5.670374e-8)) ** 0.25

        report = {
            "stop_reason": traj["stop_reason"],
            "flight_time": float(traj["t"][-1]),
            "range_km": float(traj["range_km"][-1]),
            "peak_q": float(traj['q_dot_total'][iq]),
            "peak_q_t": float(traj['t'][iq]),
            "peak_q_h": float(traj['h_km'][iq]),
            "T_surf_max_K": T_surf_max_K,
            "T_surf_max_C": T_surf_max_K - 273.15,
            "peak_nz": float(traj['n_z'][inz]),
            "peak_qdyn": float(traj['q_dyn'][iqdyn] / 1000),
            "max_mach": float(traj['Mach'][0]),
            "total_Q": float(traj['Q_cumul'][-1]),
        }

        return {
            "report": report,
            **_traj_to_json(traj),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 7b: Monte Carlo
# ============================================================
class MonteCarloRequest(BaseModel):
    V0: float = 7800.0
    gamma0: float = -6.0
    h0_km: float = 120.0
    RN: float = 0.3
    BC: float = 150.0
    LD: float = 0.3
    n_runs: int = 50
    gamma_std: float = 0.3
    BC_std: float = 15.0
    LD_std: float = 0.02
    distribution: str = "normal"

@app.post("/api/trajectory/monte-carlo")
def calc_monte_carlo(req: MonteCarloRequest):
    try:
        mc = R.monte_carlo_trajectory(
            req.V0, req.gamma0, req.h0_km, req.RN, req.BC, req.LD,
            n_runs=min(req.n_runs, 100), gamma_std=req.gamma_std,
            BC_std=req.BC_std, LD_std=req.LD_std, distribution=req.distribution
        )
        return mc
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 7c: Target Duration Solver
# ============================================================
class DurationSolverRequest(BaseModel):
    target_t: float = 1200.0
    V0: float = 7800.0
    gamma0: float = -6.0
    h0_km: float = 120.0
    RN: float = 0.3
    BC: float = 150.0
    LD: float = 0.3
    bank: float = 0.0
    vary: str = "gamma"
    method: str = "bounded"

@app.post("/api/trajectory/solve-duration")
def solve_duration(req: DurationSolverRequest):
    try:
        res = R.solve_target_duration(
            req.target_t, req.V0, req.gamma0, req.h0_km, req.RN, req.BC, req.LD,
            bank_angle_deg=req.bank, vary=req.vary, method=req.method
        )
        traj = res["trajectory"]
        return {
            "is_attainable": res["is_attainable"],
            "matched_value": res["matched_value"],
            "achieved_t_s": res["achieved_t_s"],
            "error_s": res["error_s"],
            "achievable_min_s": res["achievable_min_s"],
            "achievable_max_s": res["achievable_max_s"],
            "final_params": res["final_params"],
            "trajectory": _traj_to_json(traj),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 7d: Target Heating Solver
# ============================================================
class HeatingSolverRequest(BaseModel):
    target_val: float = 100.0
    target_type: str = "q_max"
    V0: float = 7800.0
    gamma0: float = -6.0
    h0_km: float = 120.0
    RN: float = 0.3
    BC: float = 150.0
    LD: float = 0.3
    bank: float = 0.0
    vary: str = "gamma"
    method: str = "bounded"

@app.post("/api/trajectory/solve-heating")
def solve_heating(req: HeatingSolverRequest):
    try:
        res = R.solve_target_heating(
            req.target_val, req.target_type, req.V0, req.gamma0, req.h0_km,
            req.RN, req.BC, req.LD, bank_angle_deg=req.bank, vary=req.vary, method=req.method
        )
        traj = res["trajectory"]
        return {
            "is_attainable": res["is_attainable"],
            "matched_value": res["matched_value"],
            "achieved_val": res["achieved_val"],
            "error_val": res["error_val"],
            "achievable_min": res["achievable_min"],
            "achievable_max": res["achievable_max"],
            "final_params": res["final_params"],
            "trajectory": _traj_to_json(traj),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 7e: Target Surface Temperature Solver
# ============================================================
class SurfTempSolverRequest(BaseModel):
    V0: float = 7800.0
    h0_km: float = 120.0
    RN: float = 0.3
    BC: float = 150.0
    LD: float = 0.3
    T_surf_target_K: float = 1533.15
    eps: float = 0.85
    bank: float = 0.0
    planet: str = "Earth"
    method: str = "bisect"

@app.post("/api/trajectory/solve-temp")
def solve_surf_temp(req: SurfTempSolverRequest):
    try:
        res = R.solve_target_surface_temp(
            V0=req.V0, h0_km=req.h0_km, RN=req.RN, BC=req.BC, L_D=req.LD,
            T_surf_target_K=req.T_surf_target_K, eps=req.eps, bank_angle_deg=req.bank,
            planet=req.planet, method=req.method
        )
        return {
            "is_attainable": res["is_attainable"],
            "vary_param": res["vary_param"],
            "matched_value": res["matched_value"],
            "gamma_opt_deg": res["gamma_opt_deg"],
            "RN_opt_m": res["RN_opt_m"],
            "BC_opt_kgm2": res["BC_opt_kgm2"],
            "T_surf_peak_K": res["T_surf_peak_K"],
            "T_surf_peak_C": res["T_surf_peak_C"],
            "T_surf_target_K": res["T_surf_target_K"],
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 8: Parametric Sweep
# ============================================================
class SweepRequest(BaseModel):
    Vatm: float = 7800.0
    gamma: float = -6.0
    BC: float = 150.0
    RN: float = 0.3
    LD: float = 0.3
    eps: float = 0.85
    k_sg: float = 1.7415e-4
    theta: float = 30.0
    sweep_var: str = "gamma"
    sweep_start: float = -5.0
    sweep_end: float = -25.0
    sweep_steps: int = 15

@app.post("/api/parametric-sweep")
def calc_sweep(req: SweepRequest):
    try:
        sweep_vals = np.linspace(req.sweep_start, req.sweep_end, req.sweep_steps)
        q_maxs = []
        Q_totals = []
        flight_times = []
        h_peaks = []

        for val in sweep_vals:
            s_Vatm = val if req.sweep_var == "Vatm" else req.Vatm
            s_gamma = val if req.sweep_var == "gamma" else req.gamma
            s_BC = val if req.sweep_var == "BC" else req.BC

            traj = R.trajectory_integrate(V0=s_Vatm, gamma0_deg=s_gamma, h0_km=120.0,
                                          RN=req.RN, BC=s_BC, L_D=req.LD, bank_angle_deg=0.0, t_max=7200.0, dt=1.0)
            q_maxs.append(float(np.max(traj["q_dot_total"])))
            Q_totals.append(float(traj["Q_cumul"][-1]))
            flight_times.append(float(traj["t"][-1]))
            h_peaks.append(float(traj["h_km"][np.argmax(traj["q_dot_total"])]))

        # Single point calculation
        traj_sp = R.trajectory_integrate(V0=req.Vatm, gamma0_deg=req.gamma, h0_km=120.0,
                                         RN=req.RN, BC=req.BC, L_D=req.LD, t_max=7200.0, dt=0.5)
        q_max_num = float(np.max(traj_sp["q_dot_total"]))
        Q_num = float(traj_sp["Q_cumul"][-1])

        is_lifting = abs(req.LD) > 1e-4
        if not is_lifting:
            q_max_an = R.ballistic_max_heat_rate_analytical(req.Vatm, req.BC, req.gamma, req.RN, k=req.k_sg)
            Q_an = R.ballistic_heat_load_analytical(req.Vatm, req.BC, req.gamma, req.RN, k=req.k_sg)
        else:
            q_max_an = R.lifting_max_heat_rate_analytical(req.Vatm, req.BC, req.LD, req.RN)
            Q_an = R.lifting_heat_load_analytical(req.Vatm, req.BC, req.LD, req.RN)

        q_local = R.local_heating_sphere(q_max_num, req.theta)
        Q_dot_hemi = R.total_heating_rate_sphere(q_max_num, req.RN)

        return {
            "sweep_vals": _ndarray_to_list(sweep_vals),
            "q_maxs": q_maxs,
            "Q_totals": Q_totals,
            "flight_times": flight_times,
            "h_peaks": h_peaks,
            "single_point": {
                "q_max_num": q_max_num,
                "Q_num": Q_num,
                "q_max_an": float(q_max_an),
                "Q_an": float(Q_an),
                "q_local": float(q_local),
                "Q_hemi": float(Q_dot_hemi),
            },
            "traj": _traj_to_json(traj_sp),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 9: OpenFOAM / PATO Export
# ============================================================
class OpenFOAMRequest(BaseModel):
    V0: float = 7800.0
    gamma0: float = -6.0
    h0_km: float = 120.0
    RN: float = 0.3
    BC: float = 150.0
    LD: float = 0.3
    bank: float = 0.0
    dt: float = 0.5
    t_max: float = 0.0
    of_RN: float = 0.3
    eps: float = 0.85
    title: str = "Earth Reentry Trajectory (PATO)"
    dt_out: float = 10.0

@app.post("/api/openfoam")
def calc_openfoam(req: OpenFOAMRequest):
    try:
        traj = R.trajectory_integrate(
            V0=req.V0, gamma0_deg=req.gamma0, h0_km=req.h0_km,
            RN=req.RN, BC=req.BC, L_D=req.LD,
            bank_angle_deg=req.bank, t_max=req.t_max, dt=req.dt,
        )
        res = R.generate_openfoam_pato_boundary_conditions(
            traj, RN=req.of_RN, eps=req.eps, title=req.title, dt_output=req.dt_out
        )
        return {"full_text": res["full_text"], "n_rows": res["n_rows"]}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================
# Tab 10: Reference (static)
# ============================================================
@app.get("/api/reference")
def get_reference():
    presets = {}
    for name, p in R.VEHICLE_PRESETS.items():
        presets[name] = {k: float(v) if isinstance(v, (float, int, np.floating)) else v for k, v in p.items()}
    return {
        "constants": {
            "rho0": R.RHO0_SEALEVEL,
            "V0_circ": R.V0_CIRC,
            "g0": R.G0,
            "sigma_SB": R.SIGMA_SB,
            "R_earth": R.R_EARTH,
            "H_scale": R.H_SCALE,
            "R_gas": R.R_GAS,
        },
        "vehicle_presets": presets,
        "planet_configs": {
            name: {k: float(v) if isinstance(v, (float, int, np.floating, np.integer)) else v
                   for k, v in cfg.items() if k != 'composition'}
            for name, cfg in R.PLANET_CONFIGS.items()
        },
    }

# ============================================================
# Vehicle presets
# ============================================================
@app.get("/api/presets")
def get_presets():
    return {
        "space_shuttle": {"planet": "Earth", "V0": 7850.0, "gamma0": -1.25, "h0_km": 120.0, "RN": 1.0, "BC": 250.0, "LD": 1.0, "bank": 0.0},
        "apollo_11": {"planet": "Earth", "V0": 11030.0, "gamma0": -6.48, "h0_km": 122.0, "RN": 4.69, "BC": 370.0, "LD": 0.31, "bank": 180.0},
        "orion": {"planet": "Earth", "V0": 11000.0, "gamma0": -6.0, "h0_km": 120.0, "RN": 5.0, "BC": 380.0, "LD": 0.30, "bank": 180.0},
        "rlv_td": {"planet": "Earth", "V0": 1750.0, "gamma0": -2.5, "h0_km": 65.0, "RN": 0.5, "BC": 180.0, "LD": 1.2, "bank": 0.0},
        "mars_msl": {"planet": "Mars", "V0": 5800.0, "gamma0": -12.0, "h0_km": 125.0, "RN": 2.25, "BC": 135.0, "LD": 0.24, "bank": 0.0},
        "venus_pioneer": {"planet": "Venus", "V0": 11600.0, "gamma0": -30.0, "h0_km": 150.0, "RN": 0.4, "BC": 290.0, "LD": 0.0, "bank": 0.0},
    }

# ============================================================
# Serve static files
# ============================================================
static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "static")

@app.get("/")
def serve_index():
    return FileResponse(os.path.join(static_dir, "index.html"))

app.mount("/static", StaticFiles(directory=static_dir), name="static")
