"""forApollo guidance module -- ZEM/ZEV, proportional nav, Lambert, LQR, pure pursuit, Stanley."""

import ctypes
from ctypes import c_int, c_double, POINTER

from forapollo import lib as _lib

_f64p = POINTER(c_double)

# ---------------------------------------------------------------------------
# Guidance laws (forapollo_guidance_*)
# ---------------------------------------------------------------------------

_lib.forapollo_guidance_zem.argtypes = [c_int, _f64p, c_double, _f64p, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_zem.restype = None

_lib.forapollo_guidance_zev.argtypes = [c_int, _f64p, _f64p, c_double, _f64p, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_zev.restype = None

_lib.forapollo_guidance_eguidance.argtypes = [c_int, _f64p, _f64p, c_double, _f64p, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_eguidance.restype = None

_lib.forapollo_guidance_pn_pure.argtypes = [c_int, _f64p, _f64p, c_double, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_pn_pure.restype = None

_lib.forapollo_guidance_pn_aug.argtypes = [c_int, _f64p, _f64p, c_double, _f64p, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_pn_aug.restype = None

_lib.forapollo_guidance_pn_true.argtypes = [c_int, _f64p, _f64p, c_double, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_pn_true.restype = None

_lib.forapollo_guidance_lambert.argtypes = [_f64p, _f64p, c_double, c_double, _f64p, _f64p, c_int, POINTER(c_int)]
_lib.forapollo_guidance_lambert.restype = None

_lib.forapollo_guidance_lqr.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _f64p, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_lqr.restype = None

_lib.forapollo_guidance_pure_pursuit.argtypes = [_f64p, _f64p, c_double, POINTER(c_double), POINTER(c_int)]
_lib.forapollo_guidance_pure_pursuit.restype = None

_lib.forapollo_guidance_stanley.argtypes = [_f64p, _f64p, c_double, c_double, c_double, POINTER(c_double), POINTER(c_int)]
_lib.forapollo_guidance_stanley.restype = None

_lib.forapollo_guidance_traj_track.argtypes = [c_int, _f64p, _f64p, _f64p, c_double, c_double, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_traj_track.restype = None

_lib.forapollo_guidance_min_energy.argtypes = [c_int, _f64p, _f64p, c_double, _f64p, POINTER(c_int)]
_lib.forapollo_guidance_min_energy.restype = None


def _make_info() -> tuple:
    info = c_int(0)
    return info, ctypes.byref(info)


def guidance_zem(n: int, x, t_go: float, g, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_zem(n, x, t_go, g, a_cmd, info_p)
    return info.value

def guidance_zev(n: int, x, v_target, t_go: float, g, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_zev(n, x, v_target, t_go, g, a_cmd, info_p)
    return info.value

def guidance_eguidance(n: int, x, x_target, t_go: float, g, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_eguidance(n, x, x_target, t_go, g, a_cmd, info_p)
    return info.value

def guidance_pn_pure(n: int, x_p, x_t, N_gain: float, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_pn_pure(n, x_p, x_t, N_gain, a_cmd, info_p)
    return info.value

def guidance_pn_aug(n: int, x_p, x_t, N_gain: float, a_target, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_pn_aug(n, x_p, x_t, N_gain, a_target, a_cmd, info_p)
    return info.value

def guidance_pn_true(n: int, x_p, x_t, N_gain: float, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_pn_true(n, x_p, x_t, N_gain, a_cmd, info_p)
    return info.value

def guidance_lambert(r1, r2, tof: float, mu: float, v1, v2, n_rev: int = 0) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_lambert(r1, r2, tof, mu, v1, v2, n_rev, info_p)
    return info.value

def guidance_lqr(n: int, m: int, A, B, Q_cost, R_cost, x, u_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_lqr(n, m, A, B, Q_cost, R_cost, x, u_cmd, info_p)
    return info.value

def guidance_pure_pursuit(x_pos, x_lookahead, L_wheelbase: float) -> tuple[int, float]:
    info, info_p = _make_info()
    steer = c_double()
    _lib.forapollo_guidance_pure_pursuit(x_pos, x_lookahead, L_wheelbase, ctypes.byref(steer), info_p)
    return info.value, steer.value

def guidance_stanley(x_pos, path_point, path_heading: float, v: float, k_gain: float) -> tuple[int, float]:
    info, info_p = _make_info()
    steer = c_double()
    _lib.forapollo_guidance_stanley(x_pos, path_point, path_heading, v, k_gain, ctypes.byref(steer), info_p)
    return info.value, steer.value

def guidance_traj_track(n: int, x, x_ref, v_ref, K_pos: float, K_vel: float, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_traj_track(n, x, x_ref, v_ref, K_pos, K_vel, a_cmd, info_p)
    return info.value

def guidance_min_energy(n: int, x, x_target, t_go: float, a_cmd) -> int:
    info, info_p = _make_info()
    _lib.forapollo_guidance_min_energy(n, x, x_target, t_go, a_cmd, info_p)
    return info.value
