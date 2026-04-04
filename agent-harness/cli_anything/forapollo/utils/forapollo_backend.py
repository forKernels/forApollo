"""Backend wrapper for forApollo Python bindings.

Wraps real forapollo.navigation, forapollo.guidance, forapollo.orbital
functions with numpy-friendly interfaces.
"""

import json
import ctypes
from ctypes import c_int, c_double, POINTER
import numpy as np
from typing import Any


def _arr(data, dtype=np.float64):
    return np.ascontiguousarray(data, dtype=dtype)


def _f64p(arr):
    return arr.ctypes.data_as(POINTER(c_double))


def _make_info():
    info = c_int(0)
    return info, ctypes.byref(info)


# =========================================================================
# Navigation group -- Kalman filters
# =========================================================================

def kf_predict(n, x, P, F, Q):
    """Standard Kalman filter predict step.
    x: state (n,), P: covariance (n*n,), F: transition (n*n,), Q: process noise (n*n,).
    """
    from forapollo.navigation import kf_predict as _fn
    x = _arr(x).copy()
    P = _arr(P).copy()
    info = _fn(n, _f64p(x), _f64p(P), _f64p(_arr(F)), _f64p(_arr(Q)))
    return x, P, info


def kf_update(n, m, x, P, z, H, R):
    """Standard Kalman filter update step.
    z: measurement (m,), H: observation (m*n,), R: measurement noise (m*m,).
    """
    from forapollo.navigation import kf_update as _fn
    x = _arr(x).copy()
    P = _arr(P).copy()
    validity = np.ones(m, dtype=np.int32)
    info = _fn(n, m, _f64p(x), _f64p(P), _f64p(_arr(z)), _f64p(_arr(H)),
               _f64p(_arr(R)), validity.ctypes.data_as(POINTER(c_int)))
    return x, P, info


def ekf_predict(n, x, P, Q, dt, model_id=0):
    """Extended Kalman Filter predict with built-in dynamics model."""
    from forapollo.navigation import ekf_predict as _fn
    x = _arr(x).copy()
    P = _arr(P).copy()
    params = np.zeros(1, dtype=np.float64)
    info = _fn(n, _f64p(x), _f64p(P), _f64p(_arr(Q)), float(dt),
               model_id, _f64p(params), 0)
    return x, P, info


def ukf_predict(n, x, P, Q, dt, model_id=0, alpha=1e-3, beta=2.0, kappa=0.0):
    """Unscented Kalman Filter predict step."""
    from forapollo.navigation import ukf_predict as _fn
    x = _arr(x).copy()
    P = _arr(P).copy()
    params = np.zeros(1, dtype=np.float64)
    info = _fn(n, _f64p(x), _f64p(P), model_id, _f64p(params), 0,
               _f64p(_arr(Q)), float(dt), alpha, beta, kappa)
    return x, P, info


def eskf_predict(n, x_nom, dx, P, Q, dt, model_id=0):
    """Error-State Kalman Filter predict (Apollo heritage)."""
    from forapollo.navigation import eskf_predict as _fn
    x_nom = _arr(x_nom).copy()
    dx = _arr(dx).copy()
    P = _arr(P).copy()
    params = np.zeros(1, dtype=np.float64)
    info = _fn(n, _f64p(x_nom), _f64p(dx), _f64p(P),
               model_id, _f64p(params), 0, _f64p(_arr(Q)), float(dt))
    return x_nom, dx, P, info


def eskf_inject(n, x_nom, dx, P):
    """ESKF error-state injection."""
    from forapollo.navigation import eskf_inject as _fn
    x_nom = _arr(x_nom).copy()
    dx = _arr(dx).copy()
    P = _arr(P).copy()
    info = _fn(n, _f64p(x_nom), _f64p(dx), _f64p(P))
    return x_nom, P, info


def if_predict(n, eta, Y, F, Q):
    """Information filter predict step."""
    from forapollo.navigation import if_predict as _fn
    eta = _arr(eta).copy()
    Y = _arr(Y).copy()
    info = _fn(n, _f64p(eta), _f64p(Y), _f64p(_arr(F)), _f64p(_arr(Q)))
    return eta, Y, info


def particle_filter_predict(n, particles, weights, n_particles, Q, dt, model_id=0):
    """SIR Particle filter predict step."""
    particles = _arr(particles).copy()
    weights = _arr(weights).copy()
    # Delegates to forapollo_sir_predict
    return particles, weights


# =========================================================================
# Guidance group
# =========================================================================

def guidance_zem(n, x, t_go, g):
    """Zero-Effort Miss guidance (Apollo P63/P64).
    x: state (2*n: pos+vel), t_go: time-to-go, g: gravity vector.
    """
    from forapollo.guidance import guidance_zem as _fn
    x = _arr(x)
    g = _arr(g)
    a_cmd = np.zeros(n, dtype=np.float64)
    info = _fn(n, _f64p(x), float(t_go), _f64p(g), _f64p(a_cmd))
    return a_cmd, info


def guidance_zev(n, x, v_target, t_go, g):
    """Zero-Effort Velocity guidance."""
    from forapollo.guidance import guidance_zev as _fn
    x = _arr(x)
    v_target = _arr(v_target)
    g = _arr(g)
    a_cmd = np.zeros(n, dtype=np.float64)
    info = _fn(n, _f64p(x), _f64p(v_target), float(t_go), _f64p(g), _f64p(a_cmd))
    return a_cmd, info


def guidance_lambert(r1, r2, tof, mu, direction=0):
    """Lambert problem solver (compute transfer orbit)."""
    from forapollo.guidance import guidance_lambert as _fn
    r1 = _arr(r1)
    r2 = _arr(r2)
    v1 = np.zeros(3, dtype=np.float64)
    v2 = np.zeros(3, dtype=np.float64)
    info_val, info_p = _make_info()
    from forapollo import lib as _lib
    _lib.forapollo_guidance_lambert(_f64p(r1), _f64p(r2), c_double(tof), c_double(mu),
                                    _f64p(v1), _f64p(v2), c_int(direction), info_p)
    return v1, v2, info_val.value


def guidance_lqr(n, m, A, B, Q, R):
    """Linear-Quadratic Regulator. Returns gain matrix K and cost S."""
    from forapollo.guidance import guidance_lqr as _fn
    A = _arr(A)
    B = _arr(B)
    Q_mat = _arr(Q)
    R_mat = _arr(R)
    K = np.zeros(m * n, dtype=np.float64)
    S = np.zeros(n * n, dtype=np.float64)
    info_val, info_p = _make_info()
    from forapollo import lib as _lib
    _lib.forapollo_guidance_lqr(c_int(n), c_int(m), _f64p(A), _f64p(B),
                                 _f64p(Q_mat), _f64p(R_mat), _f64p(K), _f64p(S), info_p)
    return K.reshape(m, n), S.reshape(n, n), info_val.value


def guidance_pure_pursuit(position, target, lookahead):
    """Pure pursuit path tracking. Returns steering angle."""
    from forapollo.guidance import guidance_pure_pursuit as _fn
    pos = _arr(position)
    tgt = _arr(target)
    steer = c_double(0.0)
    info_val, info_p = _make_info()
    from forapollo import lib as _lib
    _lib.forapollo_guidance_pure_pursuit(_f64p(pos), _f64p(tgt), c_double(lookahead),
                                         ctypes.byref(steer), info_p)
    return float(steer.value), info_val.value


def guidance_stanley(position, path_point, heading, speed, k):
    """Stanley controller. Returns steering angle."""
    from forapollo.guidance import guidance_stanley as _fn
    pos = _arr(position)
    path = _arr(path_point)
    steer = c_double(0.0)
    info_val, info_p = _make_info()
    from forapollo import lib as _lib
    _lib.forapollo_guidance_stanley(_f64p(pos), _f64p(path), c_double(heading),
                                     c_double(speed), c_double(k),
                                     ctypes.byref(steer), info_p)
    return float(steer.value), info_val.value


def guidance_pn_pure(n, x_pursuer, x_target, N_gain):
    """Pure proportional navigation guidance."""
    from forapollo.guidance import guidance_pn_pure as _fn
    xp = _arr(x_pursuer)
    xt = _arr(x_target)
    a_cmd = np.zeros(n, dtype=np.float64)
    info = _fn(n, _f64p(xp), _f64p(xt), float(N_gain), _f64p(a_cmd))
    return a_cmd, info


def guidance_min_energy(n, x, x_target, t_go):
    """Minimum energy guidance."""
    from forapollo.guidance import guidance_min_energy as _fn
    x = _arr(x)
    xt = _arr(x_target)
    a_cmd = np.zeros(n, dtype=np.float64)
    info_val, info_p = _make_info()
    from forapollo import lib as _lib
    _lib.forapollo_guidance_min_energy(c_int(n), _f64p(x), _f64p(xt),
                                        c_double(t_go), _f64p(a_cmd), info_p)
    return a_cmd, info_val.value


# =========================================================================
# Orbital group -- coordinate transforms
# =========================================================================

def coords_eci_to_ecef(x_eci, gmst):
    """Convert ECI (J2000) to ECEF coordinates."""
    from forapollo.orbital import coords_eci_to_ecef as _fn
    x_eci = _arr(x_eci)
    x_ecef = np.zeros(3, dtype=np.float64)
    info = _fn(_f64p(x_eci), _f64p(x_ecef), float(gmst))
    return x_ecef, info


def coords_ecef_to_eci(x_ecef, gmst):
    """Convert ECEF to ECI (J2000) coordinates."""
    from forapollo.orbital import coords_ecef_to_eci as _fn
    x_ecef = _arr(x_ecef)
    x_eci = np.zeros(3, dtype=np.float64)
    info = _fn(_f64p(x_ecef), _f64p(x_eci), float(gmst))
    return x_eci, info


def coords_ecef_to_geodetic(x_ecef, a_body=6378137.0, f_body=1.0/298.257223563):
    """Convert ECEF to geodetic (lat, lon, alt) using WGS84."""
    from forapollo.orbital import coords_ecef_to_geodetic as _fn
    x_ecef = _arr(x_ecef)
    info, lat, lon, alt = _fn(_f64p(x_ecef), a_body, f_body)
    return {"lat": lat, "lon": lon, "alt": alt}, info


def coords_geodetic_to_ecef(lat, lon, alt, a_body=6378137.0, f_body=1.0/298.257223563):
    """Convert geodetic (lat, lon, alt) to ECEF."""
    from forapollo.orbital import coords_geodetic_to_ecef as _fn
    x_ecef = np.zeros(3, dtype=np.float64)
    info = _fn(float(lat), float(lon), float(alt), _f64p(x_ecef), a_body, f_body)
    return x_ecef, info
