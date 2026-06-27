"""forApollo navigation module -- Kalman filters (KF, EKF, IEKF, UKF, ESKF, SR-EKF, SR-UKF, IF), particle filters, smoothers, batch estimators."""

import ctypes
from ctypes import c_int, c_double, c_void_p, POINTER

from forapollo import lib as _lib

_f64p = POINTER(c_double)
_i32p = POINTER(c_int)
_opaque = c_void_p


def _make_info():
    info = c_int(0)
    return info, ctypes.byref(info)


# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------

_lib.fa_version.argtypes = []
_lib.fa_version.restype = ctypes.c_uint32

def version() -> int:
    return _lib.fa_version()


# ---------------------------------------------------------------------------
# Dynamics / Observation dispatch
# ---------------------------------------------------------------------------

_lib.fa_dynamics_dispatch.argtypes = [c_int, c_int, _f64p, _f64p, c_int, c_double, _f64p, c_int, _f64p, _i32p]
_lib.fa_dynamics_dispatch.restype = None

_lib.fa_dynamics_jacobian.argtypes = [c_int, c_int, _f64p, _f64p, c_int, c_double, _f64p, c_int, _f64p, _i32p]
_lib.fa_dynamics_jacobian.restype = None

_lib.fa_observe_dispatch.argtypes = [c_int, c_int, _f64p, c_int, c_double, _f64p, c_int, _f64p, _i32p]
_lib.fa_observe_dispatch.restype = None

_lib.fa_observe_jacobian.argtypes = [c_int, c_int, _f64p, c_int, c_double, _f64p, c_int, _f64p, _i32p]
_lib.fa_observe_jacobian.restype = None


def dynamics_dispatch(model_id, n, x, u, nu, t, params, np_, x_dot) -> int:
    info, ip = _make_info()
    _lib.fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np_, x_dot, ip)
    return info.value

def dynamics_jacobian(model_id, n, x, u, nu, t, params, np_, F) -> int:
    info, ip = _make_info()
    _lib.fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np_, F, ip)
    return info.value

def observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred) -> int:
    info, ip = _make_info()
    _lib.fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, ip)
    return info.value

def observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H) -> int:
    info, ip = _make_info()
    _lib.fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, ip)
    return info.value


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------

_lib.fa_propagate.argtypes = [c_int, _f64p, _f64p, c_int, _opaque, c_int, _f64p, c_int, c_double, c_int, _i32p]
_lib.fa_propagate.restype = None

_lib.fa_propagate_stm.argtypes = [c_int, _f64p, _f64p, _f64p, c_int, _opaque, _opaque, c_int, _f64p, c_int, c_double, c_int, _i32p]
_lib.fa_propagate_stm.restype = None

_lib.fa_propagate_batch.argtypes = [c_int, c_int, _f64p, _f64p, c_int, _opaque, c_int, _f64p, c_int, c_double, c_int, _i32p]
_lib.fa_propagate_batch.restype = None


def propagate(n, x, u, nu, model_id, params, np_, dt, n_steps) -> int:
    info, ip = _make_info()
    _lib.fa_propagate(n, x, u, nu, None, model_id, params, np_, dt, n_steps, ip)
    return info.value


# ---------------------------------------------------------------------------
# Kalman filter (KF)
# ---------------------------------------------------------------------------

_lib.fa_kf_predict.argtypes = [c_int, _f64p, _f64p, _f64p, _f64p, _i32p]
_lib.fa_kf_predict.restype = None

_lib.fa_kf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _f64p, _i32p, _i32p]
_lib.fa_kf_update.restype = None


def kf_predict(n, x, P, F, Q) -> int:
    info, ip = _make_info()
    _lib.fa_kf_predict(n, x, P, F, Q, ip)
    return info.value

def kf_update(n, m, x, P, z, H, R, validity) -> int:
    info, ip = _make_info()
    _lib.fa_kf_update(n, m, x, P, z, H, R, validity, ip)
    return info.value


# ---------------------------------------------------------------------------
# Extended Kalman Filter (EKF)
# ---------------------------------------------------------------------------

_lib.fa_ekf_predict.argtypes = [c_int, _f64p, _f64p, _opaque, _opaque, _f64p, c_double, c_int, _f64p, c_int, c_int, _i32p]
_lib.fa_ekf_predict.restype = None

_lib.fa_ekf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _opaque, _opaque, c_int, _f64p, _f64p, c_int, _i32p, _i32p]
_lib.fa_ekf_update.restype = None


def ekf_predict(n, x, P, Q, dt, model_id, params, np_, n_steps=1) -> int:
    info, ip = _make_info()
    _lib.fa_ekf_predict(n, x, P, None, None, Q, dt, model_id, params, np_, n_steps, ip)
    return info.value

def ekf_update(n, m, x, P, z, obs_id, R, obs_params, nop, validity) -> int:
    info, ip = _make_info()
    _lib.fa_ekf_update(n, m, x, P, z, None, None, obs_id, R, obs_params, nop, validity, ip)
    return info.value


# ---------------------------------------------------------------------------
# IEKF, UKF, ESKF
# ---------------------------------------------------------------------------

_lib.fa_iekf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _opaque, _opaque, c_int, _f64p, _f64p, c_int, c_int, c_double, _i32p, _i32p]
_lib.fa_iekf_update.restype = None

_lib.fa_ukf_predict.argtypes = [c_int, _f64p, _f64p, _opaque, c_int, _f64p, c_int, _f64p, c_double, c_double, c_double, c_double, _i32p]
_lib.fa_ukf_predict.restype = None

_lib.fa_ukf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _opaque, c_int, _f64p, _f64p, c_int, c_double, c_double, c_double, _i32p, _i32p]
_lib.fa_ukf_update.restype = None

_lib.fa_eskf_predict.argtypes = [c_int, _f64p, _f64p, _f64p, _opaque, c_int, _f64p, c_int, _f64p, c_double, c_int, _i32p]
_lib.fa_eskf_predict.restype = None

_lib.fa_eskf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _opaque, _opaque, c_int, _f64p, _f64p, c_int, _i32p, _i32p]
_lib.fa_eskf_update.restype = None

_lib.fa_eskf_inject.argtypes = [c_int, _f64p, _f64p, _f64p, _i32p]
_lib.fa_eskf_inject.restype = None


def ukf_predict(n, x, P, model_id, params, np_, Q, dt, alpha=1e-3, beta=2.0, kappa=0.0) -> int:
    info, ip = _make_info()
    _lib.fa_ukf_predict(n, x, P, None, model_id, params, np_, Q, dt, alpha, beta, kappa, ip)
    return info.value

def eskf_predict(n, x_nom, dx, P, model_id, params, np_, Q, dt, n_steps=1) -> int:
    info, ip = _make_info()
    _lib.fa_eskf_predict(n, x_nom, dx, P, None, model_id, params, np_, Q, dt, n_steps, ip)
    return info.value

def eskf_inject(n, x_nom, dx, P) -> int:
    info, ip = _make_info()
    _lib.fa_eskf_inject(n, x_nom, dx, P, ip)
    return info.value


# ---------------------------------------------------------------------------
# Square-root filters (SR-EKF, SR-UKF)
# ---------------------------------------------------------------------------

_lib.fa_srekf_predict.argtypes = [c_int, _f64p, _f64p, _opaque, _opaque, c_int, _f64p, c_int, _f64p, c_double, c_int, _i32p]
_lib.fa_srekf_predict.restype = None

_lib.fa_srekf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _opaque, _opaque, c_int, _f64p, _f64p, c_int, _i32p, _i32p]
_lib.fa_srekf_update.restype = None

_lib.fa_srukf_predict.argtypes = [c_int, _f64p, _f64p, _opaque, c_int, _f64p, c_int, _f64p, c_double, c_double, c_double, c_double, _i32p]
_lib.fa_srukf_predict.restype = None

_lib.fa_srukf_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _opaque, c_int, _f64p, _f64p, c_int, c_double, c_double, c_double, _i32p, _i32p]
_lib.fa_srukf_update.restype = None


# ---------------------------------------------------------------------------
# Information filter
# ---------------------------------------------------------------------------

_lib.fa_if_predict.argtypes = [c_int, _f64p, _f64p, _f64p, _f64p, _i32p]
_lib.fa_if_predict.restype = None

_lib.fa_if_update.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _f64p, _i32p, _i32p]
_lib.fa_if_update.restype = None


def if_predict(n, eta, Y, F, Q) -> int:
    info, ip = _make_info()
    _lib.fa_if_predict(n, eta, Y, F, Q, ip)
    return info.value

def if_update(n, m, eta, Y, z, H, R, validity) -> int:
    info, ip = _make_info()
    _lib.fa_if_update(n, m, eta, Y, z, H, R, validity, ip)
    return info.value


# ---------------------------------------------------------------------------
# Particle filter (SIR)
# ---------------------------------------------------------------------------

_lib.fa_pf_sir_predict.argtypes = [c_int, c_int, _f64p, _f64p, _opaque, c_int, _f64p, c_int, _f64p, c_double, c_int, _i32p, _i32p]
_lib.fa_pf_sir_predict.restype = None

_lib.fa_pf_sir_update.argtypes = [c_int, c_int, c_int, _f64p, _f64p, _f64p, _opaque, c_int, _f64p, _f64p, c_int, _i32p]
_lib.fa_pf_sir_update.restype = None

_lib.fa_pf_sir_resample.argtypes = [c_int, c_int, _f64p, _f64p, _i32p, _i32p]
_lib.fa_pf_sir_resample.restype = None

_lib.fa_pf_rb_update.argtypes = [c_int, c_int, c_int, c_int, _f64p, _f64p, _f64p, _opaque, c_int, _f64p, _f64p, c_int, _i32p]
_lib.fa_pf_rb_update.restype = None


# ---------------------------------------------------------------------------
# Smoothers (RTS, URTSS)
# ---------------------------------------------------------------------------

_lib.fa_rts_smooth.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _f64p, _f64p, _f64p, _i32p]
_lib.fa_rts_smooth.restype = None

_lib.fa_urtss_smooth.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _f64p, _f64p, _f64p, _i32p]
_lib.fa_urtss_smooth.restype = None


def rts_smooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth) -> int:
    info, ip = _make_info()
    _lib.fa_rts_smooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth, ip)
    return info.value


# ---------------------------------------------------------------------------
# Batch estimators (WLS, MAP)
# ---------------------------------------------------------------------------

_lib.fa_batch_wls.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, c_int, c_double, _i32p]
_lib.fa_batch_wls.restype = None

_lib.fa_batch_map.argtypes = [c_int, c_int, _f64p, _f64p, _f64p, _f64p, _f64p, _f64p, c_int, c_double, _i32p]
_lib.fa_batch_map.restype = None


def batch_wls(n, m_total, x, z_all, H_all, R_all, max_iter=10, tol=1e-8) -> int:
    info, ip = _make_info()
    _lib.fa_batch_wls(n, m_total, x, z_all, H_all, R_all, max_iter, tol, ip)
    return info.value
