// Copyright The Fantastic Planet — By David Clabaugh
//
// dispatch.zig — size-based routing between Zig-native and Fortran implementations
//
// Small state (n <= threshold) can route to Zig kernels (zero FFI overhead).
// Large state (n > threshold) routes to Fortran via extern fn.
//
// Currently ALL dispatch goes to Fortran. The Zig-native kernels in
// kernels.zig are a future optimization. This dispatch layer exists so
// the routing can be changed without modifying the export API.

const fortran = @import("fortran.zig");
const kernels = @import("kernels.zig");

// ============================================================================
// Platform-configurable thresholds
// ============================================================================

pub const Thresholds = struct {
    /// KF: n <= kf_max_zig uses Zig-native, n > uses Fortran
    kf_max_zig: i32 = 6,
    /// Dynamics: always Fortran (model catalog lives in Fortran)
    dynamics_max_zig: i32 = 0,
    /// Observations: always Fortran (model catalog lives in Fortran)
    observe_max_zig: i32 = 0,
};

pub const default_thresholds = Thresholds{};

// ============================================================================
// Dynamics dispatch
// ============================================================================

pub fn dynamicsDispatch(
    model_id: i32,
    n: i32,
    x: [*]const f64,
    u: [*]const f64,
    nu: i32,
    t: f64,
    params: [*]const f64,
    np: i32,
    x_dot: [*]f64,
    info: *i32,
) void {
    // Always Fortran — dynamics model catalog is in Fortran
    fortran.fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info);
}

pub fn dynamicsJacobian(
    model_id: i32,
    n: i32,
    x: [*]const f64,
    u: [*]const f64,
    nu: i32,
    t: f64,
    params: [*]const f64,
    np: i32,
    F: [*]f64,
    info: *i32,
) void {
    fortran.fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info);
}

// ============================================================================
// Observation dispatch
// ============================================================================

pub fn observeDispatch(
    obs_id: i32,
    n: i32,
    x: [*]const f64,
    m: i32,
    t: f64,
    obs_params: [*]const f64,
    nop: i32,
    z_pred: [*]f64,
    info: *i32,
) void {
    fortran.fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info);
}

pub fn observeJacobian(
    obs_id: i32,
    n: i32,
    x: [*]const f64,
    m: i32,
    t: f64,
    obs_params: [*]const f64,
    nop: i32,
    H: [*]f64,
    info: *i32,
) void {
    fortran.fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, info);
}

// ============================================================================
// Propagation dispatch
// ============================================================================

pub fn propagate(
    n: i32,
    x: [*]f64,
    u: [*]const f64,
    nu: i32,
    f_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    dt: f64,
    n_steps: i32,
    info: *i32,
) void {
    fortran.fa_propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, n_steps, info);
}

pub fn propagateStm(
    n: i32,
    x: [*]f64,
    phi: [*]f64,
    u: [*]const f64,
    nu: i32,
    f_ptr: ?*const anyopaque,
    df_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    dt: f64,
    n_steps: i32,
    info: *i32,
) void {
    fortran.fa_propagate_stm(n, x, phi, u, nu, f_ptr, df_ptr, model_id, params, np, dt, n_steps, info);
}

pub fn propagateBatch(
    n: i32,
    n_states: i32,
    x_batch: [*]f64,
    u: [*]const f64,
    nu: i32,
    f_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    dt: f64,
    n_steps: i32,
    info_batch: [*]i32,
) void {
    fortran.fa_propagate_batch(n, n_states, x_batch, u, nu, f_ptr, model_id, params, np, dt, n_steps, info_batch);
}

// ============================================================================
// KF dispatch
// ============================================================================

pub fn kfPredict(
    n: i32,
    x: [*]f64,
    P: [*]f64,
    F: [*]const f64,
    Q: [*]const f64,
    info: *i32,
) void {
    // TODO: n <= kf_max_zig → kernels.kfPredict (Zig-native)
    fortran.fa_kf_predict(n, x, P, F, Q, info);
}

pub fn kfUpdate(
    n: i32,
    m: i32,
    x: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    H: [*]const f64,
    R: [*]const f64,
    validity: [*]i32,
    info: *i32,
) void {
    // TODO: n <= kf_max_zig → kernels.kfUpdate (Zig-native)
    fortran.fa_kf_update(n, m, x, P, z, H, R, validity, info);
}

// ============================================================================
// EKF dispatch
// ============================================================================

pub fn ekfPredict(
    n: i32,
    x: [*]f64,
    P: [*]f64,
    f_ptr: ?*const anyopaque,
    df_ptr: ?*const anyopaque,
    Q: [*]const f64,
    dt: f64,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    n_steps: i32,
    info: *i32,
) void {
    fortran.fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, n_steps, info);
}

pub fn ekfUpdate(
    n: i32,
    m: i32,
    x: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    dh_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_ekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info);
}

// ============================================================================
// IEKF dispatch
// ============================================================================

pub fn iekfUpdate(
    n: i32,
    m: i32,
    x: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    dh_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    max_iter: i32,
    tol: f64,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_iekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, max_iter, tol, validity, info);
}

// ============================================================================
// UKF dispatch
// ============================================================================

pub fn ukfPredict(
    n: i32,
    x: [*]f64,
    P: [*]f64,
    f_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    Q: [*]const f64,
    dt: f64,
    alpha: f64,
    beta_ukf: f64,
    kappa: f64,
    info: *i32,
) void {
    fortran.fa_ukf_predict(n, x, P, f_ptr, model_id, params, np, Q, dt, alpha, beta_ukf, kappa, info);
}

pub fn ukfUpdate(
    n: i32,
    m: i32,
    x: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    alpha: f64,
    beta_ukf: f64,
    kappa: f64,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_ukf_update(n, m, x, P, z, h_ptr, obs_id, R, obs_params, nop, alpha, beta_ukf, kappa, validity, info);
}

// ============================================================================
// ESKF dispatch
// ============================================================================

pub fn eskfPredict(
    n: i32,
    x_nom: [*]f64,
    dx: [*]f64,
    P: [*]f64,
    f_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    Q: [*]const f64,
    dt: f64,
    n_steps: i32,
    info: *i32,
) void {
    fortran.fa_eskf_predict(n, x_nom, dx, P, f_ptr, model_id, params, np, Q, dt, n_steps, info);
}

pub fn eskfUpdate(
    n: i32,
    m: i32,
    x_nom: [*]const f64,
    dx: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    dh_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_eskf_update(n, m, x_nom, dx, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info);
}

pub fn eskfInject(
    n: i32,
    x_nom: [*]f64,
    dx: [*]f64,
    P: [*]f64,
    info: *i32,
) void {
    fortran.fa_eskf_inject(n, x_nom, dx, P, info);
}

// ============================================================================
// SREKF dispatch
// ============================================================================

pub fn srekfPredict(
    n: i32,
    x: [*]f64,
    S: [*]f64,
    f_ptr: ?*const anyopaque,
    df_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    Sq: [*]const f64,
    dt: f64,
    n_steps: i32,
    info: *i32,
) void {
    fortran.fa_srekf_predict(n, x, S, f_ptr, df_ptr, model_id, params, np, Sq, dt, n_steps, info);
}

pub fn srekfUpdate(
    n: i32,
    m: i32,
    x: [*]f64,
    S: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    dh_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_srekf_update(n, m, x, S, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info);
}

// ============================================================================
// SRUKF dispatch
// ============================================================================

pub fn srukfPredict(
    n: i32,
    x: [*]f64,
    S: [*]f64,
    f_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    Sq: [*]const f64,
    dt: f64,
    alpha: f64,
    beta_ukf: f64,
    kappa: f64,
    info: *i32,
) void {
    fortran.fa_srukf_predict(n, x, S, f_ptr, model_id, params, np, Sq, dt, alpha, beta_ukf, kappa, info);
}

pub fn srukfUpdate(
    n: i32,
    m: i32,
    x: [*]f64,
    S: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    alpha: f64,
    beta_ukf: f64,
    kappa: f64,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_srukf_update(n, m, x, S, z, h_ptr, obs_id, R, obs_params, nop, alpha, beta_ukf, kappa, validity, info);
}

// ============================================================================
// Information Filter dispatch
// ============================================================================

pub fn ifPredict(
    n: i32,
    eta: [*]f64,
    Ymat: [*]f64,
    F: [*]const f64,
    Q: [*]const f64,
    info: *i32,
) void {
    fortran.fa_if_predict(n, eta, Ymat, F, Q, info);
}

pub fn ifUpdate(
    n: i32,
    m: i32,
    eta: [*]f64,
    Ymat: [*]f64,
    z: [*]const f64,
    H: [*]const f64,
    R: [*]const f64,
    validity: [*]i32,
    info: *i32,
) void {
    fortran.fa_if_update(n, m, eta, Ymat, z, H, R, validity, info);
}

// ============================================================================
// Particle Filter dispatch
// ============================================================================

pub fn pfSirPredict(
    n: i32,
    n_particles: i32,
    particles: [*]f64,
    weights: [*]f64,
    f_ptr: ?*const anyopaque,
    model_id: i32,
    params: [*]const f64,
    np: i32,
    Q: [*]const f64,
    dt: f64,
    n_steps: i32,
    seed: *i32,
    info: *i32,
) void {
    fortran.fa_pf_sir_predict(n, n_particles, particles, weights, f_ptr, model_id, params, np, Q, dt, n_steps, seed, info);
}

pub fn pfSirUpdate(
    n: i32,
    m: i32,
    n_particles: i32,
    particles: [*]const f64,
    weights: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    info: *i32,
) void {
    fortran.fa_pf_sir_update(n, m, n_particles, particles, weights, z, h_ptr, obs_id, R, obs_params, nop, info);
}

pub fn pfSirResample(
    n: i32,
    n_particles: i32,
    particles: [*]f64,
    weights: [*]f64,
    seed: *i32,
    info: *i32,
) void {
    fortran.fa_pf_sir_resample(n, n_particles, particles, weights, seed, info);
}

pub fn pfRbUpdate(
    n: i32,
    m: i32,
    n_particles: i32,
    n_linear: i32,
    particles: [*]f64,
    weights: [*]f64,
    z: [*]const f64,
    h_ptr: ?*const anyopaque,
    obs_id: i32,
    R: [*]const f64,
    obs_params: [*]const f64,
    nop: i32,
    info: *i32,
) void {
    fortran.fa_pf_rb_update(n, m, n_particles, n_linear, particles, weights, z, h_ptr, obs_id, R, obs_params, nop, info);
}

// ============================================================================
// Smoother dispatch
// ============================================================================

pub fn rtsSmooth(
    n: i32,
    nsteps: i32,
    x_filt: [*]const f64,
    P_filt: [*]const f64,
    x_pred: [*]const f64,
    P_pred: [*]const f64,
    F_all: [*]const f64,
    x_smooth: [*]f64,
    P_smooth: [*]f64,
    info: *i32,
) void {
    fortran.fa_rts_smooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth, info);
}

pub fn urtssSmooth(
    n: i32,
    nsteps: i32,
    x_filt: [*]const f64,
    P_filt: [*]const f64,
    x_pred: [*]const f64,
    P_pred: [*]const f64,
    F_all: [*]const f64,
    x_smooth: [*]f64,
    P_smooth: [*]f64,
    info: *i32,
) void {
    fortran.fa_urtss_smooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth, info);
}

// ============================================================================
// Batch estimator dispatch
// ============================================================================

pub fn batchWls(
    n: i32,
    m_total: i32,
    x: [*]f64,
    z_all: [*]const f64,
    H_all: [*]const f64,
    R_all: [*]const f64,
    max_iter: i32,
    tol: f64,
    info: *i32,
) void {
    fortran.fa_batch_wls(n, m_total, x, z_all, H_all, R_all, max_iter, tol, info);
}

pub fn batchMap(
    n: i32,
    m_total: i32,
    x: [*]f64,
    x0: [*]const f64,
    P0: [*]const f64,
    z_all: [*]const f64,
    H_all: [*]const f64,
    R_all: [*]const f64,
    max_iter: i32,
    tol: f64,
    info: *i32,
) void {
    fortran.fa_batch_map(n, m_total, x, x0, P0, z_all, H_all, R_all, max_iter, tol, info);
}
