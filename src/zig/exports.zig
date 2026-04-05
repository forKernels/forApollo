// Copyright The Fantastic Planet — By David Clabaugh
//
// exports.zig — C ABI exports for Python/Rust/C++ consumers
//
// All functions prefixed forapollo_* and call fa_* Fortran kernels through
// the safety + dispatch layers. This is the public API of libforapollo.
//
// Flow: caller → forapollo_* (validate) → dispatch → fortran fa_*

const safety = @import("safety.zig");
const fortran = @import("fortran.zig");
const dispatch = @import("dispatch.zig");
pub const cuda = @import("cuda.zig");

// Force the compiler to pull in all referenced modules so link errors
// surface at build time rather than at dlopen time.
comptime {
    _ = @import("kernels.zig");
}

// ============================================================================
// Version
// ============================================================================

/// Returns the library version as a packed u32: 0xMMmmpp (major.minor.patch).
pub export fn forapollo_version() callconv(.c) u32 {
    return 0x000100; // 0.1.0
}

// ============================================================================
// Dynamics
// ============================================================================

pub export fn forapollo_dynamics_dispatch(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.dynamicsDispatch(model_id, n, x, u, nu, t, params, np, x_dot, info);
}

pub export fn forapollo_dynamics_jacobian(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.dynamicsJacobian(model_id, n, x, u, nu, t, params, np, F, info);
}

// ============================================================================
// Observation
// ============================================================================

pub export fn forapollo_observe_dispatch(
    obs_id: i32,
    n: i32,
    x: [*]const f64,
    m: i32,
    t: f64,
    obs_params: [*]const f64,
    nop: i32,
    z_pred: [*]f64,
    info: *i32,
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.observeDispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info);
}

pub export fn forapollo_observe_jacobian(
    obs_id: i32,
    n: i32,
    x: [*]const f64,
    m: i32,
    t: f64,
    obs_params: [*]const f64,
    nop: i32,
    H: [*]f64,
    info: *i32,
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.observeJacobian(obs_id, n, x, m, t, obs_params, nop, H, info);
}

// ============================================================================
// Propagation
// ============================================================================

pub export fn forapollo_propagate(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, n_steps, info);
}

pub export fn forapollo_propagate_stm(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.propagateStm(n, x, phi, u, nu, f_ptr, df_ptr, model_id, params, np, dt, n_steps, info);
}

pub export fn forapollo_propagate_batch(
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
) callconv(.c) void {
    if (n <= 0 or n_states <= 0) {
        // Set all info codes to invalid input
        var i: i32 = 0;
        while (i < n_states) : (i += 1) {
            info_batch[@intCast(i)] = 3;
        }
        return;
    }
    dispatch.propagateBatch(n, n_states, x_batch, u, nu, f_ptr, model_id, params, np, dt, n_steps, info_batch);
}

// ============================================================================
// Kalman Filter
// ============================================================================

pub export fn forapollo_kf_predict(
    n: i32,
    x: [*]f64,
    P: [*]f64,
    F: [*]const f64,
    Q: [*]const f64,
    info: *i32,
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.kfPredict(n, x, P, F, Q, info);
}

pub export fn forapollo_kf_update(
    n: i32,
    m: i32,
    x: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    H: [*]const f64,
    R: [*]const f64,
    validity: [*]i32,
    info: *i32,
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.kfUpdate(n, m, x, P, z, H, R, validity, info);
}

// ============================================================================
// Extended Kalman Filter
// ============================================================================

pub export fn forapollo_ekf_predict(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.ekfPredict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, n_steps, info);
}

pub export fn forapollo_ekf_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.ekfUpdate(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info);
}

// ============================================================================
// Iterated Extended Kalman Filter
// ============================================================================

pub export fn forapollo_iekf_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.iekfUpdate(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, max_iter, tol, validity, info);
}

// ============================================================================
// Unscented Kalman Filter
// ============================================================================

pub export fn forapollo_ukf_predict(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.ukfPredict(n, x, P, f_ptr, model_id, params, np, Q, dt, alpha, beta_ukf, kappa, info);
}

pub export fn forapollo_ukf_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.ukfUpdate(n, m, x, P, z, h_ptr, obs_id, R, obs_params, nop, alpha, beta_ukf, kappa, validity, info);
}

// ============================================================================
// Error-State Kalman Filter
// ============================================================================

pub export fn forapollo_eskf_predict(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.eskfPredict(n, x_nom, dx, P, f_ptr, model_id, params, np, Q, dt, n_steps, info);
}

pub export fn forapollo_eskf_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.eskfUpdate(n, m, x_nom, dx, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info);
}

pub export fn forapollo_eskf_inject(
    n: i32,
    x_nom: [*]f64,
    dx: [*]f64,
    P: [*]f64,
    info: *i32,
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.eskfInject(n, x_nom, dx, P, info);
}

// ============================================================================
// Square-Root Extended Kalman Filter
// ============================================================================

pub export fn forapollo_srekf_predict(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.srekfPredict(n, x, S, f_ptr, df_ptr, model_id, params, np, Sq, dt, n_steps, info);
}

pub export fn forapollo_srekf_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.srekfUpdate(n, m, x, S, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info);
}

// ============================================================================
// Square-Root Unscented Kalman Filter
// ============================================================================

pub export fn forapollo_srukf_predict(
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
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.srukfPredict(n, x, S, f_ptr, model_id, params, np, Sq, dt, alpha, beta_ukf, kappa, info);
}

pub export fn forapollo_srukf_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.srukfUpdate(n, m, x, S, z, h_ptr, obs_id, R, obs_params, nop, alpha, beta_ukf, kappa, validity, info);
}

// ============================================================================
// Information Filter
// ============================================================================

pub export fn forapollo_if_predict(
    n: i32,
    eta: [*]f64,
    Ymat: [*]f64,
    F: [*]const f64,
    Q: [*]const f64,
    info: *i32,
) callconv(.c) void {
    if (n <= 0) {
        info.* = 3;
        return;
    }
    dispatch.ifPredict(n, eta, Ymat, F, Q, info);
}

pub export fn forapollo_if_update(
    n: i32,
    m: i32,
    eta: [*]f64,
    Ymat: [*]f64,
    z: [*]const f64,
    H: [*]const f64,
    R: [*]const f64,
    validity: [*]i32,
    info: *i32,
) callconv(.c) void {
    if (n <= 0 or m <= 0) {
        info.* = 3;
        return;
    }
    dispatch.ifUpdate(n, m, eta, Ymat, z, H, R, validity, info);
}

// ============================================================================
// Particle Filter (SIR)
// ============================================================================

pub export fn forapollo_pf_sir_predict(
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
) callconv(.c) void {
    if (n <= 0 or n_particles <= 0) {
        info.* = 3;
        return;
    }
    dispatch.pfSirPredict(n, n_particles, particles, weights, f_ptr, model_id, params, np, Q, dt, n_steps, seed, info);
}

pub export fn forapollo_pf_sir_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0 or n_particles <= 0) {
        info.* = 3;
        return;
    }
    dispatch.pfSirUpdate(n, m, n_particles, particles, weights, z, h_ptr, obs_id, R, obs_params, nop, info);
}

pub export fn forapollo_pf_sir_resample(
    n: i32,
    n_particles: i32,
    particles: [*]f64,
    weights: [*]f64,
    seed: *i32,
    info: *i32,
) callconv(.c) void {
    if (n <= 0 or n_particles <= 0) {
        info.* = 3;
        return;
    }
    dispatch.pfSirResample(n, n_particles, particles, weights, seed, info);
}

pub export fn forapollo_pf_rb_update(
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
) callconv(.c) void {
    if (n <= 0 or m <= 0 or n_particles <= 0) {
        info.* = 3;
        return;
    }
    dispatch.pfRbUpdate(n, m, n_particles, n_linear, particles, weights, z, h_ptr, obs_id, R, obs_params, nop, info);
}

// ============================================================================
// Smoothers
// ============================================================================

pub export fn forapollo_rts_smooth(
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
) callconv(.c) void {
    if (n <= 0 or nsteps <= 0) {
        info.* = 3;
        return;
    }
    dispatch.rtsSmooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth, info);
}

pub export fn forapollo_urtss_smooth(
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
) callconv(.c) void {
    if (n <= 0 or nsteps <= 0) {
        info.* = 3;
        return;
    }
    dispatch.urtssSmooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth, info);
}

// ============================================================================
// Batch Estimators
// ============================================================================

pub export fn forapollo_batch_wls(
    n: i32,
    m_total: i32,
    x: [*]f64,
    z_all: [*]const f64,
    H_all: [*]const f64,
    R_all: [*]const f64,
    max_iter: i32,
    tol: f64,
    info: *i32,
) callconv(.c) void {
    if (n <= 0 or m_total <= 0) {
        info.* = 3;
        return;
    }
    dispatch.batchWls(n, m_total, x, z_all, H_all, R_all, max_iter, tol, info);
}

pub export fn forapollo_batch_map(
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
) callconv(.c) void {
    if (n <= 0 or m_total <= 0) {
        info.* = 3;
        return;
    }
    dispatch.batchMap(n, m_total, x, x0, P0, z_all, H_all, R_all, max_iter, tol, info);
}

// ============================================================================
// Guidance
// ============================================================================

pub export fn forapollo_guidance_zem(n: i32, x: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidanceZem(n, x, t_go, g, a_cmd, info);
}

pub export fn forapollo_guidance_zev(n: i32, x: [*]const f64, v_target: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidanceZev(n, x, v_target, t_go, g, a_cmd, info);
}

pub export fn forapollo_guidance_eguidance(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidanceEguidance(n, x, x_target, t_go, g, a_cmd, info);
}

pub export fn forapollo_guidance_pn_pure(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidancePnPure(n, x_p, x_t, N_gain, a_cmd, info);
}

pub export fn forapollo_guidance_pn_aug(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_target: [*]const f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidancePnAug(n, x_p, x_t, N_gain, a_target, a_cmd, info);
}

pub export fn forapollo_guidance_pn_true(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidancePnTrue(n, x_p, x_t, N_gain, a_cmd, info);
}

pub export fn forapollo_guidance_lambert(r1: [*]const f64, r2: [*]const f64, tof: f64, mu: f64, v1: [*]f64, v2: [*]f64, n_rev: i32, info: *i32) callconv(.c) void {
    if (tof <= 0.0 or mu <= 0.0) { info.* = 3; return; }
    dispatch.guidanceLambert(r1, r2, tof, mu, v1, v2, n_rev, info);
}

pub export fn forapollo_guidance_lqr(n: i32, m: i32, A: [*]const f64, B: [*]const f64, Q_cost: [*]const f64, R_cost: [*]const f64, x: [*]const f64, u_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n <= 0 or m <= 0) { info.* = 3; return; }
    dispatch.guidanceLqr(n, m, A, B, Q_cost, R_cost, x, u_cmd, info);
}

pub export fn forapollo_guidance_pure_pursuit(x_pos: [*]const f64, x_lookahead: [*]const f64, L_wheelbase: f64, steer_cmd: *f64, info: *i32) callconv(.c) void {
    if (L_wheelbase <= 0.0) { info.* = 3; return; }
    dispatch.guidancePurePursuit(x_pos, x_lookahead, L_wheelbase, steer_cmd, info);
}

pub export fn forapollo_guidance_stanley(x_pos: [*]const f64, path_point: [*]const f64, path_heading: f64, v: f64, k_gain: f64, steer_cmd: *f64, info: *i32) callconv(.c) void {
    dispatch.guidanceStanley(x_pos, path_point, path_heading, v, k_gain, steer_cmd, info);
}

pub export fn forapollo_guidance_traj_track(n: i32, x: [*]const f64, x_ref: [*]const f64, v_ref: [*]const f64, K_pos: f64, K_vel: f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidanceTrajTrack(n, x, x_ref, v_ref, K_pos, K_vel, a_cmd, info);
}

pub export fn forapollo_guidance_min_energy(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, a_cmd: [*]f64, info: *i32) callconv(.c) void {
    if (n < 6) { info.* = 3; return; }
    dispatch.guidanceMinEnergy(n, x, x_target, t_go, a_cmd, info);
}

// ============================================================================
// Coordinates
// ============================================================================

pub export fn forapollo_coords_transform(from_id: i32, to_id: i32, n: i32, x_in: [*]const f64, x_out: [*]f64, t: f64, params: [*]const f64, np: i32, info: *i32) callconv(.c) void {
    if (n < 3) { info.* = 3; return; }
    dispatch.coordsTransform(from_id, to_id, n, x_in, x_out, t, params, np, info);
}

pub export fn forapollo_coords_eci_to_ecef(x_eci: [*]const f64, x_ecef: [*]f64, gmst: f64, info: *i32) callconv(.c) void {
    dispatch.coordsEciToEcef(x_eci, x_ecef, gmst, info);
}

pub export fn forapollo_coords_ecef_to_eci(x_ecef: [*]const f64, x_eci: [*]f64, gmst: f64, info: *i32) callconv(.c) void {
    dispatch.coordsEcefToEci(x_ecef, x_eci, gmst, info);
}

pub export fn forapollo_coords_ecef_to_geodetic(x_ecef: [*]const f64, lat: *f64, lon: *f64, alt: *f64, a_body: f64, f_body: f64, info: *i32) callconv(.c) void {
    if (a_body <= 0.0) { info.* = 3; return; }
    dispatch.coordsEcefToGeodetic(x_ecef, lat, lon, alt, a_body, f_body, info);
}

pub export fn forapollo_coords_geodetic_to_ecef(lat: f64, lon: f64, alt: f64, x_ecef: [*]f64, a_body: f64, f_body: f64, info: *i32) callconv(.c) void {
    if (a_body <= 0.0) { info.* = 3; return; }
    dispatch.coordsGeodeticToEcef(lat, lon, alt, x_ecef, a_body, f_body, info);
}

pub export fn forapollo_coords_cart_to_keplerian(rv: [*]const f64, mu: f64, oe: [*]f64, info: *i32) callconv(.c) void {
    if (mu <= 0.0) { info.* = 3; return; }
    dispatch.coordsCartToKeplerian(rv, mu, oe, info);
}

pub export fn forapollo_coords_keplerian_to_cart(oe: [*]const f64, mu: f64, rv: [*]f64, info: *i32) callconv(.c) void {
    if (mu <= 0.0) { info.* = 3; return; }
    dispatch.coordsKeplerianToCart(oe, mu, rv, info);
}

// ============================================================================
// Astrodynamics
// ============================================================================

pub export fn forapollo_astro_kepler_solve(mean_anom: f64, ecc: f64, ecc_anom: *f64, info: *i32) callconv(.c) void {
    if (ecc < 0.0 or ecc >= 1.0) { info.* = 3; return; }
    dispatch.astroKeplerSolve(mean_anom, ecc, ecc_anom, info);
}

pub export fn forapollo_astro_vis_viva(r: f64, a: f64, mu: f64, v: *f64, info: *i32) callconv(.c) void {
    if (r <= 0.0 or mu <= 0.0) { info.* = 3; return; }
    dispatch.astroVisViva(r, a, mu, v, info);
}

pub export fn forapollo_astro_period(a: f64, mu: f64, T: *f64, info: *i32) callconv(.c) void {
    if (a <= 0.0 or mu <= 0.0) { info.* = 3; return; }
    dispatch.astroPeriod(a, mu, T, info);
}

pub export fn forapollo_astro_hohmann(r1: f64, r2: f64, mu: f64, dv1: *f64, dv2: *f64, tof: *f64, info: *i32) callconv(.c) void {
    if (r1 <= 0.0 or r2 <= 0.0 or mu <= 0.0) { info.* = 3; return; }
    dispatch.astroHohmann(r1, r2, mu, dv1, dv2, tof, info);
}

pub export fn forapollo_astro_planetary_mu(body_id: i32, mu: *f64, info: *i32) callconv(.c) void {
    dispatch.astroPlanetaryMu(body_id, mu, info);
}

// ============================================================================
// Environment
// ============================================================================

pub export fn forapollo_environ_atmosphere_us76(h: f64, rho: *f64, T_atm: *f64, p_atm: *f64, info: *i32) callconv(.c) void {
    if (h < 0.0 or h > 86.0) { info.* = 3; return; }
    dispatch.environAtmosphereUs76(h, rho, T_atm, p_atm, info);
}

pub export fn forapollo_environ_gravity_j2(r_vec: [*]const f64, mu: f64, J2: f64, R_eq: f64, g_vec: [*]f64, info: *i32) callconv(.c) void {
    if (mu <= 0.0 or R_eq <= 0.0) { info.* = 3; return; }
    dispatch.environGravityJ2(r_vec, mu, J2, R_eq, g_vec, info);
}

pub export fn forapollo_environ_geodesic_vincenty(lat1: f64, lon1: f64, lat2: f64, lon2: f64, a_body: f64, f_body: f64, dist: *f64, az1: *f64, az2: *f64, info: *i32) callconv(.c) void {
    if (a_body <= 0.0) { info.* = 3; return; }
    dispatch.environGeodesicVincenty(lat1, lon1, lat2, lon2, a_body, f_body, dist, az1, az2, info);
}

// ============================================================================
// Time
// ============================================================================

pub export fn forapollo_time_gmst(ut1_jd: f64, gmst_rad: *f64, info: *i32) callconv(.c) void {
    dispatch.timeGmst(ut1_jd, gmst_rad, info);
}

pub export fn forapollo_time_cal_to_jd(year: i32, month: i32, day: i32, hour: i32, minute: i32, second: f64, jd: *f64, info: *i32) callconv(.c) void {
    dispatch.timeCalToJd(year, month, day, hour, minute, second, jd, info);
}

pub export fn forapollo_time_utc_to_tai(utc_jd: f64, tai_jd: *f64, info: *i32) callconv(.c) void {
    dispatch.timeUtcToTai(utc_jd, tai_jd, info);
}
