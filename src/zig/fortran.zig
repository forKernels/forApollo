// Copyright The Fantastic Planet — By David Clabaugh
//
// fortran.zig — extern fn declarations for all fa_* Fortran bind(C) symbols
// These are resolved at link time against libforapollo_fortran.a
//
// Convention:
//   Fortran `value` attribute → pass-by-value (i32, f64)
//   Fortran without `value`  → pass-by-pointer ([*]f64, *i32, [*]i32)
//   Fortran `type(c_funptr), value` → ?*const anyopaque

// ============================================================================
// Dynamics (forapollo_dynamics.f90)
// ============================================================================

pub extern "c" fn fa_dynamics_dispatch(
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
) void;

pub extern "c" fn fa_dynamics_jacobian(
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
) void;

// ============================================================================
// Observation (forapollo_observe.f90)
// ============================================================================

pub extern "c" fn fa_observe_dispatch(
    obs_id: i32,
    n: i32,
    x: [*]const f64,
    m: i32,
    t: f64,
    obs_params: [*]const f64,
    nop: i32,
    z_pred: [*]f64,
    info: *i32,
) void;

pub extern "c" fn fa_observe_jacobian(
    obs_id: i32,
    n: i32,
    x: [*]const f64,
    m: i32,
    t: f64,
    obs_params: [*]const f64,
    nop: i32,
    H: [*]f64,
    info: *i32,
) void;

// ============================================================================
// Propagation (forapollo_propagate.f90)
// ============================================================================

pub extern "c" fn fa_propagate(
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
) void;

pub extern "c" fn fa_propagate_stm(
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
) void;

pub extern "c" fn fa_propagate_batch(
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
) void;

// ============================================================================
// Kalman Filter (forapollo_estimate.f90)
// ============================================================================

pub extern "c" fn fa_kf_predict(
    n: i32,
    x: [*]f64,
    P: [*]f64,
    F: [*]const f64,
    Q: [*]const f64,
    info: *i32,
) void;

pub extern "c" fn fa_kf_update(
    n: i32,
    m: i32,
    x: [*]f64,
    P: [*]f64,
    z: [*]const f64,
    H: [*]const f64,
    R: [*]const f64,
    validity: [*]i32,
    info: *i32,
) void;

// ============================================================================
// Extended Kalman Filter
// ============================================================================

pub extern "c" fn fa_ekf_predict(
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
) void;

pub extern "c" fn fa_ekf_update(
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
) void;

// ============================================================================
// Iterated Extended Kalman Filter
// ============================================================================

pub extern "c" fn fa_iekf_update(
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
) void;

// ============================================================================
// Unscented Kalman Filter
// ============================================================================

pub extern "c" fn fa_ukf_predict(
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
) void;

pub extern "c" fn fa_ukf_update(
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
) void;

// ============================================================================
// Error-State Kalman Filter
// ============================================================================

pub extern "c" fn fa_eskf_predict(
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
) void;

pub extern "c" fn fa_eskf_update(
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
) void;

pub extern "c" fn fa_eskf_inject(
    n: i32,
    x_nom: [*]f64,
    dx: [*]f64,
    P: [*]f64,
    info: *i32,
) void;

// ============================================================================
// Square-Root Extended Kalman Filter
// ============================================================================

pub extern "c" fn fa_srekf_predict(
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
) void;

pub extern "c" fn fa_srekf_update(
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
) void;

// ============================================================================
// Square-Root Unscented Kalman Filter
// ============================================================================

pub extern "c" fn fa_srukf_predict(
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
) void;

pub extern "c" fn fa_srukf_update(
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
) void;

// ============================================================================
// Information Filter
// ============================================================================

pub extern "c" fn fa_if_predict(
    n: i32,
    eta: [*]f64,
    Ymat: [*]f64,
    F: [*]const f64,
    Q: [*]const f64,
    info: *i32,
) void;

pub extern "c" fn fa_if_update(
    n: i32,
    m: i32,
    eta: [*]f64,
    Ymat: [*]f64,
    z: [*]const f64,
    H: [*]const f64,
    R: [*]const f64,
    validity: [*]i32,
    info: *i32,
) void;

// ============================================================================
// Particle Filter (SIR)
// ============================================================================

pub extern "c" fn fa_pf_sir_predict(
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
) void;

pub extern "c" fn fa_pf_sir_update(
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
) void;

pub extern "c" fn fa_pf_sir_resample(
    n: i32,
    n_particles: i32,
    particles: [*]f64,
    weights: [*]f64,
    seed: *i32,
    info: *i32,
) void;

pub extern "c" fn fa_pf_rb_update(
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
) void;

// ============================================================================
// Smoothers
// ============================================================================

pub extern "c" fn fa_rts_smooth(
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
) void;

pub extern "c" fn fa_urtss_smooth(
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
) void;

// ============================================================================
// Batch Estimators
// ============================================================================

pub extern "c" fn fa_batch_wls(
    n: i32,
    m_total: i32,
    x: [*]f64,
    z_all: [*]const f64,
    H_all: [*]const f64,
    R_all: [*]const f64,
    max_iter: i32,
    tol: f64,
    info: *i32,
) void;

pub extern "c" fn fa_batch_map(
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
) void;
