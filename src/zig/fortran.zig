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

// ============================================================================
// Guidance (forapollo_guidance.f90)
// ============================================================================

// Zero-effort guidance (Apollo heritage)
pub extern "c" fn fa_guidance_zem(n: i32, x: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_zev(n: i32, x: [*]const f64, v_target: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_eguidance(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void;

// Proportional navigation
pub extern "c" fn fa_guidance_pn_pure(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_pn_aug(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_target: [*]const f64, a_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_pn_true(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_cmd: [*]f64, info: *i32) void;

// Polynomial guidance
pub extern "c" fn fa_guidance_gravity_turn(t: f64, t_burn: f64, x0: [*]const f64, g: f64, thrust_mag: f64, mass: f64, x_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_linear_tangent(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, g: [*]const f64, a_mag: f64, a_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_peg(n: i32, x: [*]const f64, x_target: [*]const f64, v_exhaust: f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void;

// Optimal control
pub extern "c" fn fa_guidance_lqr(n: i32, m: i32, A: [*]const f64, B: [*]const f64, Q_cost: [*]const f64, R_cost: [*]const f64, x: [*]const f64, u_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_ilqr(n: i32, m: i32, N_horizon: i32, x_traj: [*]f64, u_traj: [*]f64, Q_cost: [*]const f64, R_cost: [*]const f64, Qf: [*]const f64, x_ref: [*]const f64, dt: f64, max_iter: i32, tol: f64, info: *i32) void;
pub extern "c" fn fa_guidance_ddp(n: i32, m: i32, N_horizon: i32, x_traj: [*]f64, u_traj: [*]f64, Q_cost: [*]const f64, R_cost: [*]const f64, Qf: [*]const f64, x_ref: [*]const f64, dt: f64, max_iter: i32, tol: f64, alpha: f64, info: *i32) void;

// MPC
pub extern "c" fn fa_guidance_mpc_shooting(n: i32, m: i32, N_horizon: i32, x0: [*]const f64, u_traj: [*]f64, Q_cost: [*]const f64, R_cost: [*]const f64, Qf: [*]const f64, x_ref: [*]const f64, dt: f64, max_iter: i32, tol: f64, info: *i32) void;
pub extern "c" fn fa_guidance_mpc_collocation(n: i32, m: i32, N_horizon: i32, x_traj: [*]f64, u_traj: [*]f64, Q_cost: [*]const f64, R_cost: [*]const f64, Qf: [*]const f64, x_ref: [*]const f64, dt: f64, max_iter: i32, tol: f64, info: *i32) void;

// Targeting
pub extern "c" fn fa_guidance_lambert(r1: [*]const f64, r2: [*]const f64, tof: f64, mu: f64, v1: [*]f64, v2: [*]f64, n_rev: i32, info: *i32) void;
pub extern "c" fn fa_guidance_single_shooting(n: i32, m: i32, x0: [*]const f64, x_target: [*]const f64, u_guess: [*]const f64, dt: f64, N_steps: i32, max_iter: i32, tol: f64, u_result: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_multi_shooting(n: i32, m: i32, N_seg: i32, x_nodes: [*]f64, u_nodes: [*]f64, dt_seg: f64, max_iter: i32, tol: f64, info: *i32) void;
pub extern "c" fn fa_guidance_diffcorr(n: i32, m: i32, x0: [*]f64, x_target: [*]const f64, x_free_idx: [*]const i32, n_free: i32, dt: f64, max_iter: i32, tol: f64, info: *i32) void;

// Path following
pub extern "c" fn fa_guidance_pure_pursuit(x_pos: [*]const f64, x_lookahead: [*]const f64, L_wheelbase: f64, steer_cmd: *f64, info: *i32) void;
pub extern "c" fn fa_guidance_stanley(x_pos: [*]const f64, path_point: [*]const f64, path_heading: f64, v: f64, k_gain: f64, steer_cmd: *f64, info: *i32) void;
pub extern "c" fn fa_guidance_traj_track(n: i32, x: [*]const f64, x_ref: [*]const f64, v_ref: [*]const f64, K_pos: f64, K_vel: f64, a_cmd: [*]f64, info: *i32) void;

// Energy-optimal
pub extern "c" fn fa_guidance_min_fuel(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, v_exhaust: f64, a_max: f64, a_cmd: [*]f64, info: *i32) void;
pub extern "c" fn fa_guidance_min_energy(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, a_cmd: [*]f64, info: *i32) void;

// ============================================================================
// Coordinates (forapollo_coords.f90)
// ============================================================================

pub extern "c" fn fa_coords_transform(from_id: i32, to_id: i32, n: i32, x_in: [*]const f64, x_out: [*]f64, t: f64, params: [*]const f64, np: i32, info: *i32) void;
pub extern "c" fn fa_coords_rotation(from_id: i32, to_id: i32, t: f64, params: [*]const f64, np: i32, R: [*]f64, info: *i32) void;
pub extern "c" fn fa_coords_eci_to_ecef(x_eci: [*]const f64, x_ecef: [*]f64, gmst: f64, info: *i32) void;
pub extern "c" fn fa_coords_ecef_to_eci(x_ecef: [*]const f64, x_eci: [*]f64, gmst: f64, info: *i32) void;
pub extern "c" fn fa_coords_ecef_to_geodetic(x_ecef: [*]const f64, lat: *f64, lon: *f64, alt: *f64, a_body: f64, f_body: f64, info: *i32) void;
pub extern "c" fn fa_coords_geodetic_to_ecef(lat: f64, lon: f64, alt: f64, x_ecef: [*]f64, a_body: f64, f_body: f64, info: *i32) void;
pub extern "c" fn fa_coords_ecef_to_ned(x_ecef: [*]const f64, x_ned: [*]f64, lat_ref: f64, lon_ref: f64, alt_ref: f64, a_body: f64, f_body: f64, info: *i32) void;
pub extern "c" fn fa_coords_ned_to_ecef(x_ned: [*]const f64, x_ecef: [*]f64, lat_ref: f64, lon_ref: f64, alt_ref: f64, a_body: f64, f_body: f64, info: *i32) void;
pub extern "c" fn fa_coords_ecef_to_enu(x_ecef: [*]const f64, x_enu: [*]f64, lat_ref: f64, lon_ref: f64, alt_ref: f64, a_body: f64, f_body: f64, info: *i32) void;
pub extern "c" fn fa_coords_enu_to_ecef(x_enu: [*]const f64, x_ecef: [*]f64, lat_ref: f64, lon_ref: f64, alt_ref: f64, a_body: f64, f_body: f64, info: *i32) void;
pub extern "c" fn fa_coords_eci_to_lvlh(x_eci: [*]const f64, x_lvlh: [*]f64, r_ref: [*]const f64, v_ref: [*]const f64, info: *i32) void;
pub extern "c" fn fa_coords_lvlh_to_eci(x_lvlh: [*]const f64, x_eci: [*]f64, r_ref: [*]const f64, v_ref: [*]const f64, info: *i32) void;
pub extern "c" fn fa_coords_cart_to_keplerian(rv: [*]const f64, mu: f64, oe: [*]f64, info: *i32) void;
pub extern "c" fn fa_coords_keplerian_to_cart(oe: [*]const f64, mu: f64, rv: [*]f64, info: *i32) void;
pub extern "c" fn fa_coords_cart_to_equinoctial(rv: [*]const f64, mu: f64, eq: [*]f64, info: *i32) void;
pub extern "c" fn fa_coords_equinoctial_to_cart(eq: [*]const f64, mu: f64, rv: [*]f64, info: *i32) void;
pub extern "c" fn fa_coords_rsw_matrix(r: [*]const f64, v: [*]const f64, R_rsw: [*]f64, info: *i32) void;
pub extern "c" fn fa_coords_vnc_matrix(r: [*]const f64, v: [*]const f64, R_vnc: [*]f64, info: *i32) void;

// ============================================================================
// Astrodynamics (forapollo_astro.f90)
// ============================================================================

pub extern "c" fn fa_astro_kepler_solve(mean_anom: f64, ecc: f64, ecc_anom: *f64, info: *i32) void;
pub extern "c" fn fa_astro_true_anomaly(ecc_anom: f64, ecc: f64, nu: *f64, info: *i32) void;
pub extern "c" fn fa_astro_mean_anomaly(ecc_anom: f64, ecc: f64, mean_anom: *f64, info: *i32) void;
pub extern "c" fn fa_astro_vis_viva(r: f64, a: f64, mu: f64, v: *f64, info: *i32) void;
pub extern "c" fn fa_astro_period(a: f64, mu: f64, T: *f64, info: *i32) void;
pub extern "c" fn fa_astro_soi(a_orbit: f64, m_body: f64, m_central: f64, r_soi: *f64, info: *i32) void;
pub extern "c" fn fa_astro_hohmann(r1: f64, r2: f64, mu: f64, dv1: *f64, dv2: *f64, tof: *f64, info: *i32) void;
pub extern "c" fn fa_astro_bielliptic(r1: f64, r2: f64, r_int: f64, mu: f64, dv1: *f64, dv2: *f64, dv3: *f64, tof: *f64, info: *i32) void;
pub extern "c" fn fa_astro_eclipse_conical(pos_sat: [*]const f64, pos_sun: [*]const f64, pos_body: [*]const f64, rad_body: f64, shadow: *f64, info: *i32) void;
pub extern "c" fn fa_astro_ground_track(r_eci: [*]const f64, gmst: f64, lat: *f64, lon: *f64, info: *i32) void;
pub extern "c" fn fa_astro_planetary_mu(body_id: i32, mu: *f64, info: *i32) void;
pub extern "c" fn fa_astro_planetary_radius(body_id: i32, radius: *f64, info: *i32) void;
pub extern "c" fn fa_astro_stumpff_c2(psi: f64, c2: *f64, info: *i32) void;
pub extern "c" fn fa_astro_stumpff_c3(psi: f64, c3: *f64, info: *i32) void;
pub extern "c" fn fa_astro_universal_kepler(r0: f64, vr0: f64, alpha: f64, dt: f64, mu: f64, chi: *f64, info: *i32) void;

// ============================================================================
// Environment (forapollo_environ.f90)
// ============================================================================

pub extern "c" fn fa_environ_atmosphere_exp(h: f64, rho0: f64, h_scale: f64, rho: *f64, info: *i32) void;
pub extern "c" fn fa_environ_atmosphere_us76(h: f64, rho: *f64, T_atm: *f64, p_atm: *f64, info: *i32) void;
pub extern "c" fn fa_environ_gravity_pointmass(r_vec: [*]const f64, mu: f64, g_vec: [*]f64, info: *i32) void;
pub extern "c" fn fa_environ_gravity_j2(r_vec: [*]const f64, mu: f64, J2: f64, R_eq: f64, g_vec: [*]f64, info: *i32) void;
pub extern "c" fn fa_environ_gravity_j4(r_vec: [*]const f64, mu: f64, J2: f64, J4: f64, R_eq: f64, g_vec: [*]f64, info: *i32) void;
pub extern "c" fn fa_environ_srp(r_sat: [*]const f64, r_sun: [*]const f64, A_over_m: f64, Cr: f64, a_srp: [*]f64, info: *i32) void;
pub extern "c" fn fa_environ_geodesic_vincenty(lat1: f64, lon1: f64, lat2: f64, lon2: f64, a_body: f64, f_body: f64, dist: *f64, az1: *f64, az2: *f64, info: *i32) void;
pub extern "c" fn fa_environ_geodesic_haversine(lat1: f64, lon1: f64, lat2: f64, lon2: f64, R_body: f64, dist: *f64, info: *i32) void;
pub extern "c" fn fa_environ_magnetic_dipole(r_vec: [*]const f64, m_dipole: [*]const f64, B_vec: [*]f64, info: *i32) void;

// ============================================================================
// Time (forapollo_time.f90)
// ============================================================================

pub extern "c" fn fa_time_jd_to_mjd(jd: f64, mjd: *f64, info: *i32) void;
pub extern "c" fn fa_time_mjd_to_jd(mjd: f64, jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_unix_to_jd(unix_sec: f64, jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_jd_to_unix(jd: f64, unix_sec: *f64, info: *i32) void;
pub extern "c" fn fa_time_utc_to_tai(utc_jd: f64, tai_jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_tai_to_utc(tai_jd: f64, utc_jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_tai_to_tt(tai_jd: f64, tt_jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_tt_to_tai(tt_jd: f64, tai_jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_tt_to_tdb(tt_jd: f64, tdb_jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_utc_to_gps(utc_jd: f64, gps_sec: *f64, info: *i32) void;
pub extern "c" fn fa_time_gps_to_utc(gps_sec: f64, utc_jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_gmst(ut1_jd: f64, gmst_rad: *f64, info: *i32) void;
pub extern "c" fn fa_time_cal_to_jd(year: i32, month: i32, day: i32, hour: i32, minute: i32, second: f64, jd: *f64, info: *i32) void;
pub extern "c" fn fa_time_jd_to_cal(jd: f64, year: *i32, month: *i32, day: *i32, hour: *i32, minute: *i32, second: *f64, info: *i32) void;
pub extern "c" fn fa_time_leap_seconds(utc_jd: f64, dt_ls: *f64, info: *i32) void;
