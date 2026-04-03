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

// ============================================================================
// Guidance dispatch — all forward to Fortran
// ============================================================================

pub fn guidanceZem(n: i32, x: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_zem(n, x, t_go, g, a_cmd, info);
}
pub fn guidanceZev(n: i32, x: [*]const f64, v_target: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_zev(n, x, v_target, t_go, g, a_cmd, info);
}
pub fn guidanceEguidance(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, g: [*]const f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_eguidance(n, x, x_target, t_go, g, a_cmd, info);
}
pub fn guidancePnPure(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_pn_pure(n, x_p, x_t, N_gain, a_cmd, info);
}
pub fn guidancePnAug(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_target: [*]const f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_pn_aug(n, x_p, x_t, N_gain, a_target, a_cmd, info);
}
pub fn guidancePnTrue(n: i32, x_p: [*]const f64, x_t: [*]const f64, N_gain: f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_pn_true(n, x_p, x_t, N_gain, a_cmd, info);
}
pub fn guidanceLambert(r1: [*]const f64, r2: [*]const f64, tof: f64, mu: f64, v1: [*]f64, v2: [*]f64, n_rev: i32, info: *i32) void {
    fortran.fa_guidance_lambert(r1, r2, tof, mu, v1, v2, n_rev, info);
}
pub fn guidanceLqr(n: i32, m: i32, A: [*]const f64, B: [*]const f64, Q_cost: [*]const f64, R_cost: [*]const f64, x: [*]const f64, u_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_lqr(n, m, A, B, Q_cost, R_cost, x, u_cmd, info);
}
pub fn guidancePurePursuit(x_pos: [*]const f64, x_lookahead: [*]const f64, L_wheelbase: f64, steer_cmd: *f64, info: *i32) void {
    fortran.fa_guidance_pure_pursuit(x_pos, x_lookahead, L_wheelbase, steer_cmd, info);
}
pub fn guidanceStanley(x_pos: [*]const f64, path_point: [*]const f64, path_heading: f64, v: f64, k_gain: f64, steer_cmd: *f64, info: *i32) void {
    fortran.fa_guidance_stanley(x_pos, path_point, path_heading, v, k_gain, steer_cmd, info);
}
pub fn guidanceTrajTrack(n: i32, x: [*]const f64, x_ref: [*]const f64, v_ref: [*]const f64, K_pos: f64, K_vel: f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_traj_track(n, x, x_ref, v_ref, K_pos, K_vel, a_cmd, info);
}
pub fn guidanceMinEnergy(n: i32, x: [*]const f64, x_target: [*]const f64, t_go: f64, a_cmd: [*]f64, info: *i32) void {
    fortran.fa_guidance_min_energy(n, x, x_target, t_go, a_cmd, info);
}

// ============================================================================
// Coordinates dispatch — all forward to Fortran
// ============================================================================

pub fn coordsTransform(from_id: i32, to_id: i32, n: i32, x_in: [*]const f64, x_out: [*]f64, t: f64, params: [*]const f64, np: i32, info: *i32) void {
    fortran.fa_coords_transform(from_id, to_id, n, x_in, x_out, t, params, np, info);
}
pub fn coordsEciToEcef(x_eci: [*]const f64, x_ecef: [*]f64, gmst: f64, info: *i32) void {
    fortran.fa_coords_eci_to_ecef(x_eci, x_ecef, gmst, info);
}
pub fn coordsEcefToEci(x_ecef: [*]const f64, x_eci: [*]f64, gmst: f64, info: *i32) void {
    fortran.fa_coords_ecef_to_eci(x_ecef, x_eci, gmst, info);
}
pub fn coordsEcefToGeodetic(x_ecef: [*]const f64, lat: *f64, lon: *f64, alt: *f64, a_body: f64, f_body: f64, info: *i32) void {
    fortran.fa_coords_ecef_to_geodetic(x_ecef, lat, lon, alt, a_body, f_body, info);
}
pub fn coordsGeodeticToEcef(lat: f64, lon: f64, alt: f64, x_ecef: [*]f64, a_body: f64, f_body: f64, info: *i32) void {
    fortran.fa_coords_geodetic_to_ecef(lat, lon, alt, x_ecef, a_body, f_body, info);
}
pub fn coordsCartToKeplerian(rv: [*]const f64, mu: f64, oe: [*]f64, info: *i32) void {
    fortran.fa_coords_cart_to_keplerian(rv, mu, oe, info);
}
pub fn coordsKeplerianToCart(oe: [*]const f64, mu: f64, rv: [*]f64, info: *i32) void {
    fortran.fa_coords_keplerian_to_cart(oe, mu, rv, info);
}

// ============================================================================
// Astrodynamics dispatch — all forward to Fortran
// ============================================================================

pub fn astroKeplerSolve(mean_anom: f64, ecc: f64, ecc_anom: *f64, info: *i32) void {
    fortran.fa_astro_kepler_solve(mean_anom, ecc, ecc_anom, info);
}
pub fn astroVisViva(r: f64, a: f64, mu: f64, v: *f64, info: *i32) void {
    fortran.fa_astro_vis_viva(r, a, mu, v, info);
}
pub fn astroPeriod(a: f64, mu: f64, T: *f64, info: *i32) void {
    fortran.fa_astro_period(a, mu, T, info);
}
pub fn astroHohmann(r1: f64, r2: f64, mu: f64, dv1: *f64, dv2: *f64, tof: *f64, info: *i32) void {
    fortran.fa_astro_hohmann(r1, r2, mu, dv1, dv2, tof, info);
}
pub fn astroPlanetaryMu(body_id: i32, mu: *f64, info: *i32) void {
    fortran.fa_astro_planetary_mu(body_id, mu, info);
}

// ============================================================================
// Environment dispatch — all forward to Fortran
// ============================================================================

pub fn environAtmosphereUs76(h: f64, rho: *f64, T_atm: *f64, p_atm: *f64, info: *i32) void {
    fortran.fa_environ_atmosphere_us76(h, rho, T_atm, p_atm, info);
}
pub fn environGravityJ2(r_vec: [*]const f64, mu: f64, J2: f64, R_eq: f64, g_vec: [*]f64, info: *i32) void {
    fortran.fa_environ_gravity_j2(r_vec, mu, J2, R_eq, g_vec, info);
}
pub fn environGeodesicVincenty(lat1: f64, lon1: f64, lat2: f64, lon2: f64, a_body: f64, f_body: f64, dist: *f64, az1: *f64, az2: *f64, info: *i32) void {
    fortran.fa_environ_geodesic_vincenty(lat1, lon1, lat2, lon2, a_body, f_body, dist, az1, az2, info);
}

// ============================================================================
// Time dispatch — all forward to Fortran
// ============================================================================

pub fn timeGmst(ut1_jd: f64, gmst_rad: *f64, info: *i32) void {
    fortran.fa_time_gmst(ut1_jd, gmst_rad, info);
}
pub fn timeCalToJd(year: i32, month: i32, day: i32, hour: i32, minute: i32, second: f64, jd: *f64, info: *i32) void {
    fortran.fa_time_cal_to_jd(year, month, day, hour, minute, second, jd, info);
}
pub fn timeUtcToTai(utc_jd: f64, tai_jd: *f64, info: *i32) void {
    fortran.fa_time_utc_to_tai(utc_jd, tai_jd, info);
}
