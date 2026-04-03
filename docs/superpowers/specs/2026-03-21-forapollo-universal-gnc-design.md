# forApollo — Universal State Estimation, Navigation & Guidance Engine

**Date:** 2026-03-21
**Status:** Approved
**Author:** David Clabaugh / The Fantastic Planet

---

## 1. Vision

forApollo is a universal state estimation, navigation, and guidance engine. The core insight: Kalman filtering, trajectory optimization, and guidance laws are domain-agnostic — the math is identical whether the state vector represents a spacecraft position, a stock price, a robot joint angle, or a training loss curve.

Born from Apollo AGC flight software (Colossus/Luminary), generalized for everything.

**"The state vector doesn't care what it represents."**

---

## 2. Architecture

### 2.1 Zortran Stack

```
Python (ctypes) → Zig safety layer → Fortran 2008 kernels
```

- Python calls `forapollo_*` (Zig C ABI exports)
- Zig validates, bounds-checks, then calls `fa_*` (Fortran `bind(C)`)
- Fortran does pure math, OpenMP for batch operations
- No language besides Zig and Fortran touches Fortran

### 2.2 Three-Layer Design

| Layer | Files | Purpose |
|-------|-------|---------|
| **Engine** | `estimate`, `propagate`, `guidance`, `coords` | Domain-agnostic core. Works with any state vector. |
| **Models** | `dynamics`, `observe` | Pluggable catalog of built-in dynamics & measurement models with analytic Jacobians. |
| **Domain** | `astro`, `environ`, `time` | Space-specific utilities. Non-space users never touch these. |

### 2.3 Fortran Source Files (9 total)

All in `src/fortran/`. All symbols use `fa_` prefix. All use `bind(C, name="fa_*")`.

```
ENGINE:
  forapollo_estimate.f90    — 16 estimator algorithms
  forapollo_propagate.f90   — state propagation, STM, batch (OpenMP)
  forapollo_guidance.f90    — 23 guidance/control laws
  forapollo_coords.f90      — coordinate frame catalog & transforms

MODELS:
  forapollo_dynamics.f90    — 17 built-in dynamics models + dispatch
  forapollo_observe.f90     — 16 built-in measurement models + dispatch

DOMAIN:
  forapollo_astro.f90       — ephemeris, orbital elements, transfers
  forapollo_environ.f90     — atmosphere, gravity, SRP, geodesy
  forapollo_time.f90        — time systems (future: moves to forTime)
```

### 2.4 Zig Layer

```
src/zig/
  dispatch.zig   — size-based routing (small state → Zig, large → Fortran)
  fortran.zig    — extern fn declarations for fa_* symbols
  kernels.zig    — Zig-native small-state implementations
  exports.zig    — C ABI exports (forapollo_* for Python/Rust/C++)
  safety.zig     — bounds checking, error mapping, validation
```

### 2.5 Python Bindings

```
python/forapollo/
  __init__.py       — top-level imports
  _lib.py           — ctypes.CDLL loader, symbol declarations
  estimators.py     — EKF, UKF, ESKF, ParticleFilter, etc.
  dynamics.py       — model ID constants + @custom decorator
  sensors.py        — observation ID constants + custom wrapper
  guidance.py       — guidance law wrappers
  coords.py         — frame transforms
  astro.py          — orbital mechanics utilities
  environ.py        — atmosphere, gravity
  time.py           — time conversions (future: wraps forTime)
```

Call chain: `Python ctypes → forapollo_* (Zig) → fa_* (Fortran)`

Python never calls Fortran directly. The `@custom` decorator wraps Python functions as C callbacks via `ctypes.CFUNCTYPE` for user-supplied dynamics/measurements.

---

## 3. Dependencies

forApollo links prebuilt `.a` archives + Zig modules from the forKernels ecosystem. No upstream `.f90` in this repo.

| Dependency | What forApollo uses | Prefix |
|------------|-------------------|--------|
| **forMath** | Linear algebra, ODE solvers (RK4/DP5/DP8), quaternions, Lie groups SO(3)/SE(3), special functions, random, numerical differentiation | `fm_` |
| **forFFT** | Spectral methods, gravity harmonic expansion | `ff_` |
| **forOpt** | Heavy optimization for MPC, trajectory targeting | `fo_` |
| **forTernary** | Three-valued logic for sensor gating, mode detection, estimator health | `fk_ternary_` |
| **forGraph** | Graph search (A*, Dijkstra) for path planning, multi-target assignment, mission phase DAG | `fgr_` |
| **forTime** | Time system conversions (future, when forTime ships) | `ft_` |

### Dependency Rules

1. NEVER add another library's `.f90` files to this build. Link their prebuilt `.a` archive.
2. NEVER rewrite Zig wrappers that exist upstream. `@import("formath")`.
3. If you need a symbol: check the prebuilt archive, import their Zig module, do NOT copy source.
4. The only `.f90` files in this repo are THIS library's kernels.
5. Default build links prebuilt deps. Source build opt-in: `-Ddev=true`.

### Distribution Model

The public repo ships stripped prebuilt `.a` archives of upstream dependencies containing ONLY the symbols forApollo calls. No source, no docs, just the compiled objects and minimal Zig extern declarations.

```
deps/
├── formath/
│   ├── lib/libformath_*.a     ← committed, stripped to forApollo's subset
│   └── zig/formath.zig        ← minimal extern declarations only
├── forfft/
│   ├── lib/libforfft.a
│   └── zig/forfft.zig
├── foropt/
│   ├── lib/libforopt.a
│   └── zig/foropt.zig
├── forternary/
│   ├── lib/libforternary.a
│   └── zig/forternary.zig
└── forgraph/
    ├── lib/libforgraph.a
    └── zig/forgraph.zig
```

Multi-platform prebuilts of forApollo itself:

```
prebuilt/
├── aarch64-linux-gnu/libforapollo.a      ← Jetson / RPi / Cloud
└── aarch64-apple-darwin/libforapollo.a   ← Mac development
```

---

## 4. Ecosystem Relationships

### forApollo provides (downstream consumers):

| Consumer | What they call |
|----------|---------------|
| **forNav** | Universal EKF/UKF with robotics dynamics presets (replaces forNav's hardcoded filters over time) |
| **forSim** | Rigid body state integration, collision prediction |
| **forEdge** | Bundled sensor fusion for edge robotics |
| **forCV** | Visual navigation measurement models |

### Layered relationship with forNav (Approach C):

- forApollo exports generic `fa_ekf_predict(n, x, P, f_ptr, Q)` with function pointers
- forNav wraps with robotics-specific dynamics (`fn_imu_dynamics`, `fn_wheel_odometry`)
- forNav migrates to calling forApollo over time, no breaking changes

---

## 5. Engine Layer — Estimators

### 5.1 Interface Pattern

Every estimator uses consistent signatures:

```fortran
subroutine fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, info) &
    bind(C, name="fa_ekf_predict")
    use iso_c_binding
    integer(c_int), value :: n, model_id, np
    real(c_double), intent(inout) :: x(n), P(n*n)
    type(c_funptr), value :: f_ptr      ! dynamics (null → use model_id)
    type(c_funptr), value :: df_ptr     ! Jacobian (null → auto finite-diff or built-in analytic)
    real(c_double), intent(in) :: Q(n*n), dt, params(np)
    integer(c_int), intent(out) :: info  ! 0=ok, 1=diverged, 2=singular, 3=invalid
end subroutine
```

### 5.2 Hybrid Dispatch

- If `f_ptr /= c_null_funptr`: use user-supplied function pointer
- If `f_ptr == c_null_funptr`: dispatch to built-in catalog via `model_id`
- If `df_ptr == c_null_funptr` and built-in model: use analytic Jacobian
- If `df_ptr == c_null_funptr` and custom model: fallback to forMath `fm_numdiff_*`

### 5.3 Ternary Measurement Gating

Every update step passes measurements through a forTernary validity gate:

```
Sensor → Innovation chi-squared test → forTernary gate
  +1: fuse normally
   0: uncertain — use prediction only
  -1: reject & flag
```

Update signatures return validity per measurement:

```fortran
subroutine fa_ekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, &
                          obs_params, nop, validity, info) &
    bind(C, name="fa_ekf_update")
    use iso_c_binding
    integer(c_int), value :: n, m, obs_id, nop
    real(c_double), intent(inout) :: x(n), P(n*n)       ! flat row-major
    real(c_double), intent(in) :: z(m), R(m*m)           ! flat row-major
    type(c_funptr), value :: h_ptr       ! measurement fn (null → use obs_id)
    type(c_funptr), value :: dh_ptr      ! Jacobian (null → auto)
    real(c_double), intent(in) :: obs_params(nop)
    integer(c_int), intent(out) :: validity(m)  ! -1/0/+1 per measurement (forTernary)
    integer(c_int), intent(out) :: info         ! 0=ok, 1=diverged, 2=singular, 3=invalid
end subroutine
```

Ternary also used for:

| Location | Decision |
|----------|----------|
| Estimator health | Covariance trace → healthy/degrading/diverged |
| Mode detection | Maneuver detector → coasting/uncertain/maneuvering |
| Multi-sensor fusion | Sensor consistency → agree/ambiguous/conflict |

### 5.4 Estimator Catalog (16 algorithms)

**Kalman family:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_kf_predict`, `fa_kf_update` | Linear Kalman filter |
| `fa_ekf_predict`, `fa_ekf_update` | Extended Kalman filter |
| `fa_iekf_update` | Iterated EKF (re-linearizes, `max_iter` param) |
| `fa_eskf_predict`, `fa_eskf_update`, `fa_eskf_inject` | Error-State KF (Apollo heritage) |
| `fa_ukf_predict`, `fa_ukf_update` | Unscented KF (sigma-point) |

**Square-root variants:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_srekf_predict`, `fa_srekf_update` | Square-root EKF (Cholesky factor of P) |
| `fa_srukf_predict`, `fa_srukf_update` | Square-root UKF |

**Information filters:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_if_predict`, `fa_if_update` | Information filter (inverse covariance) |
| `fa_eif_predict`, `fa_eif_update` | Extended information filter |

**Particle filters:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_pf_sir_predict`, `fa_pf_sir_update`, `fa_pf_sir_resample` | SIR (bootstrap) |
| `fa_pf_aux_update` | Auxiliary particle filter |
| `fa_pf_rb_update` | Rao-Blackwellized (particle + Kalman hybrid) |

**Smoothers:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_rts_smooth` | RTS backward smoother |
| `fa_urtss_smooth` | Unscented smoother |

**Batch:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_batch_wls` | Weighted least squares |
| `fa_batch_map` | Maximum a posteriori |

### 5.5 Propagation

| Symbol | Purpose |
|--------|---------|
| `fa_propagate` | Propagate state x(n) by dt using dynamics f (pointer or model_id) |
| `fa_propagate_stm` | Propagate state + state transition matrix (n + n*n) |
| `fa_propagate_batch` | OpenMP batch: propagate N states in parallel |

Internally dispatches to forMath ODE solvers (RK4/DP5/DP8/DOP853) via `@import("formath")`. `integrator_id` parameter selects solver. forApollo never reimplements Runge-Kutta.

---

## 6. Engine Layer — Guidance

### 6.1 Guidance Catalog (23 algorithms)

All take state `x(n)`, target `x_target(n)`, return commanded `u(nu)`.

**Zero-effort:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_zem_zev` | Zero-effort-miss/velocity (Apollo P63/P64 core) |
| `fa_e_guidance` | E-guidance |

**Proportional navigation:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_pn_pure` | Pure proportional navigation |
| `fa_pn_augmented` | Augmented PN |
| `fa_pn_true` | True PN |

**Polynomial guidance:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_gravity_turn` | Gravity turn |
| `fa_linear_tangent` | Linear tangent steering |
| `fa_peg` | Powered explicit guidance |

**Optimal control:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_lqr_solve` | Linear-quadratic regulator |
| `fa_ilqr_solve` | Iterative LQR (nonlinear) |
| `fa_ddp_solve` | Differential dynamic programming |

**Model predictive control:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_mpc_shooting` | MPC via shooting method |
| `fa_mpc_collocation` | MPC via direct collocation |

**Targeting:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_target_single_shoot` | Single shooting |
| `fa_target_multi_shoot` | Multiple shooting |
| `fa_target_lambert` | Lambert-based targeting |
| `fa_target_diffcorr` | Differential correction |

**Path following:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_pursuit_pure` | Pure pursuit |
| `fa_pursuit_stanley` | Stanley controller |
| `fa_track_trajectory` | Trajectory tracking |

**Energy-optimal:**

| Symbol | Algorithm |
|--------|-----------|
| `fa_transfer_minfuel` | Minimum-fuel transfer |
| `fa_transfer_mintime` | Minimum-time transfer |
| `fa_transfer_minenergy` | Minimum-energy transfer |

### 6.2 forGraph Integration

| Use case | Algorithm |
|----------|-----------|
| Waypoint sequencing | A* / Dijkstra via `@import("forgraph")` |
| Multi-target assignment | Bipartite matching |
| Mission phase sequencing | DAG traversal |

---

## 7. Engine Layer — Coordinates

forApollo defines the frame *catalog* (what ECI→ECEF means). forMath handles the rotation *math* (quaternion multiply, SO(3) exp/log).

### 7.1 Frame IDs

```fortran
integer(c_int), parameter :: FA_FRAME_ECI         = 1   ! J2000
integer(c_int), parameter :: FA_FRAME_ICRF        = 2
integer(c_int), parameter :: FA_FRAME_ECEF        = 10
integer(c_int), parameter :: FA_FRAME_MOON_FIXED  = 11
integer(c_int), parameter :: FA_FRAME_BODY_FIXED  = 12
integer(c_int), parameter :: FA_FRAME_NED         = 20
integer(c_int), parameter :: FA_FRAME_ENU         = 21
integer(c_int), parameter :: FA_FRAME_LVLH        = 22
integer(c_int), parameter :: FA_FRAME_PQW         = 30  ! Perifocal
integer(c_int), parameter :: FA_FRAME_RSW         = 31  ! Radial-along-cross
integer(c_int), parameter :: FA_FRAME_VNC         = 32  ! Velocity-normal-conormal
integer(c_int), parameter :: FA_FRAME_GEODETIC    = 40  ! lat/lon/alt
integer(c_int), parameter :: FA_FRAME_TOPO        = 50  ! az/el/range
integer(c_int), parameter :: FA_FRAME_CAMERA      = 60  ! pinhole
integer(c_int), parameter :: FA_FRAME_CUSTOM      = 99  ! user-defined
```

### 7.2 Transform Interface

```fortran
! Convert between any two frames
subroutine fa_frame_transform(from_id, to_id, n, r_in, r_out, &
                               ref_params, nrp, t, info) &
    bind(C, name="fa_frame_transform")
    ! ref_params: epoch, body orientation, station location, etc.
end subroutine

! Get rotation matrix between frames (for batch transforms)
subroutine fa_frame_rotation(from_id, to_id, ref_params, nrp, t, R33, info) &
    bind(C, name="fa_frame_rotation")
    real(c_double), intent(out) :: R33(9)  ! 3x3 DCM, flat
end subroutine

! Geodetic conversions
subroutine fa_geodetic_from_ecef(r_ecef, a, f, lat, lon, alt, info) &
    bind(C, name="fa_geodetic_from_ecef")
    ! a = semi-major axis (6378137.0 for WGS84), f = flattening
end subroutine
```

---

## 8. Models Layer — Dynamics Catalog

### 8.1 Model IDs

```fortran
! Orbital (1-9)
integer(c_int), parameter :: FA_DYN_KEPLER         = 1
integer(c_int), parameter :: FA_DYN_J2             = 2
integer(c_int), parameter :: FA_DYN_CR3BP          = 3
integer(c_int), parameter :: FA_DYN_DRAG           = 4

! Rigid body (10-19)
integer(c_int), parameter :: FA_DYN_RIGIDBODY_6DOF = 10

! Ground vehicles (20-29)
integer(c_int), parameter :: FA_DYN_BICYCLE        = 20
integer(c_int), parameter :: FA_DYN_ACKERMANN      = 21
integer(c_int), parameter :: FA_DYN_DIFFDRIVE      = 22

! Aerial (30-39)
integer(c_int), parameter :: FA_DYN_QUADROTOR      = 30
integer(c_int), parameter :: FA_DYN_FIXEDWING      = 31

! Tracking (40-49)
integer(c_int), parameter :: FA_DYN_CONST_VEL      = 40
integer(c_int), parameter :: FA_DYN_CONST_ACCEL    = 41
integer(c_int), parameter :: FA_DYN_CONST_TURN     = 42

! Stochastic (50-59)
integer(c_int), parameter :: FA_DYN_GBM            = 50
integer(c_int), parameter :: FA_DYN_OU             = 51

! Scalar (60-69)
integer(c_int), parameter :: FA_DYN_DOUBLE_INT     = 60
integer(c_int), parameter :: FA_DYN_SPRING_MASS    = 61
```

### 8.2 State Vector Conventions

| Model | n | State layout | Control u (nu) | Params |
|-------|---|-------------|----------------|--------|
| Kepler | 6 | [rx,ry,rz, vx,vy,vz] | — | mu |
| J2 | 6 | same | — | mu, J2, R_eq |
| CR3BP | 6 | [x,y,z, vx,vy,vz] rotating | — | mu_ratio |
| Drag | 7 | [r(3), v(3), Cd*A/m] | — | mu, rho0, h_scale |
| Rigid body | 13 | [r(3), v(3), q(4), omega(3)] | [F(3), T(3)] | mass, I(3) |
| Bicycle | 4 | [x,y, theta, v] | [accel, steer_rate] | L (wheelbase) |
| Ackermann | 5 | [x,y, theta, v, steer] | [accel, steer_rate] | L |
| Diff drive | 4 | [x,y, theta, v] | [v_left, v_right] | wheel_sep |
| Quadrotor | 12 | [r(3), v(3), angles(3), rates(3)] | [thrust, tau_x, tau_y, tau_z] (nu=4) | mass, I(3) |
| Fixed-wing | 12 | [r(3), v(3), angles(3), rates(3)] | [thrust, ail, elev, rud] (nu=4) | mass, S, Cl, Cd |
| Const-vel | 2*d | [pos(d), vel(d)] | — (nu=0) | — (d inferred from n/2) |
| Const-accel | 3*d | [pos(d), vel(d), acc(d)] | — (nu=0) | — (d inferred from n/3) |
| Const-turn | 5 | [x,y, vx,vy, omega] | — | — |
| GBM | 1 | [S] | — | mu_drift, sigma |
| O-U | 1 | [X] | — | theta, mu, sigma |
| Double-int | 2*d | [pos(d), vel(d)] | [force(d)] | mass |
| Spring-mass | 2 | [x, v] | [F_ext] | k, c, m |

### 8.3 Dispatch + Jacobians

```fortran
! Dispatch dynamics by model_id
subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
    bind(C, name="fa_dynamics_dispatch")
    integer(c_int), value :: model_id, n, nu, np
    real(c_double), intent(in) :: x(n), u(nu), t, params(np)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info  ! 0=ok, 3=invalid model_id

! Analytic Jacobian for built-in model
subroutine fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info) &
    bind(C, name="fa_dynamics_jacobian")
    integer(c_int), value :: model_id, n, nu, np
    real(c_double), intent(in) :: x(n), u(nu), t, params(np)
    real(c_double), intent(out) :: F(n*n)  ! n-by-n, flat row-major
    integer(c_int), intent(out) :: info    ! 0=ok, 3=invalid model_id
```

Numbering gaps intentional: room for future models per category.

---

## 9. Models Layer — Measurement Catalog

### 9.1 Observation IDs

```fortran
integer(c_int), parameter :: FA_OBS_POSITION       = 1
integer(c_int), parameter :: FA_OBS_RANGE           = 2
integer(c_int), parameter :: FA_OBS_BEARING         = 3
integer(c_int), parameter :: FA_OBS_RANGE_BEARING   = 4
integer(c_int), parameter :: FA_OBS_VELOCITY        = 10
integer(c_int), parameter :: FA_OBS_DOPPLER         = 11
integer(c_int), parameter :: FA_OBS_MAGNETOMETER    = 20
integer(c_int), parameter :: FA_OBS_STARTRACKER     = 21
integer(c_int), parameter :: FA_OBS_SUNSENSOR       = 22
integer(c_int), parameter :: FA_OBS_ACCEL           = 30
integer(c_int), parameter :: FA_OBS_GYRO            = 31
integer(c_int), parameter :: FA_OBS_RADAR           = 40
integer(c_int), parameter :: FA_OBS_SCALAR          = 50
integer(c_int), parameter :: FA_OBS_PINHOLE         = 60
integer(c_int), parameter :: FA_OBS_REL_POS         = 70
integer(c_int), parameter :: FA_OBS_REL_VEL         = 71
```

### 9.2 Dispatch + Jacobians

```fortran
subroutine fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info) &
    bind(C, name="fa_observe_dispatch")
    integer(c_int), value :: obs_id, n, m, nop
    real(c_double), intent(in) :: x(n), t, obs_params(nop)
    real(c_double), intent(out) :: z_pred(m)
    integer(c_int), intent(out) :: info  ! 0=ok, 3=invalid obs_id

subroutine fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, info) &
    bind(C, name="fa_observe_jacobian")
    integer(c_int), value :: obs_id, n, m, nop
    real(c_double), intent(in) :: x(n), t, obs_params(nop)
    real(c_double), intent(out) :: H(m*n)  ! m-by-n, flat row-major
    integer(c_int), intent(out) :: info    ! 0=ok, 3=invalid obs_id
```

### 9.3 Observation Parameters

| Model | m | obs_params |
|-------|---|-----------|
| Position | 2-3 | — |
| Range | 1 | station position (3) |
| Bearing | 1 | station position (3) |
| Range+bearing | 2 | station position (3) |
| Velocity | 2-3 | — |
| Doppler | 1 | station pos+vel (6) |
| Magnetometer | 3 | reference field in inertial (3), body quaternion in state |
| Star tracker | 3 | star unit vector in inertial (3) |
| Sun sensor | 3 | sun unit vector in inertial (3) |
| Accel | 3 | gravity in inertial (3), body quaternion in state |
| Gyro | 3 | body angular rate from state |
| Radar | 3 | station position (3) |
| Scalar | 1 | state index to observe |
| Pinhole | 2 | fx, fy, cx, cy (4) intrinsics + body-to-camera quaternion (4) |
| Rel position | 3 | second state x2(n2) |
| Rel velocity | 3 | second state x2(n2) |

---

## 10. Domain Layer

### 10.1 forapollo_astro.f90

| Symbol | Purpose |
|--------|---------|
| `fa_ephem_planet` | JPL ephemeris (DE405/DE430) — planet pos/vel at epoch |
| `fa_planet_mu` | Gravitational parameter by planet ID |
| `fa_planet_radius` | Mean/equatorial radius by planet ID |
| `fa_eclipse_check` | Shadow geometry (umbra/penumbra) |
| `fa_ground_track` | Lat/lon/alt from inertial state |
| `fa_hohmann` | Hohmann transfer delta-v |
| `fa_bielliptic` | Bi-elliptic transfer delta-v |
| `fa_elements_from_rv` | Cartesian → classical elements (a,e,i,Om,w,nu) |
| `fa_rv_from_elements` | Classical elements → Cartesian |
| `fa_elements_equinoctial` | Classical ↔ equinoctial conversions |
| `fa_vis_viva` | Velocity at radius |
| `fa_period` | Orbital period |
| `fa_soi` | Sphere of influence radius |

### 10.2 forapollo_environ.f90

| Symbol | Purpose |
|--------|---------|
| `fa_atmos_us76` | US Standard Atmosphere 1976 |
| `fa_atmos_exp` | Exponential atmosphere |
| `fa_gravity_j2` | J2 perturbation |
| `fa_gravity_j4` | J2+J4 perturbation |
| `fa_gravity_harmonics` | Spherical harmonic gravity (degree/order N, uses forFFT) |
| `fa_srp` | Solar radiation pressure |
| `fa_geodesic_vincenty` | Geodesic distance (Vincenty) |
| `fa_geodesic_karney` | Geodesic distance (Karney) |

### 10.3 forapollo_time.f90 (future: migrates to forTime)

| Symbol | Purpose |
|--------|---------|
| `fa_time_convert` | Convert between time systems by ID pair |
| `fa_utc_to_tai` | UTC → TAI (leap seconds) |
| `fa_tai_to_tt` | TAI → TT (+32.184s) |
| `fa_tt_to_tdb` | TT → TDB (relativistic) |
| `fa_utc_to_gps` | UTC → GPS time |
| `fa_jd_from_cal` | Calendar → Julian date |
| `fa_mjd_from_jd` | JD → MJD |
| `fa_unix_from_utc` | UTC → Unix epoch |
| `fa_leap_seconds` | Leap second table lookup |

When forTime ships with `ft_` prefix, this file is removed and replaced by `@import("fortime")`.

---

## 11. Build System

```bash
zig build                          # links prebuilt deps, compiles forApollo
zig build test                     # run tests
zig build -Duse-prebuilt=true      # use prebuilt forApollo objects
zig build -Dgenerate-prebuilt=true # regenerate prebuilt objects
zig build -Ddev=true               # source-build deps (development only)
zig build -Dtarget=aarch64-linux-gnu  # cross-compile for Jetson/RPi
```

### Build outputs

```
zig-out/lib/libforapollo.a         ← static library
zig-out/lib/libforapollo.so        ← shared library (for Python ctypes)
zig-out/lib/libforapollo.dylib     ← macOS shared library
```

### Link order (linker command-line, leftmost searched first)

```
-lforapollo -lforgraph -lforternary -lforopt -lformath_* -lforfft
```

forApollo FIRST (depends on everything) → forGraph → forTernary → forOpt → forMath → forFFT LAST (depended on by forMath, no deps itself)

---

## 12. OpenMP Strategy

- Single-point estimator operations (predict/update for one state): sequential
- Batch variants (`fa_propagate_batch`, `fa_pf_sir_*` particle operations): `!$omp parallel do`
- Dispatch threshold: ~100 elements for Fortran batch path
- Particle filter resampling: inherently sequential (systematic resampling), but weight computation is parallel

---

## 13. Error Handling

All `info` output codes:

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Diverged (covariance trace exceeds threshold) |
| 2 | Singular matrix (P or R not invertible) |
| 3 | Invalid input (negative dimension, null required pointer) |
| 4 | Convergence failure (IEKF max_iter, targeting) |
| 5 | Integration failure (ODE solver step rejected) |

Zig safety layer catches code 3 before calling Fortran. Fortran reports 0-5. Zig maps to descriptive error in Python.

---

## 13.1 Matrix Storage Convention

**All matrices are flat 1D arrays in row-major order.** This applies everywhere:

- `P(n*n)` — covariance matrix, flat row-major
- `Q(n*n)` — process noise, flat row-major
- `R(m*m)` — measurement noise, flat row-major
- `F(n*n)` — dynamics Jacobian, flat row-major
- `H(m*n)` — measurement Jacobian, flat row-major

Row-major is chosen for C ABI interop (Python/numpy, C++, Rust all expect row-major). Fortran internally may use column-major, but the `bind(C)` interface always presents row-major to callers. The Zig safety layer handles any necessary transpose.

---

## 13.2 Thread Safety

All `fa_*` routines are **stateless and reentrant**. No module-level mutable state. Multiple filters can run concurrently on different state vectors without locking. OpenMP parallelism is used only within batch operations, not across individual estimator calls.

---

## 14. Heritage & Sources

| Source | Algorithms drawn |
|--------|-----------------|
| Apollo AGC (Colossus/Luminary) | Kepler propagation, Lambert targeting, ESKF navigation, powered descent guidance, attitude maneuvers |
| Fortran Astrodynamics Toolkit (jacobwilliams, BSD) | Orbital element conversions, coordinate transforms, time systems |
| NASA CFDTOOLS | Coordinate transforms, grid utilities |
| SPICE (JPL/NAIF) | Ephemeris concepts, reference frames, time system architecture |
| Lee & Wright 2014 | Block-tridiagonal solver |

All reimplemented in modern Fortran 2008+ with `iso_c_binding`. No legacy code vendored.

---

## 15. Use Cases

| Domain | What forApollo does |
|--------|-------------------|
| Spacecraft | Orbit determination, powered descent, rendezvous |
| Robotics | Sensor fusion, path planning, trajectory optimization |
| Drones | GPS/IMU/baro fusion, waypoint guidance |
| Autonomous vehicles | Multi-sensor tracking, lane following, MPC |
| Radar tracking | Multi-target tracking, track-before-detect |
| Finance | Price process estimation, regime detection |
| ML training | Loss curve estimation, LR guidance, convergence prediction |
| IoT/Edge | Sensor filtering, anomaly detection, state monitoring |

---

Copyright The Fantastic Planet — By David Clabaugh
