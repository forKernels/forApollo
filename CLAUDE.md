# forApollo

Universal state estimation, navigation, and guidance engine. Part of the forKernels ecosystem.

The same Kalman filter that guided Apollo to the Moon tracks a portfolio's drift. The same Lambert solver that computes transfer orbits finds optimal paths through any state space. The math doesn't care about the domain.

**"The state vector doesn't care what it represents."**

## Architecture

Zortran three-layer pattern (forApollo is the **canonical reference model** for the [Zortran Wiring Contract](../Wintermute/docs/architecture/zortran-wiring-contract.html)):

- **Layer 1 — Fortran kernel (internal).** Pure math in `src/fortran/*.f90`. `bind(C, name="fa_*")` uses the **internal-only** `fa_` prefix that **no consumer is permitted to bind**.
- **Layer 2 — Producer Zig export (the sole public C-ABI).** `src/zig/exports.zig` declares `export fn fapo_op(...) callconv(.c)`. It validates, marshals, dispatches, and **calls the Layer-1 kernel via `extern`**. This Zig compiles into `libforapollo.a`.
- **Layer 3 — Consumer (Python / Wintermute / siblings).** Links the prebuilt `.a` and binds the **Zig export** `fapo_*` only — never an `fa_*` Fortran symbol. Compiles zero forApollo Zig.

Fortran does the math behind a closed door; Zig owns the boundary. No C code.

```
src/fortran/          ← Fortran kernels (.f90), all bind(C, name="fa_*")
src/zig/              ← Zig dispatch, safe wrappers, exports
deps/                 ← Prebuilt upstream .a archives + Zig modules (committed)
prebuilt/             ← Prebuilt forApollo .a per platform (committed)
python/               ← Python bindings (ctypes → Zig safety layer)
```

> **Delivery (CANON, updated 2026-07-02):** `zig build` assembles `prebuilt/<target>/libforapollo.a` via `tools/wrap.zig` — the Zig exports object + this repo's Fortran `.o` are `ld -r` combined and everything except `fapo_*` is internalized (`fa_*` and the `mat_*/vec_*` helpers become local). The wrap has a built-in gate: any public non-`forapollo_` symbol fails the build. `<target>` is one of the short names `macos | thor | linX86 | winX86`.

### Functional Tiers (distinct from the wiring layers above)

| Tier | Files | Purpose |
|-------|-------|---------|
| **Engine** | `estimate`, `propagate`, `guidance`, `coords` | Domain-agnostic core. Any state vector. |
| **Models** | `dynamics`, `observe` | Pluggable catalog of built-in dynamics & measurement models with analytic Jacobians. |
| **Domain** | `astro`, `environ`, `time` | Space-specific utilities. Non-space users never touch these. |

### Call Chain

```
Consumer (Python ctypes / Wintermute / sibling) → fapo_* (Layer 2: Zig export, sole public C-ABI) → fa_* (Layer 1: Fortran kernels, internal)
```

No consumer ever binds an `fa_*` symbol — `fa_*` is internal-only and lives inside `libforapollo.a`. Every call enters through a Zig `export fn fapo_*`, which validates, bounds-checks, marshals, and then calls the Fortran kernel via `extern`. The math stays in Fortran (no pure-Zig reimplementation — Anti-Pattern A); Zig owns the public boundary (Fortran never exposes `bind(C, name="fapo_*")` — Anti-Pattern B).

## Naming Conventions

- Fortran symbols (Layer 1, **internal-only** — no consumer binds these): `fa_` prefix (e.g., `fa_ekf_predict`, `fa_lambert_solve`)
- Fortran files: `forapollo_*.f90`
- Zig exports (Layer 2, **the sole public C-ABI**): `fapo_*` prefix — what Python/Wintermute/siblings link
- NEVER the retired `fk_` prefix; forApollo's internal prefix is `fa_`, its public prefix is `fapo_` (short-prefix canon 2026-07-02; `forapollo_` retired)
- Module structure: one `.f90` per thematic group, matching its Zig wrapper file

## Matrix Convention

All matrices are **flat 1D arrays in row-major order** at the C ABI boundary:
- `P(n*n)` — covariance, flat row-major
- `Q(n*n)` — process noise, flat row-major
- `F(n*n)` — dynamics Jacobian, flat row-major
- `H(m*n)` — measurement Jacobian, flat row-major

## Fortran Source Files (9 total)

All in `src/fortran/`. All `bind(C, name="fa_*")`.

### ENGINE (domain-agnostic):

**forapollo_estimate.f90** — 16 estimator algorithms:
- KF, EKF, IEKF, ESKF (error-state — Apollo heritage), UKF
- SR-EKF, SR-UKF (square-root, numerically stable)
- IF, EIF (information filters — multi-sensor fusion)
- SIR particle, auxiliary particle, Rao-Blackwellized particle
- RTS smoother, unscented smoother (offline/batch)
- Batch WLS, batch MAP (orbit determination)

**forapollo_propagate.f90** — state propagation:
- `fa_propagate` — advance state by dt using dynamics f (pointer or model_id)
- `fa_propagate_stm` — propagate state + state transition matrix
- `fa_propagate_batch` — OpenMP batch: propagate N states in parallel
- Dispatches to forMath ODE solvers (RK4/DP5/DP8/DOP853). Never reimplements RK.

**forapollo_guidance.f90** — 23 guidance/control laws:
- Zero-effort: ZEM/ZEV (Apollo P63/P64), E-guidance
- Proportional nav: pure, augmented, true PN
- Polynomial: gravity turn, linear tangent, PEG
- Optimal control: LQR, iLQR, DDP
- MPC: shooting, direct collocation
- Targeting: single/multi shooting, Lambert, differential correction
- Path following: pure pursuit, Stanley, trajectory tracking
- Energy-optimal: min-fuel, min-time, min-energy

**forapollo_coords.f90** — coordinate frame catalog:
- Inertial: ECI (J2000), ICRF, generic
- Rotating: ECEF, Moon-fixed, body-fixed
- Local: NED, ENU, LVLH
- Orbital: PQW, RSW, VNC
- Geodetic: WGS84, generic ellipsoid
- Topocentric, pinhole camera, user-defined
- Calls forMath for rotation math (quaternions, SO(3) exp/log)

### MODELS (pluggable catalog):

**forapollo_dynamics.f90** — 17 built-in dynamics models:
- Orbital: Kepler, J2, CR3BP, atmospheric drag
- Rigid body: 6-DOF quaternion with forces/torques
- Ground: bicycle, Ackermann, differential drive
- Aerial: quadrotor 12-state, fixed-wing 6-DOF
- Tracking: constant-velocity, constant-accel, constant-turn
- Stochastic: geometric Brownian motion, Ornstein-Uhlenbeck
- Scalar: double integrator, spring-mass-damper
- Each model ships with analytic Jacobian

**forapollo_observe.f90** — 16 built-in measurement models:
- Position (direct, range-only, bearing-only, range+bearing)
- Velocity (direct, Doppler)
- Attitude (magnetometer, star tracker, sun sensor)
- Inertial (accelerometer, gyroscope)
- Radar (range + azimuth + elevation)
- Scalar, pinhole camera, relative position/velocity
- Each model ships with analytic Jacobian

### DOMAIN (space-specific):

**forapollo_astro.f90** — astrodynamics utilities:
- JPL ephemeris (DE405/DE430), planetary constants
- Eclipse geometry, ground track
- Hohmann/bi-elliptic transfers
- Orbital element conversions (classical, equinoctial, Cartesian)
- Vis-viva, period, sphere of influence

**forapollo_environ.f90** — environment models:
- US Standard Atmosphere 1976, exponential atmosphere
- J2/J4/spherical harmonic gravity (uses forFFT)
- Solar radiation pressure
- Geodesics (Vincenty, Karney)

**forapollo_time.f90** — precision time systems (future: migrates to forTime):
- UTC, TAI, TT, TDB, GPS, MJD, JD, Unix epoch
- Leap second table, relativistic corrections

## Hybrid Dispatch (function pointer + built-in catalog)

```fortran
! Custom dynamics — user supplies function pointer
call fa_ekf_predict(n, x, P, my_f_ptr, my_df_ptr, Q, dt, 0, params, np, info)

! Built-in — pass null, select by model ID
call fa_ekf_predict(n, x, P, c_null_funptr, c_null_funptr, Q, dt, FA_DYN_KEPLER, params, np, info)
```

Built-in models provide free analytic Jacobians. Custom models fall back to forMath numerical differentiation.

## Ternary Measurement Gating (forTernary)

Every measurement update passes through a three-valued validity gate:
- +1: fuse normally
-  0: uncertain — prediction only
- -1: reject & flag

Also used for: estimator health, mode detection, multi-sensor consistency.

## Kernel Style

```fortran
subroutine fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, info) &
    bind(C, name="fa_ekf_predict")
    use iso_c_binding
    integer(c_int), value :: n, model_id, np
    real(c_double), intent(inout) :: x(n), P(n*n)  ! flat row-major
    type(c_funptr), value :: f_ptr      ! dynamics (null → use model_id)
    type(c_funptr), value :: df_ptr     ! Jacobian (null → auto)
    real(c_double), intent(in) :: Q(n*n), dt, params(np)
    integer(c_int), intent(out) :: info  ! 0=ok, 1=diverged, 2=singular, 3=invalid
    ! Domain-agnostic: works for 3D position, 6D orbit, 100D portfolio
end subroutine
```

## Dependencies

Links prebuilt `.a` archives + Zig modules. No upstream `.f90` in this repo.

| Dependency | What forApollo uses | Prefix |
|------------|-------------------|--------|
| **forMath** | Linear algebra, ODE solvers, quaternions, Lie groups, special functions, random, numdiff | `fm_` |
| **forFFT** | Spectral methods, gravity harmonic expansion | `ff_` |
| **forOpt** | Heavy optimization for MPC, trajectory targeting | `fo_` |
| **forTernary** | Three-valued logic for sensor gating, mode detection, estimator health | `fk_ternary_` |
| **forGraph** | Graph search for path planning, multi-target assignment, mission phase DAG | `fgr_` |
| **forTime** | Time system conversions (future, when forTime ships) | `ft_` |

### DEPENDENCY RULES — READ BEFORE WRITING CODE

1. NEVER add another library's .f90 files to this build. Link their prebuilt .a archive.
2. NEVER rewrite Zig wrappers that exist in upstream libraries. Import them: `@import("formath")`.
3. If you need a symbol from forMath/forFFT/forOpt/etc: check prebuilt archive, import Zig module, do NOT copy source or write new extern fn.
4. The only .f90 files in this repo are THIS library's kernels.
5. Default build links prebuilt deps. Source build opt-in: `-Ddev=true`.

## Build

```bash
zig build                          # build library (links prebuilt deps)
zig build test                     # run tests
zig build -Duse-prebuilt=true      # use prebuilt forApollo objects
zig build -Dgenerate-prebuilt=true # regenerate prebuilt objects
zig build -Ddev=true               # source-build deps (development only)
zig build -Dtarget=aarch64-linux-gnu  # cross-compile for Jetson/RPi
```

## Error Codes

All `info` outputs: 0=ok, 1=diverged, 2=singular, 3=invalid input, 4=convergence failure, 5=integration failure.

## Thread Safety

All `fa_*` routines are stateless and reentrant. No module-level mutable state.

## OpenMP

Batch variants use `!$omp parallel do`. Single-point estimation is sequential. Dispatch threshold ~100 elements.

## Edge Deployment

All kernels run on Jetson (aarch64) with zero external dependencies. The same EKF that tracks a spacecraft tracks a robot arm joint — different f(x) and h(x), identical algorithm.

## Sources & Heritage

- Apollo AGC flight software (Colossus/Luminary — guidance, navigation, ESKF)
- Fortran Astrodynamics Toolkit (jacobwilliams — orbital mechanics, BSD)
- NASA CFDTOOLS (coordinate transforms)
- SPICE concepts (JPL/NAIF — ephemeris, reference frames, time)
- Lee & Wright 2014 (block-tridiagonal solver)

All reimplemented in modern Fortran 2008+ with iso_c_binding. No legacy code vendored.

Copyright The Fantastic Planet — By David Clabaugh
