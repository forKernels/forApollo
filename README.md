<div align="center">

# 🚀 forApollo

### *The math that landed on the Moon — generalized for everything.*

[![Fortran](https://img.shields.io/badge/Fortran-2008%2B-734F96?style=for-the-badge&logo=fortran&logoColor=white)](https://fortran-lang.org)
[![Zig](https://img.shields.io/badge/Zig-0.15.x-F7A41D?style=for-the-badge&logo=zig&logoColor=white)](https://ziglang.org)
[![License](https://img.shields.io/badge/License-Proprietary-red?style=for-the-badge)](#license)
[![Part of](https://img.shields.io/badge/forKernels-Ecosystem-00D4AA?style=for-the-badge)](#forkernels-ecosystem)

<br>

*A universal state estimation, navigation, and guidance engine.*
*Born from Apollo. Built for everything.*

---

<img width="700" alt="forApollo heritage" src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/98/Aldrin_Apollo_11_original.jpg/800px-Aldrin_Apollo_11_original.jpg">

<sub>Buzz Aldrin on the lunar surface, July 20, 1969. The guidance computer that got him there ran at 2 MHz with 74KB of memory. The algorithms inside it are timeless.</sub>

</div>

---

## 🌍 The Story

On July 20, 1969, the Apollo Guidance Computer — a machine slower than a modern calculator — landed two humans on the Moon. It did this with **74 kilobytes of memory**, a **2 MHz clock**, and some of the most elegant algorithms ever written.

Those algorithms didn't care that they were running on a spacecraft.

The **Kalman filter** that tracked Apollo's position? It's the same math that fuses GPS and IMU data on your phone. The **Lambert solver** that computed transfer orbits? It's the same boundary-value problem that plans robot arm trajectories. The **powered descent guidance** that landed Eagle on the Sea of Tranquility? It's the same optimal control that lands a drone on a rooftop.

**The state vector doesn't care what it represents.**

A position in orbit. A stock price. A robot joint angle. A training loss curve. They're all just numbers — and they all need the same thing: *estimation* (where am I?), *prediction* (where will I be?), and *guidance* (how do I get where I want to go?).

**forApollo** takes the flight-proven algorithms from the Apollo Guidance Computer and generalizes them into a universal engine. Same math. Any domain. Any state space.

---

## 🏗️ Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    YOUR APPLICATION                          │
│         Python · Rust · C++ · JavaScript · Zig              │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│               ZIG SAFETY LAYER                              │
│    Bounds checking · Error mapping · Size-based dispatch    │
│    @import("formath") · @import("forternary")              │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│            FORTRAN 2008 KERNELS (fa_*)                      │
│    Pure math · OpenMP parallel · ISO_C_BINDING              │
│    Zero dependencies · Edge-deployable                      │
└─────────────────────────────────────────────────────────────┘
```

### Three Layers, One Engine

| Layer | Files | Purpose |
|-------|-------|---------|
| 🔧 **Engine** | `estimate` · `propagate` · `guidance` · `coords` | Domain-agnostic core. Works with any state vector. |
| 🧩 **Models** | `dynamics` · `observe` | Pluggable catalog of dynamics & measurement models. |
| 🌌 **Domain** | `astro` · `environ` · `time` | Space-specific utilities. Non-space users never touch these. |

---

## 🧠 The Engine (Domain-Agnostic)

The engine doesn't know if it's tracking a spacecraft, a drone, or a stock price. It just sees state vectors, covariance matrices, and dynamics functions.

### Estimators — *"Where am I?"*

| Algorithm | Heritage | What it does |
|-----------|----------|-------------|
| 🟢 **KF** | Kalman 1960 | Optimal for linear systems |
| 🟢 **EKF** | Apollo AGC | Jacobian linearization — the workhorse |
| 🟢 **IEKF** | — | Re-linearizes at updated state for tighter estimates |
| 🔴 **ESKF** | **Apollo flight software** | Error-state filter — *what actually flew to the Moon* |
| 🟣 **UKF** | Julier & Uhlmann 1997 | Sigma-point — no Jacobians needed |
| 🔵 **SR-EKF / SR-UKF** | — | Square-root variants — numerically bulletproof |
| 🟡 **IF / EIF** | — | Information form — natural for multi-sensor fusion |
| 🟠 **SIR Particle** | Gordon et al. 1993 | Bootstrap — handles any nonlinearity |
| 🟠 **Auxiliary Particle** | Pitt & Shephard 1999 | Better proposal distribution |
| 🟠 **Rao-Blackwellized** | — | Hybrid: particles for nonlinear, Kalman for linear substates |
| ⚪ **RTS Smoother** | Rauch et al. 1965 | Backward pass — optimal offline estimation |
| ⚪ **URTSS** | — | Unscented smoother |
| 🔷 **Batch WLS** | — | Weighted least squares (orbit determination) |
| 🔷 **Batch MAP** | — | Maximum a posteriori with prior |

### Guidance — *"How do I get there?"*

| Category | Algorithms |
|----------|-----------|
| 🎯 **Zero-effort** | ZEM/ZEV (Apollo powered descent), E-guidance |
| 🎯 **Proportional nav** | Pure PN, Augmented PN, True PN |
| 🎯 **Polynomial** | Gravity turn, Linear tangent steering, PEG |
| 🎯 **Optimal control** | LQR, iLQR, DDP |
| 🎯 **Model predictive** | MPC (shooting), MPC (collocation) |
| 🎯 **Targeting** | Single/multi shooting, Lambert, Differential correction |
| 🎯 **Path following** | Pure pursuit, Stanley, Trajectory tracking |
| 🎯 **Energy-optimal** | Min-fuel, Min-time, Min-energy transfers |

### Ternary Gating — *"Can I trust this sensor?"*

Every measurement passes through a **forTernary** validity gate before fusion:

```
Sensor → Ternary Gate → +1 (fuse) · 0 (hold, prediction only) · -1 (reject & flag)
```

No more binary thresholds on uncertain data. Three-valued logic propagates uncertainty honestly.

---

## 🧩 The Models (Pluggable Catalog)

Use a built-in model by ID, or supply your own function pointer. Built-in models ship with **analytic Jacobians** — the EKF gets exact derivatives for free.

### Dynamics — *"How does my system evolve?"*

| Domain | Models | State dim |
|--------|--------|-----------|
| 🌌 **Orbital** | Keplerian · J2 · CR3BP · Atmospheric drag | 6-7 |
| 🤖 **Rigid body** | 6-DOF quaternion with forces/torques | 13 |
| 🚗 **Ground** | Bicycle · Ackermann · Differential drive | 4-5 |
| 🚁 **Aerial** | Quadrotor 12-state · Fixed-wing 6-DOF | 12 |
| 📡 **Tracking** | Constant-velocity · Constant-accel · Constant-turn | 2d-5 |
| 📈 **Stochastic** | Geometric Brownian motion · Ornstein-Uhlenbeck | 1 |
| ⚙️ **Scalar** | Double integrator · Spring-mass-damper | 2-2d |

### Measurements — *"What do my sensors see?"*

| Category | Models |
|----------|--------|
| 📍 **Position** | Direct position · Range-only · Bearing-only · Range+bearing |
| 💨 **Velocity** | Direct velocity · Doppler (range-rate) |
| 🧭 **Attitude** | Magnetometer · Star tracker · Sun sensor |
| 📐 **Inertial** | Accelerometer · Gyroscope |
| 📡 **Radar/LiDAR** | Range + azimuth + elevation |
| 🔢 **Scalar** | Direct scalar observation |
| 📷 **Image** | Pinhole camera pixel coordinates |
| 🔗 **Relative** | Relative position/velocity (rendezvous, formation) |

### Custom Models

```fortran
! Custom dynamics — supply your own function pointer
call fa_ekf_predict(n, x, P, my_f_ptr, my_df_ptr, Q, dt, 0, params, np, info)

! Built-in model — pass null, select by ID (free analytic Jacobians)
call fa_ekf_predict(n, x, P, c_null_funptr, c_null_funptr, Q, dt, FA_DYN_KEPLER, params, np, info)
```

---

## 🌌 The Domain Layer (Space-Specific)

For orbital mechanics users. Everyone else can ignore this entirely.

| Module | Contents |
|--------|----------|
| 🪐 **Astro** | JPL ephemeris (DE405/430) · Planetary constants · Eclipse geometry · Hohmann/bi-elliptic transfers · Orbital element conversions · Vis-viva |
| 🌍 **Environment** | US Standard Atmosphere 1976 · Exponential atmosphere · J2/J4/spherical harmonic gravity · Solar radiation pressure · Vincenty/Karney geodesics |
| ⏱️ **Time** | UTC · TAI · TT · TDB · GPS · MJD · JD · Unix epoch · Leap seconds · Relativistic corrections (TDB-TT) |

---

## 🔌 Coordinate Frames

forApollo defines *what* to rotate. forMath handles *how* (quaternions, SO(3), Lie groups).

| Category | Frames |
|----------|--------|
| **Inertial** | ECI (J2000) · ICRF · Generic inertial |
| **Rotating** | ECEF · Moon-fixed · Body-fixed |
| **Local** | NED · ENU · LVLH |
| **Orbital** | Perifocal (PQW) · RSW · VNC |
| **Geodetic** | WGS84 · Generic ellipsoid |
| **Topocentric** | Azimuth/elevation/range |
| **Camera** | Pinhole body-to-camera |
| **Generic** | User-defined rotation/quaternion |

---

## 📦 Dependencies

forApollo links prebuilt archives from the forKernels ecosystem. **No upstream `.f90` files in this repo.**

| Dependency | What forApollo uses |
|------------|-------------------|
| **forMath** | Linear algebra · ODE solvers · Quaternions · Lie groups · Special functions · Random · Numerical differentiation |
| **forFFT** | Spectral methods · Gravity harmonic expansion |
| **forOpt** | Heavy optimization for MPC · Trajectory targeting |
| **forTernary** | Three-valued logic for sensor gating · Mode detection · Estimator health |
| **forGraph** | Graph search for path planning · Multi-target assignment · Mission phase DAG |

```
deps/
├── formath/    { lib/*.a + zig/ }
├── forfft/     { lib/*.a + zig/ }
├── foropt/     { lib/*.a + zig/ }
├── forternary/ { lib/*.a + zig/ }
└── forgraph/   { lib/*.a + zig/ }
```

---

## 🔨 Build

```bash
zig build                          # build library (links prebuilt deps)
zig build test                     # run tests
zig build -Duse-prebuilt=true      # use prebuilt Fortran objects
zig build -Dgenerate-prebuilt=true # regenerate prebuilt objects
zig build -Ddev=true               # source-build deps (development only)
```

### Cross-compilation

```bash
zig build -Dtarget=aarch64-linux-gnu   # Jetson Orin / Raspberry Pi / Cloud
zig build                              # macOS (development)
```

---

## 🎯 Use Cases

The same engine powers all of these:

| Domain | What forApollo does |
|--------|-------------------|
| 🚀 **Spacecraft** | Orbit determination · Powered descent · Rendezvous targeting |
| 🤖 **Robotics** | Sensor fusion · Path planning · Trajectory optimization |
| 🚁 **Drones** | GPS/IMU/barometer fusion · Waypoint guidance · Obstacle avoidance |
| 🚗 **Autonomous vehicles** | Multi-sensor tracking · Lane following · MPC |
| 📡 **Radar tracking** | Multi-target tracking · Track-before-detect |
| 📈 **Finance** | Price process estimation · Regime detection · Mean-reversion tracking |
| 🧠 **ML training** | Loss curve estimation · Learning rate guidance · Convergence prediction |
| ⚡ **IoT/Edge** | Sensor filtering · Anomaly detection · State monitoring |

---

## 📜 Heritage & Sources

forApollo's algorithms are drawn from flight-proven and peer-reviewed sources:

| Source | What we took |
|--------|-------------|
| 🌙 **Apollo AGC** ([Colossus/Luminary](https://github.com/chrislgarry/Apollo-11)) | Kepler propagation · Lambert targeting · Powered descent guidance · ESKF navigation · Attitude maneuvers |
| 📚 **Fortran Astrodynamics Toolkit** (jacobwilliams, BSD) | Orbital element conversions · Coordinate transforms · Time systems |
| 🛰️ **NASA CFDTOOLS** | Coordinate transforms · Grid utilities |
| 🔭 **SPICE** (JPL/NAIF) | Ephemeris concepts · Reference frames · Time system architecture |
| 📖 **Lee & Wright 2014** | Block-tridiagonal solver for spectral codes |

All reimplemented in modern **Fortran 2008+** with `iso_c_binding`. No legacy code vendored.

---

## 🏛️ forKernels Ecosystem

forApollo is part of the **forKernels** high-performance computing ecosystem:

```
Fortran (pure math) → Zig (safety + dispatch) → Any language (user API)
```

| Library | Domain | Relationship |
|---------|--------|-------------|
| **forMath** | Mathematics | forApollo's foundation — linalg, ODE, quaternions |
| **forFFT** | Spectral methods | Gravity harmonic expansion, spectral analysis |
| **forOpt** | Optimization | MPC solvers, trajectory targeting |
| **forTernary** | Three-valued logic | Sensor gating, mode detection, estimator health |
| **forGraph** | Graph algorithms | Path planning, multi-target assignment |
| **forNav** | Robotics navigation | Consumes forApollo for universal estimation |
| **forSim** | Physics simulation | Uses forApollo for rigid body state integration |
| **forCV** | Computer vision | Feeds measurements to forApollo's estimators |
| **forEdge** | Edge robotics | Bundles forApollo + forNav + forCV for deployment |

---

<div align="center">

### *"The math doesn't care what the state vector represents."*

<sub>Copyright © The Fantastic Planet — By David Clabaugh</sub>

<sub>Built with 🧊 Zig + 🔬 Fortran | Zero dependencies | Edge-deployable | OpenMP parallel</sub>

</div>
