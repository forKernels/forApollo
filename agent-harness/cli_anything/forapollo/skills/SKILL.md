# forApollo CLI-Anything Skill

## Identity
CLI agent harness for **forApollo** -- universal state estimation, guidance, and navigation engine. The same Kalman filter that guided Apollo to the Moon.

## Capabilities
- **nav**: Kalman filter (KF predict/update), Extended KF (EKF), Unscented KF (UKF), Error-State KF (ESKF -- Apollo heritage), Information filter (IF), particle filter
- **guidance**: Zero-Effort Miss/Velocity (ZEM/ZEV -- Apollo P63/P64), Lambert problem solver (transfer orbits), LQR controller, pure pursuit, Stanley controller, pure/augmented proportional navigation, minimum energy guidance
- **orbital**: ECI-to-ECEF / ECEF-to-ECI coordinate transforms, ECEF-to-geodetic (WGS84), geodetic-to-ECEF

## Usage
```bash
# Guidance: Apollo-style Zero-Effort Miss
cli-anything-forapollo guidance zem -n 3 0,0,1000,0,0,-50 --tgo 20 0,0,-9.81 --json

# Lambert problem: compute transfer orbit
cli-anything-forapollo guidance lambert 7000000,0,0 0,8000000,0 --tof 3600 --json

# Coordinate transform
cli-anything-forapollo orbital ecef-to-geodetic 6378137,0,0 --json

# Navigation: EKF predict
cli-anything-forapollo nav ekf-predict -n 6 0,0,0,1,0,0 <P_flat> <Q_flat> --dt 0.1 --json

# REPL mode (default)
cli-anything-forapollo
```

## Backend
Wraps `forapollo.navigation`, `forapollo.guidance`, `forapollo.orbital` Python bindings calling real `forapollo_*` C ABI symbols from the Zig/Fortran kernel library. 16 estimator algorithms, 23 guidance laws, coordinate frame catalog.

## State Vector Convention
All state vectors and matrices are flat 1D arrays in row-major order at the CLI boundary. Covariance P is n*n flat, transition F is n*n flat, etc.

## Error Codes
All operations return info: 0=ok, 1=diverged, 2=singular, 3=invalid input, 4=convergence failure.
