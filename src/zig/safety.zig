// Copyright The Fantastic Planet — By David Clabaugh
//
// safety.zig — bounds checking, error mapping, validation
// Called by exports.zig before dispatching to Fortran.
//
// This layer enforces preconditions that the Fortran kernels assume but
// do not always check. Catching bad inputs here produces clear error codes
// instead of segfaults or silent corruption.

pub const Error = error{
    InvalidDimension,
    NullPointer,
    SingularMatrix,
    Diverged,
    ConvergenceFailed,
    IntegrationFailed,
    InvalidInput,
};

/// Validate state and measurement dimensions.
/// n must be > 0; m must be >= 0.
pub fn validateDimensions(n: i32, m: i32) Error!void {
    if (n <= 0) return Error.InvalidDimension;
    if (m < 0) return Error.InvalidDimension;
}

/// Validate that a required pointer is non-null.
pub fn validatePointer(ptr: ?[*]f64) Error!void {
    if (ptr == null) return Error.NullPointer;
}

/// Validate a const pointer is non-null.
pub fn validateConstPointer(ptr: ?[*]const f64) Error!void {
    if (ptr == null) return Error.NullPointer;
}

/// Validate a mutable i32 pointer is non-null.
pub fn validateIntPointer(ptr: ?[*]i32) Error!void {
    if (ptr == null) return Error.NullPointer;
}

/// Map Fortran info codes to Zig errors.
/// info = 0 means success; nonzero maps to specific error types.
pub fn mapInfoCode(info: i32) Error!void {
    return switch (info) {
        0 => {},
        1 => Error.Diverged,
        2 => Error.SingularMatrix,
        3 => Error.InvalidInput,
        4 => Error.ConvergenceFailed,
        5 => Error.IntegrationFailed,
        else => Error.InvalidInput,
    };
}

/// Convert a Zig error back to a Fortran-style info code for C ABI return.
pub fn errorToInfoCode(err: Error) i32 {
    return switch (err) {
        Error.InvalidDimension => 3,
        Error.NullPointer => 3,
        Error.SingularMatrix => 2,
        Error.Diverged => 1,
        Error.ConvergenceFailed => 4,
        Error.IntegrationFailed => 5,
        Error.InvalidInput => 3,
    };
}

/// Validate estimator predict arguments: n > 0, x/P/Q non-null.
pub fn validateEstimatorPredict(n: i32, x: ?[*]f64, P: ?[*]f64, Q: ?[*]const f64) Error!void {
    try validateDimensions(n, 0);
    try validatePointer(x);
    try validatePointer(P);
    try validateConstPointer(Q);
}

/// Validate estimator update arguments: n > 0, m > 0, x/P/z/R non-null.
pub fn validateEstimatorUpdate(
    n: i32,
    m: i32,
    x: ?[*]f64,
    P: ?[*]f64,
    z: ?[*]const f64,
    R: ?[*]const f64,
) Error!void {
    try validateDimensions(n, m);
    if (m <= 0) return Error.InvalidDimension;
    try validatePointer(x);
    try validatePointer(P);
    try validateConstPointer(z);
    try validateConstPointer(R);
}

/// Validate propagation arguments: n > 0, x non-null.
pub fn validatePropagate(n: i32, x: ?[*]f64) Error!void {
    try validateDimensions(n, 0);
    try validatePointer(x);
}

/// Validate batch propagation: n > 0, n_states > 0, x_batch non-null.
pub fn validateBatchPropagate(n: i32, n_states: i32, x_batch: ?[*]f64) Error!void {
    try validateDimensions(n, 0);
    if (n_states <= 0) return Error.InvalidDimension;
    try validatePointer(x_batch);
}

/// Validate particle filter arguments.
pub fn validateParticleFilter(n: i32, n_particles: i32, particles: ?[*]f64, weights: ?[*]f64) Error!void {
    try validateDimensions(n, 0);
    if (n_particles <= 0) return Error.InvalidDimension;
    try validatePointer(particles);
    try validatePointer(weights);
}

/// Validate smoother arguments.
pub fn validateSmoother(n: i32, nsteps: i32) Error!void {
    try validateDimensions(n, 0);
    if (nsteps <= 0) return Error.InvalidDimension;
}

/// Validate batch estimator arguments.
pub fn validateBatchEstimator(n: i32, m_total: i32, x: ?[*]f64) Error!void {
    try validateDimensions(n, 0);
    if (m_total <= 0) return Error.InvalidDimension;
    try validatePointer(x);
}

/// Validate guidance arguments: n >= 6 for state-based guidance.
pub fn validateGuidance(n: i32) Error!void {
    if (n < 6) return Error.InvalidDimension;
}

/// Validate time-to-go is positive and nonzero.
pub fn validateTimeToGo(t_go: f64) Error!void {
    if (t_go == 0.0 or t_go != t_go) return Error.InvalidInput; // NaN check
}

/// Validate positive scalar (mu, radius, etc).
pub fn validatePositiveScalar(val: f64) Error!void {
    if (val <= 0.0 or val != val) return Error.InvalidInput;
}

/// Validate coordinate frame ID is in known range.
pub fn validateFrameId(id: i32) Error!void {
    // Valid ranges: 1-2 (inertial), 10-12 (rotating), 20-22 (local),
    // 30-32 (orbital), 40-41 (geodetic), 50-51 (other)
    if (id <= 0) return Error.InvalidInput;
}
