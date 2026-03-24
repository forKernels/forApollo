// Copyright The Fantastic Planet — By David Clabaugh
//
// kernels.zig — Zig-native implementations for small state dimensions
//
// These avoid FFI overhead for the most common tracking problems (n <= 6).
// Currently a placeholder — all dispatch goes through Fortran via dispatch.zig.
//
// Future optimization plan:
//   - Zig-native KF predict/update for n=4 (2D position + velocity tracking)
//   - Zig-native KF predict/update for n=6 (3D position + velocity tracking)
//   - Zig-native KF predict/update for n=9 (3D pos + vel + accel)
//
// The math is identical to the Fortran but avoids the C ABI call overhead,
// and can use Zig's comptime for fixed-size matrix operations.
//
// When implemented, dispatch.zig will route small-n calls here instead
// of through fortran.zig. The export API in exports.zig does not change.

// Placeholder to ensure the module compiles and can be imported by dispatch.zig.
// Zig-native estimator kernels will be implemented here as the project matures.
pub fn placeholder() void {
    // Zig-native estimators will be implemented here.
    // For now, all calls go through fortran.zig.
}
