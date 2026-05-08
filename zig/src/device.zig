//! GPU device detection for forKernels (shared pattern)
//! Checks fc_rt_is_available() from forCUDA at runtime.
//! Caches result — detection runs once, dispatch is zero-cost after.

const std = @import("std");

pub const Device = enum {
    cpu,
    gpu,
};

/// forCUDA runtime detection — returns 1 if GPU is initialized and available
extern fn fc_rt_is_available() callconv(.c) c_int;
/// Stub: always returns 0 when forCUDA is not linked (CPU-only builds)
extern fn fc_rt_is_stub() callconv(.c) c_int;

const c_int = std.c.c_int;

var detected: ?Device = null;

/// Detect GPU availability. Cached after first call.
pub fn detect() Device {
    if (detected) |d| return d;
    const available = fc_rt_is_available();
    const d: Device = if (available == 1) .gpu else .cpu;
    detected = d;
    return d;
}

/// Force a specific device (for testing or override)
pub fn setDevice(d: Device) void {
    detected = d;
}

/// Reset detection cache (re-detect on next call)
pub fn reset() void {
    detected = null;
}

/// Returns true if GPU is available
pub inline fn hasGpu() bool {
    return detect() == .gpu;
}

