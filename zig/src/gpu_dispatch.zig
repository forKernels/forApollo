//! forApollo GPU Dispatch — batched EKF on GPU
const device = @import("device.zig");
const gpu = @import("bindings_gpu.zig");
const c_int = @import("std").c.c_int;

pub fn tryEkfPredictBatch(x: *anyopaque, P: *anyopaque, F: *anyopaque, Q: *anyopaque, n_filters: i32, state_dim: i32) bool {
    if (!device.hasGpu() or n_filters < 16) return false;
    gpu.fa_ekf_predict_batch_gpu(x, P, F, Q, @intCast(n_filters), @intCast(state_dim));
    return true;
}
