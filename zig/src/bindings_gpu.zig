//! GPU kernel bindings for forApollo (nvfortran CUDA Fortran, SM_110)
const c_int = @import("std").c.c_int;
pub extern "C" fn fa_ekf_predict_batch_gpu(d_x: *anyopaque, d_P: *anyopaque, d_F: *anyopaque, d_Q: *anyopaque, n_filters: c_int, state_dim: c_int) void;
