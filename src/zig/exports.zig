// Copyright The Fantastic Planet — By David Clabaugh
//
// exports.zig — C ABI surface for Python/Rust/C++ consumers
//
// The public C-ABI of libforapollo is the `fa_*` surface. The Fortran kernels
// (src/fortran/*.f90, all `bind(C, name="fa_*")`) ARE the public symbols and
// ship directly in the archive. There is no second `forapollo_*` alias layer:
// consumers bind `fa_*` directly. The dispatch / safety helpers below remain
// importable by Zig-source consumers, and the `-Dgpu` fast-path adds one extra
// `fa_*` entry point.
//
// Flow: caller → fa_* (Fortran kernel) [via dispatch/safety glue for Zig users]

const safety = @import("safety.zig");
const fortran = @import("fortran.zig");
// Public re-export so Zig-source consumers can call
// `@import("forapollo").dispatch.<wrapper>(...)`.
pub const dispatch = @import("dispatch.zig");
pub const cuda = @import("cuda.zig");

const build_opts = @import("build_opts");

// Force the compiler to pull in all referenced modules so link errors
// surface at build time rather than at dlopen time.
comptime {
    _ = @import("kernels.zig");
}

// ============================================================================
// GPU fast-path (opt-in -Dgpu, Thor/Blackwell SM_110)
// ============================================================================
//
// The GPU dispatch layer lives in the sibling module tree (zig/src/*): it calls
// forCUDA's fc_rt_* runtime to detect a device, then launches the nvfortran
// CUDA-Fortran kernel fa_ekf_predict_batch_gpu. It is wired ONLY under -Dgpu:
// the @export below is gated by a comptime-known flag, so the default CPU build
// never analyses gpu_dispatch/device/bindings_gpu — no fa_*_gpu kernel and no
// forCUDA fc_rt_* externs leak into the lean, self-contained CPU archive.
const gpu_layer = struct {
    const gpu_dispatch = @import("gpu_dispatch");

    /// Batched EKF predict with GPU routing. Routes to forCUDA if a device is
    /// present and the batch is large enough; sets used_gpu=1 on the GPU path,
    /// else 0 (caller falls back to the CPU fa_ekf_predict path). Exported as
    /// fa_ekf_predict_batch — distinct from the raw on-device kernel
    /// fa_ekf_predict_batch_gpu it dispatches to.
    fn ekfPredictBatchGpu(
        x: *anyopaque,
        P: *anyopaque,
        F: *anyopaque,
        Q: *anyopaque,
        n_filters: i32,
        state_dim: i32,
        used_gpu: *i32,
    ) callconv(.c) void {
        used_gpu.* = if (gpu_dispatch.tryEkfPredictBatch(x, P, F, Q, n_filters, state_dim)) 1 else 0;
    }
};

comptime {
    if (build_opts.gpu) {
        @export(&gpu_layer.ekfPredictBatchGpu, .{
            .name = "fa_ekf_predict_batch",
            .linkage = .strong,
        });
    }
}

// ============================================================================
// Version
// ============================================================================

/// Returns the library version as a packed u32: 0xMMmmpp (major.minor.patch).
/// The one non-kernel symbol with no Fortran counterpart, so it lives here.
pub export fn fa_version() callconv(.c) u32 {
    return 0x000100; // 0.1.0
}
