// forApollo — Zig Build System (Stage 2)
// Copyright The Fantastic Planet — By David Clabaugh
//
// Per FORKERNELS_BUILD_STANDARD (CANON 2026-06-19):
//   - Output: prebuilt/<branch>/libforapollo.a  (committed, lean)
//   - Branches/targets: macos | thor | linX86 | winX86  (short names only)
//   - Sibling resolution: ../sibling/prebuilt/<branch>/
//   - Public C-ABI prefix: fapo_ (Zig exports, sole public surface; short-prefix
//     canon 2026-07-02, forapollo_ retired); fa_ is the INTERNAL Fortran bind(C)
//     prefix — internalized by tools/wrap.zig
//
// Stage 1: Makefile compiles Fortran -> libforapollo_fortran.a
// Stage 2: This file links Fortran objects + deps -> static + shared libraries
//
// Usage:
//   zig build                         # build static library
//   zig build shared                  # build shared library (needs all deps)
//   zig build -Duse-prebuilt=true     # use prebuilt/libforapollo_fortran.a
//   zig build -Ddev=true              # enable debug assertions and logging

const std = @import("std");

// Branch name = delivery dir name. Each target branch hardcodes its own.
const BRANCH_NAME = "winX86";

// forKernels branch-name delivery canon — winX86/linX86/thor/macos.
fn getTargetName(t: std.Target) []const u8 {
    return switch (t.os.tag) {
        .macos => "macos",
        .linux => switch (t.cpu.arch) {
            .aarch64 => "thor",
            .x86_64 => "linX86",
            else => "unknown",
        },
        .windows => "winX86",
        else => "unknown",
    };
}

// ---------------------------------------------------------------------------
// Dependency library names (linked by name, resolved via search paths)
// ---------------------------------------------------------------------------

/// forMath component libraries needed by forApollo
const formath_libs = [_][]const u8{
    "formath_linalg",
    "formath_ode",
    "formath_quaternion",
    "formath_liegroups",
    "formath_optimize",
    "formath_special",
    "formath_random",
    "formath_numdiff",
    "formath_fft",
};

/// Other upstream dependency libraries
const upstream_libs = [_][]const u8{
    "forfft",
    "foropt",
    "forternary",
    "forGraph",
    "fortime",
};

// ---------------------------------------------------------------------------
// Helper: add sibling library search paths
// ---------------------------------------------------------------------------

fn addSiblingPaths(step: *std.Build.Step.Compile, b: *std.Build, target_name: []const u8) void {
    const siblings = [_]struct { name: []const u8, path: []const u8 }{
        .{ .name = "formath", .path = "../forMath" },
        .{ .name = "forfft", .path = "../forFFT" },
        .{ .name = "foropt", .path = "../forOpt" },
        .{ .name = "forternary", .path = "../forTernary" },
        .{ .name = "forgraph", .path = "../forGraph" },
        .{ .name = "fortime", .path = "../forTime" },
        .{ .name = "forcuda", .path = "../forCUDA" },
    };

    for (siblings) |dep| {
        const sibling = dep.path;
        // Branch-name delivery (forKernels standard): each sibling is on its
        // matching target branch and delivers to zig-out/<branch>/lib.
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/zig-out/" ++ BRANCH_NAME ++ "/lib", .{sibling}) });
        // forMath uses src/zig/zig-out/<branch>/lib for its module libs.
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/src/zig/zig-out/" ++ BRANCH_NAME ++ "/lib", .{sibling}) });
        // Platform-name fallback for siblings not yet on the delivery standard.
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/prebuilt/{s}/lib", .{ sibling, target_name }) });
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/zig-out/{s}/lib", .{ sibling, target_name }) });
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/zig-out/lib", .{sibling}) });
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/prebuilt/lib", .{sibling}) });
    }
}

// ---------------------------------------------------------------------------
// Helper: link all dependencies onto a compile step
// ---------------------------------------------------------------------------

fn linkDeps(
    b: *std.Build,
    step: *std.Build.Step.Compile,
    fortran_obj: []const u8,
    target_name: []const u8,
) void {
    // Fortran kernel archive (stage-1 output)
    step.addObjectFile(b.path(fortran_obj));

    // Sibling search paths
    addSiblingPaths(step, b, target_name);

    // Link forMath component libraries
    for (formath_libs) |lib_name| {
        step.linkSystemLibrary(lib_name);
    }

    // Link other upstream libraries
    for (upstream_libs) |lib_name| {
        step.linkSystemLibrary(lib_name);
    }

    // System runtime libraries
    const tgt = step.rootModuleTarget();
    const is_macos = tgt.os.tag == .macos;
    const is_linux = tgt.os.tag == .linux;
    const is_windows = tgt.os.tag == .windows;
    if (is_macos) {
        step.addLibraryPath(.{ .cwd_relative = "/opt/homebrew/lib/gcc/current" });
        step.addLibraryPath(.{ .cwd_relative = "/opt/homebrew/lib/gcc/15" });
    } else if (is_linux) {
        step.addLibraryPath(.{ .cwd_relative = b.fmt("deps/{s}", .{target_name}) });
        step.addLibraryPath(.{ .cwd_relative = b.fmt("deps/syslibs/{s}", .{target_name}) });
        step.addLibraryPath(.{ .cwd_relative = "/usr/lib" });
        step.addLibraryPath(.{ .cwd_relative = "/usr/lib/aarch64-linux-gnu" });
        step.addLibraryPath(.{ .cwd_relative = "/usr/lib/x86_64-linux-gnu" });
        step.addLibraryPath(.{ .cwd_relative = "/usr/lib/gcc/aarch64-linux-gnu/13" });
        step.addLibraryPath(.{ .cwd_relative = "/usr/lib/gcc/x86_64-linux-gnu/13" });
    } else if (is_windows) {
        // Windows (MSYS2 UCRT64): gfortran runtime — GCC 15.2.0
        step.addLibraryPath(.{ .cwd_relative = "C:/msys64/ucrt64/lib" });
        step.addLibraryPath(.{ .cwd_relative = "C:/msys64/ucrt64/lib/gcc/x86_64-w64-mingw32/15.2.0" });
        // MINGW64 fallback
        step.addLibraryPath(.{ .cwd_relative = "C:/msys64/mingw64/lib" });
        step.addLibraryPath(.{ .cwd_relative = "C:/msys64/mingw64/lib/gcc/x86_64-w64-mingw32/15.2.0" });
    }
    step.linkSystemLibrary("gfortran");
    step.linkSystemLibrary("gomp");
    if (is_windows) {
        step.linkSystemLibrary("quadmath");
        // libgcc_s.a is an import library; addObjectFile because LLD can't process it via -lgcc_s
        step.addObjectFile(.{ .cwd_relative = "C:/msys64/ucrt64/lib/libgcc_s.a" });
        step.linkSystemLibrary("gcc_eh");
    }
    step.linkLibC();
}

// ---------------------------------------------------------------------------
// Build entry point
// ---------------------------------------------------------------------------

pub fn build(b: *std.Build) void {
    // ========================================================================
    // forKernels per-target delivery: zig-out/<branch>/{lib,bin,include}
    // Branch name selects the delivery dir. Downstream repos consume forApollo
    // via, e.g., ../forApollo/zig-out/winX86/lib/forapollo.lib
    //
    // NOTE: addInstallFile (and friends) resolve `.prefix` against
    // `b.install_path`, not `b.install_prefix` — see std/Build.zig
    // getInstallPath. install_prefix is bookkeeping only; install_path is
    // the field every install step actually reads.
    // ========================================================================
    b.install_path = b.pathFromRoot("zig-out/" ++ BRANCH_NAME);
    b.install_prefix = b.install_path;
    b.exe_dir = b.pathJoin(&.{ b.install_path, "bin" });
    b.lib_dir = b.pathJoin(&.{ b.install_path, "lib" });
    b.h_dir = b.pathJoin(&.{ b.install_path, "include" });

    const target = b.standardTargetOptions(.{});
    const optimize = b.option(std.builtin.OptimizeMode, "optimize", "Optimization mode (delivery default: ReleaseFast — Debug materializes 'undefined' as real bytes and shipped a 68MB forCV archive)") orelse .ReleaseFast;
    const target_name = getTargetName(target.result);

    // Canonical delivery (CANON 2026-06-19): prebuilt/<branch>/lib<name>.a
    // <branch> = macos | thor | linX86 | winX86. Short names, NO lib/ segment.
    b.install_path = b.pathFromRoot("prebuilt");
    b.install_prefix = b.install_path;
    b.lib_dir = b.pathFromRoot(b.fmt("prebuilt/{s}", .{target_name}));

    // -----------------------------------------------------------------------
    // Build options
    // -----------------------------------------------------------------------

    const use_prebuilt = b.option(
        bool,
        "use-prebuilt",
        "Link against prebuilt/libforapollo_fortran.a instead of build/",
    ) orelse false;

    const dev = b.option(
        bool,
        "dev",
        "Enable debug assertions and verbose logging",
    ) orelse false;

    // GPU (Thor/Blackwell) kernels are an opt-in extra. They require nvfortran +
    // CUDA runtime and pull undefined CUDA symbols into the archive, which breaks
    // the lean, self-contained CPU delivery. Default OFF; the canonical
    // prebuilt/<branch>/libforapollo.a is CPU-only. fa_ekf_predict_batch_gpu is
    // not referenced by the Zig public ABI, so omitting it leaves no unresolved
    // forapollo_/fa_ symbols.
    const gpu = b.option(
        bool,
        "gpu",
        "Compile and bundle nvfortran CUDA kernels (Thor/Blackwell, opt-in)",
    ) orelse false;

    // Build-time options passed into Zig source
    const build_opts = b.addOptions();
    build_opts.addOption(bool, "dev", dev);
    build_opts.addOption(bool, "use_prebuilt", use_prebuilt);
    build_opts.addOption(bool, "gpu", gpu);

    // -----------------------------------------------------------------------
    // Resolve Fortran archive path (per-target)
    // -----------------------------------------------------------------------

    const fortran_archive: []const u8 = if (use_prebuilt)
        b.fmt("prebuilt/{s}/lib/libforapollo_fortran.a", .{target_name})
    else
        b.fmt("build/{s}/lib/libforapollo_fortran.a", .{target_name});

    // -----------------------------------------------------------------------
    // Stage 1: Build Fortran kernels via make (if not using prebuilt)
    // -----------------------------------------------------------------------

    // Pass TARGET so Stage 1 (Makefile) writes its per-target Fortran objects to
    // the same dir Stage 2 reads (prebuilt/<target>/obj) — even when cross-
    // compiling (-Dtarget=...). All 4 targets build without clobbering.
    const make_step = b.addSystemCommand(&.{ "make", "lib", b.fmt("TARGET={s}", .{target_name}) });

    // -----------------------------------------------------------------------
    // Static library: libforapollo.a — assembled by tools/wrap.zig.
    // Zig exports object + this repo's Fortran .o are dup/MOD-localized,
    // ld -r combined, then EVERYTHING except fapo_* is internalized
    // (fa_* is internal-only per the wiring contract — the wrap enforces it
    // with a gate that fails the build on any leak). Sibling deps (forMath,
    // forCUDA, ...) remain external per the no-bundled-deps canon.
    // -----------------------------------------------------------------------

    // GPU dispatch layer (zig/src/*) is a separate module tree. Expose it under
    // the name "gpu_dispatch" so exports.zig can pull it in via @import without
    // crossing module-path boundaries. Lazily analyzed: only referenced under
    // -Dgpu, so the default CPU build never code-gens it (stays lean).
    const gpu_dispatch_mod = b.createModule(.{
        .root_source_file = b.path("zig/src/gpu_dispatch.zig"),
        .target = target,
        .optimize = optimize,
    });

    const root_module = b.createModule(.{
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    root_module.addOptions("build_opts", build_opts);
    root_module.addImport("gpu_dispatch", gpu_dispatch_mod);

    const exports_obj = b.addObject(.{
        .name = "forapollo_exports",
        .root_module = root_module,
    });

    const wrap_exe = b.addExecutable(.{
        .name = "forapollo-wrap",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tools/wrap.zig"),
            .target = b.graph.host,
            .optimize = .ReleaseSafe,
        }),
    });

    const wrap_run = b.addRunArtifact(wrap_exe);
    wrap_run.has_side_effects = true; // writes prebuilt/<target>/libforapollo.a
    wrap_run.addArg(b.pathFromRoot(b.fmt("prebuilt/{s}/libforapollo.a", .{target_name})));
    wrap_run.addArg(b.pathFromRoot(b.fmt(".zig-cache/wrap/{s}", .{target_name})));
    wrap_run.addArg("*fapo_*");
    wrap_run.addFileArg(exports_obj.getEmittedBin());

    // Fortran objects always live at prebuilt/<target>/obj — the Makefile
    // (Stage 1) writes them there per-target, so the make-output path always
    // equals this Stage-2 read path for the SAME target.
    const fortran_obj_dir = b.fmt("prebuilt/{s}/obj", .{target_name});
    const fortran_kernel_basenames = [_][]const u8{
        "forapollo_dynamics", "forapollo_observe",  "forapollo_estimate",
        "forapollo_propagate", "forapollo_guidance", "forapollo_coords",
        "forapollo_astro",     "forapollo_environ",  "forapollo_time",
    };
    for (fortran_kernel_basenames) |name| {
        wrap_run.addArg(b.pathFromRoot(b.fmt("{s}/{s}.o", .{ fortran_obj_dir, name })));
    }

    if (!use_prebuilt) {
        wrap_run.step.dependOn(&make_step.step);
    }

    // Compile GPU kernels via nvfortran (opt-in via -Dgpu, Thor/Blackwell only)
    if (gpu and target.result.cpu.arch == .aarch64 and target.result.os.tag == .linux) {
        const nvfortran_path = "/opt/nvidia/hpc_sdk/Linux_aarch64/26.3/compilers/bin/nvfortran";
        const repo_gpu_sources = [_][]const u8{  "forapollo_ekf_batch_gpu", };
        for (repo_gpu_sources) |gpu_src| {
            const compile = b.addSystemCommand(&.{
                nvfortran_path, "-c", "-cuda", "-gpu=cc110", "-O3", "-Mfree", "-fPIC",
            });
            compile.addFileArg(b.path(b.fmt("src/gpu/{s}.cuf", .{gpu_src})));
            compile.addArg("-o");
            const obj = compile.addOutputFileArg(b.fmt("{s}.o", .{gpu_src}));
            wrap_run.addFileArg(obj);
        }
        // forCUDA runtime (fc_rt_*) stays an undefined ref in the archive;
        // the consumer's final link resolves it (no-bundled-deps canon).
    }

    b.getInstallStep().dependOn(&wrap_run.step);

    // -----------------------------------------------------------------------
    // Shared library: zig build shared  (links Fortran + all deps)
    // -----------------------------------------------------------------------

    const shared_step = b.step("shared", "Build shared library (requires all sibling deps)");
    {
        const shared_module = b.createModule(.{
            .root_source_file = b.path("src/zig/exports.zig"),
            .target = target,
            .optimize = optimize,
        });
        shared_module.addOptions("build_opts", build_opts);
        shared_module.addImport("gpu_dispatch", gpu_dispatch_mod);

        const shared_lib = b.addLibrary(.{
            .linkage = .dynamic,
            .name = "forapollo",
            .root_module = shared_module,
            .version = .{ .major = 0, .minor = 1, .patch = 0 },
        });
        linkDeps(b, shared_lib, fortran_archive, target_name);

        if (!use_prebuilt) {
            shared_lib.step.dependOn(&make_step.step);
        }

        const install = b.addInstallArtifact(shared_lib, .{
        });
        shared_step.dependOn(&install.step);
    }

    // -----------------------------------------------------------------------
    // Test step
    // -----------------------------------------------------------------------

    const test_module = b.createModule(.{
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });
    test_module.addOptions("build_opts", build_opts);
    test_module.addImport("gpu_dispatch", gpu_dispatch_mod);

    const unit_tests = b.addTest(.{
        .root_module = test_module,
    });
    linkDeps(b, unit_tests, fortran_archive, target_name);

    if (!use_prebuilt) {
        unit_tests.step.dependOn(&make_step.step);
    }

    const run_tests = b.addRunArtifact(unit_tests);
    const test_step = b.step("test", "Run Zig unit tests");
    test_step.dependOn(&run_tests.step);
}
