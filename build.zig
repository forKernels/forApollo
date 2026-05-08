// forApollo — Zig Build System (Stage 2)
// Copyright The Fantastic Planet — By David Clabaugh
//
// Per FORKERNELS delivery standard (per-target branch):
//   - This is the `winX86` branch; delivery dir matches branch name.
//   - Output: zig-out/winX86/{lib,bin,include}/forapollo.lib
//   - Sibling resolution: ../<sibling>/zig-out/winX86/lib/
//     (each sibling repo is also on its own winX86 branch)
//   - Prefix: forapollo_
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

fn getTargetName(t: std.Target) []const u8 {
    return switch (t.os.tag) {
        .macos => switch (t.cpu.arch) {
            .aarch64 => "macos-arm64",
            else => "macos-unknown",
        },
        .linux => switch (t.cpu.arch) {
            .aarch64 => "linux-arm64",
            .x86_64 => "linux-x86_64",
            else => "linux-unknown",
        },
        .windows => switch (t.cpu.arch) {
            .x86_64 => "windows-x86_64",
            else => "windows-unknown",
        },
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
    const optimize = b.standardOptimizeOption(.{});
    const target_name = getTargetName(target.result);

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

    // Build-time options passed into Zig source
    const build_opts = b.addOptions();
    build_opts.addOption(bool, "dev", dev);
    build_opts.addOption(bool, "use_prebuilt", use_prebuilt);

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

    const make_step = b.addSystemCommand(&.{ "make", "lib" });

    // -----------------------------------------------------------------------
    // Static library: libforapollo.a (Zig-only, no external linking)
    // -----------------------------------------------------------------------

    const root_module = b.createModule(.{
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });
    root_module.addOptions("build_opts", build_opts);

    const static_lib = b.addLibrary(.{
        .linkage = .static,
        .name = "forapollo",
        .root_module = root_module,
    });
    static_lib.linkLibC();

    if (!use_prebuilt) {
        static_lib.step.dependOn(&make_step.step);
    }

    {
        const install = b.addInstallArtifact(static_lib, .{});
        b.getInstallStep().dependOn(&install.step);
    }

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

        const install = b.addInstallArtifact(shared_lib, .{});
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
