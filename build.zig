// forApollo — Zig Build System (Stage 2)
// Copyright The Fantastic Planet — By David Clabaugh
//
// Per FORKERNELS_BUILD_STANDARD:
//   - Output: zig-out/{target}/lib/libforapollo.a
//   - Targets: linux-arm64, linux-x86_64, macos-arm64, windows-x86_64
//   - Sibling resolution: prebuilt/{target}/lib/ -> ../sibling/zig-out/{target}/lib/
//     -> ../sibling/zig-out/lib/ (fallback)
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
        // Standard paths
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/prebuilt/{s}/lib", .{ sibling, target_name }) });
        step.addLibraryPath(.{ .cwd_relative = b.fmt("{s}/zig-out/{s}/lib", .{ sibling, target_name }) });
        // Fallback for non-standardized repos
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
    const is_macos = step.rootModuleTarget().os.tag == .macos;
    const is_linux = step.rootModuleTarget().os.tag == .linux;
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
    }
    step.linkSystemLibrary("gfortran");
    step.linkSystemLibrary("gomp");
    step.linkLibC();
}

// ---------------------------------------------------------------------------
// Build entry point
// ---------------------------------------------------------------------------

pub fn build(b: *std.Build) void {
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
        const install = b.addInstallArtifact(static_lib, .{
            .dest_dir = .{ .override = .{ .custom = b.fmt("{s}/lib", .{target_name}) } },
        });
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

        const install = b.addInstallArtifact(shared_lib, .{
            .dest_dir = .{ .override = .{ .custom = b.fmt("{s}/lib", .{target_name}) } },
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
