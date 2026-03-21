// forApollo — Zig Build System (Stage 2)
// Copyright The Fantastic Planet — By David Clabaugh
//
// Stage 1: Makefile compiles Fortran → libforapollo_fortran.a
// Stage 2: This file links Fortran objects + deps → static + shared libraries
//
// Usage:
//   zig build                         # standard build
//   zig build -Duse-prebuilt=true     # use prebuilt/libforapollo_fortran.a
//   zig build -Dgenerate-prebuilt=true # copy output to prebuilt/
//   zig build -Ddev=true              # enable debug assertions and logging

const std = @import("std");

// ---------------------------------------------------------------------------
// Dependency archive lists
// ---------------------------------------------------------------------------

/// All formath component archives (order matters for static linking)
const formath_archives = [_][]const u8{
    "deps/formath/lib/libformath_linalg.a",
    "deps/formath/lib/libformath_ode.a",
    "deps/formath/lib/libformath_quaternion.a",
    "deps/formath/lib/libformath_liegroups.a",
    "deps/formath/lib/libformath_optimize.a",
    "deps/formath/lib/libformath_special.a",
    "deps/formath/lib/libformath_random.a",
    "deps/formath/lib/libformath_numdiff.a",
    "deps/formath/lib/libformath_fft.a",
};

/// Single-archive upstream dependencies
const upstream_archives = [_][]const u8{
    "deps/forfft/lib/libforfft.a",
    "deps/foropt/lib/libforopt.a",
    "deps/forternary/lib/libforternary.a",
    "deps/forgraph/lib/libforgraph.a",
};

// ---------------------------------------------------------------------------
// Helper: link all dependency archives + system libs onto a compile step
// ---------------------------------------------------------------------------

/// linkDeps attaches every Fortran archive and required system library to
/// `step`. Call this for both the static library, the shared library, and
/// the test runner so the logic lives in exactly one place.
fn linkDeps(b: *std.Build, step: *std.Build.Step.Compile, fortran_obj: []const u8) void {
    // Fortran kernel archive (stage-1 output)
    step.addObjectFile(b.path(fortran_obj));

    // formath component archives
    for (formath_archives) |archive| {
        step.addObjectFile(b.path(archive));
    }

    // Other upstream archives
    for (upstream_archives) |archive| {
        step.addObjectFile(b.path(archive));
    }

    // System runtime libraries
    step.linkSystemLibrary("gfortran"); // Fortran runtime
    step.linkSystemLibrary("gomp"); // OpenMP runtime
    step.linkLibC();
}

// ---------------------------------------------------------------------------
// Build entry point
// ---------------------------------------------------------------------------

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // -----------------------------------------------------------------------
    // Build options
    // -----------------------------------------------------------------------

    const use_prebuilt = b.option(
        bool,
        "use-prebuilt",
        "Link against prebuilt/libforapollo_fortran.a instead of build/",
    ) orelse false;

    const generate_prebuilt = b.option(
        bool,
        "generate-prebuilt",
        "Copy the Fortran archive to prebuilt/ after building",
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
    // Resolve Fortran archive path
    // -----------------------------------------------------------------------

    const fortran_archive: []const u8 = if (use_prebuilt)
        "prebuilt/libforapollo_fortran.a"
    else
        "build/libforapollo_fortran.a";

    // Optionally emit a warning so the developer notices the prebuilt flag
    if (generate_prebuilt) {
        std.debug.print(
            "[forapollo] -Dgenerate-prebuilt=true: after building, copy " ++
                "build/libforapollo_fortran.a → prebuilt/\n",
            .{},
        );
    }

    // -----------------------------------------------------------------------
    // Static library: libforapollo.a
    // -----------------------------------------------------------------------

    const static_lib = b.addStaticLibrary(.{
        .name = "forapollo",
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });
    static_lib.root_module.addOptions("build_opts", build_opts);
    linkDeps(b, static_lib, fortran_archive);
    b.installArtifact(static_lib);

    // -----------------------------------------------------------------------
    // Shared library: libforapollo.so / .dylib  (Python ctypes target)
    // -----------------------------------------------------------------------

    const shared_lib = b.addSharedLibrary(.{
        .name = "forapollo",
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
        .version = .{ .major = 0, .minor = 1, .patch = 0 },
    });
    shared_lib.root_module.addOptions("build_opts", build_opts);
    linkDeps(b, shared_lib, fortran_archive);
    b.installArtifact(shared_lib);

    // -----------------------------------------------------------------------
    // Test step
    // -----------------------------------------------------------------------

    const unit_tests = b.addTest(.{
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });
    unit_tests.root_module.addOptions("build_opts", build_opts);
    linkDeps(b, unit_tests, fortran_archive);

    const run_tests = b.addRunArtifact(unit_tests);
    const test_step = b.step("test", "Run Zig unit tests");
    test_step.dependOn(&run_tests.step);
}
