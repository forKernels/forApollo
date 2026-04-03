# forApollo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a universal state estimation, navigation, and guidance engine in Fortran 2008 + Zig + Python, following the Zortran pattern.

**Architecture:** 9 Fortran kernel files (engine/models/domain layers), 5 Zig safety/dispatch files, 10 Python ctypes binding files. Dependencies on forMath, forFFT, forOpt, forTernary, forGraph via prebuilt `.a` archives and Zig `@import`. All matrices flat row-major at C ABI boundary. All symbols `fa_*` prefix in Fortran, `forapollo_*` in Zig exports.

**Tech Stack:** Fortran 2008 (gfortran), Zig 0.15.x, Python 3.10, OpenMP, ctypes

**Spec:** `docs/superpowers/specs/2026-03-21-forapollo-universal-gnc-design.md`

**CRITICAL RULES:**
- No upstream `.f90` in this repo. Link prebuilt `.a` archives only.
- No Zig wrapper reimplementation. `@import("formath")` for upstream wrappers.
- No placeholder Fortran kernels. Every function has at least a real implementation with comments.
- No fake benchmarks. Only measured results.
- All matrices flat 1D row-major at `bind(C)` boundary.
- Copyright The Fantastic Planet — By David Clabaugh on all files.

---

## Phase Overview

This plan is organized into 6 phases. Each phase produces working, testable software:

| Phase | What | Fortran Files | Zig | Python | Tests |
|-------|------|--------------|-----|--------|-------|
| 1 | Scaffold + Build System | — | build.zig, stub Zig | — | build compiles |
| 2 | Models Layer | dynamics, observe | fortran.zig (externs) | — | Fortran unit tests |
| 3 | Engine Core | estimate, propagate | dispatch, kernels, safety | — | EKF end-to-end |
| 4 | Engine Extended | guidance, coords | exports.zig | — | guidance tests |
| 5 | Domain Layer | astro, environ, time | — | — | orbital tests |
| 6 | Python Bindings | — | — | all python/ | Python integration tests |

Phases 2-5 each follow TDD: write test → implement kernel → verify → commit.

---

## File Map

### Files to Create

```
forApollo/
├── src/
│   ├── fortran/
│   │   ├── forapollo_estimate.f90      (Phase 3) ~2000 lines — 16 estimators
│   │   ├── forapollo_propagate.f90     (Phase 3) ~400 lines — propagate, STM, batch
│   │   ├── forapollo_guidance.f90      (Phase 4) ~1500 lines — 23 guidance laws
│   │   ├── forapollo_coords.f90        (Phase 4) ~800 lines — frame transforms
│   │   ├── forapollo_dynamics.f90      (Phase 2) ~1200 lines — 17 dynamics models + Jacobians
│   │   ├── forapollo_observe.f90       (Phase 2) ~800 lines — 16 measurement models + Jacobians
│   │   ├── forapollo_astro.f90         (Phase 5) ~1000 lines — ephemeris, elements, transfers
│   │   ├── forapollo_environ.f90       (Phase 5) ~600 lines — atmosphere, gravity, geodesy
│   │   └── forapollo_time.f90          (Phase 5) ~400 lines — time systems
│   └── zig/
│       ├── fortran.zig                 (Phase 2) — extern fn declarations for all fa_*
│       ├── safety.zig                  (Phase 3) — bounds checking, error mapping
│       ├── dispatch.zig                (Phase 3) — size-based Zig/Fortran routing
│       ├── kernels.zig                 (Phase 3) — Zig-native small-state impls
│       └── exports.zig                 (Phase 4) — C ABI exports (forapollo_*)
├── python/
│   └── forapollo/
│       ├── __init__.py                 (Phase 6)
│       ├── _lib.py                     (Phase 6) — ctypes loader
│       ├── estimators.py               (Phase 6)
│       ├── dynamics.py                 (Phase 6)
│       ├── sensors.py                  (Phase 6)
│       ├── guidance.py                 (Phase 6)
│       ├── coords.py                   (Phase 6)
│       ├── astro.py                    (Phase 6)
│       ├── environ.py                  (Phase 6)
│       └── time.py                     (Phase 6)
├── tests/
│   ├── fortran/
│   │   ├── test_dynamics.f90           (Phase 2)
│   │   ├── test_observe.f90            (Phase 2)
│   │   ├── test_estimate.f90           (Phase 3)
│   │   ├── test_propagate.f90          (Phase 3)
│   │   ├── test_guidance.f90           (Phase 4)
│   │   ├── test_coords.f90             (Phase 4)
│   │   ├── test_astro.f90              (Phase 5)
│   │   ├── test_environ.f90            (Phase 5)
│   │   └── test_time.f90               (Phase 5)
│   └── python/
│       ├── test_estimators.py          (Phase 6)
│       ├── test_dynamics.py            (Phase 6)
│       └── test_integration.py         (Phase 6)
├── deps/                               (Phase 1) — prebuilt upstream archives
│   ├── formath/
│   │   ├── lib/                        ← .a archives
│   │   └── zig/                        ← Zig modules
│   ├── forfft/
│   ├── foropt/
│   ├── forternary/
│   └── forgraph/
├── build.zig                           (Phase 1)
├── Makefile                            (Phase 1)
├── setup.py                            (Phase 6)
└── requirements.txt                    (Phase 6)
```

### Files to Modify

```
CLAUDE.md        — already updated (done during brainstorming)
README.md        — already updated (done during brainstorming)
.gitignore       — already updated (done during brainstorming)
```

---

## Phase 1: Scaffold + Build System

### Task 1.1: Create Directory Structure

**Files:**
- Create: `src/fortran/` (directory)
- Create: `src/zig/` (directory)
- Create: `tests/fortran/` (directory)
- Create: `tests/python/` (directory)
- Create: `python/forapollo/` (directory)
- Create: `deps/formath/lib/` `deps/formath/zig/` (directories)
- Create: `deps/forfft/lib/` `deps/forfft/zig/` (directories)
- Create: `deps/foropt/lib/` `deps/foropt/zig/` (directories)
- Create: `deps/forternary/lib/` `deps/forternary/zig/` (directories)
- Create: `deps/forgraph/lib/` `deps/forgraph/zig/` (directories)
- Create: `prebuilt/aarch64-linux-gnu/` `prebuilt/aarch64-apple-darwin/` (directories)

- [ ] **Step 1: Create all directories**

```bash
mkdir -p src/fortran src/zig tests/fortran tests/python python/forapollo
mkdir -p deps/formath/{lib,zig} deps/forfft/{lib,zig} deps/foropt/{lib,zig}
mkdir -p deps/forternary/{lib,zig} deps/forgraph/{lib,zig}
mkdir -p prebuilt/{aarch64-linux-gnu,aarch64-apple-darwin}
```

- [ ] **Step 2: Verify structure**

```bash
find . -type d | grep -v .git | sort
```

Expected: all directories listed above exist.

- [ ] **Step 3: Add .gitkeep to empty dirs that need tracking**

```bash
touch deps/formath/lib/.gitkeep deps/formath/zig/.gitkeep
touch deps/forfft/lib/.gitkeep deps/forfft/zig/.gitkeep
touch deps/foropt/lib/.gitkeep deps/foropt/zig/.gitkeep
touch deps/forternary/lib/.gitkeep deps/forternary/zig/.gitkeep
touch deps/forgraph/lib/.gitkeep deps/forgraph/zig/.gitkeep
touch prebuilt/aarch64-linux-gnu/.gitkeep
touch prebuilt/aarch64-apple-darwin/.gitkeep
```

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "scaffold: create forApollo directory structure"
```

### Task 1.2: Create Makefile (Fortran compilation — Stage 1)

**Files:**
- Create: `Makefile`

The Makefile compiles Fortran sources to `.o` objects. Zig handles linking in Stage 2.

- [ ] **Step 1: Write Makefile**

```makefile
# forApollo — Fortran Compilation (Stage 1)
# Copyright The Fantastic Planet — By David Clabaugh
#
# Usage:
#   make all         — compile all Fortran kernels
#   make clean       — remove compiled objects
#   make prebuilt    — copy objects to prebuilt/ for current platform
#
# Stage 2 (Zig) links these objects + deps/*.a into final library.

FC = gfortran
MARCH ?= native
FFLAGS = -O3 -march=$(MARCH) -ftree-vectorize -fPIC -cpp \
         -std=f2008 -fall-intrinsics \
         -Wall -Wextra -Wimplicit-interface \
         -fopenmp

SRC_DIR = src/fortran
OBJ_DIR = build/obj

# Kernel sources (order matters for module dependencies)
SOURCES = \
    $(SRC_DIR)/forapollo_dynamics.f90 \
    $(SRC_DIR)/forapollo_observe.f90 \
    $(SRC_DIR)/forapollo_propagate.f90 \
    $(SRC_DIR)/forapollo_estimate.f90 \
    $(SRC_DIR)/forapollo_guidance.f90 \
    $(SRC_DIR)/forapollo_coords.f90 \
    $(SRC_DIR)/forapollo_astro.f90 \
    $(SRC_DIR)/forapollo_environ.f90 \
    $(SRC_DIR)/forapollo_time.f90

OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SOURCES))

# Platform detection for prebuilt target
UNAME_M := $(shell uname -m)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    PLATFORM = aarch64-apple-darwin
else
    PLATFORM = aarch64-linux-gnu
endif

.PHONY: all clean prebuilt

all: $(OBJ_DIR) $(OBJECTS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJ_DIR)

# Archive into static library
lib: all
	ar rcs build/libforapollo_fortran.a $(OBJECTS)

prebuilt: lib
	cp build/libforapollo_fortran.a prebuilt/$(PLATFORM)/libforapollo_fortran.a

clean:
	rm -rf build/

# Test compilation (links with deps)
test: all
	$(FC) $(FFLAGS) tests/fortran/test_dynamics.f90 $(OBJECTS) -o build/test_dynamics
	./build/test_dynamics
```

- [ ] **Step 2: Verify Makefile parses**

```bash
make -n all 2>&1 | head -5
```

Expected: prints gfortran commands (will fail since no .f90 files exist yet — that's fine).

- [ ] **Step 3: Commit**

```bash
git add Makefile && git commit -m "build: add Makefile for Fortran compilation (stage 1)"
```

### Task 1.3: Create build.zig (Zig linking — Stage 2)

**Files:**
- Create: `build.zig`
- Create: `build.zig.zon`

- [ ] **Step 1: Write build.zig.zon**

```zig
// Copyright The Fantastic Planet — By David Clabaugh
.{
    .name = .@"forapollo",
    .version = .{ 0, 1, 0 },
    .fingerprint = 0xfa_a9_01_10_00_00_00_00_00_00_00_00_00_00_00_00,
    .minimum_zig_version = .{ 0, 15, 0 },
    .paths = .{ "build.zig", "build.zig.zon", "src/" },
}
```

- [ ] **Step 2: Write build.zig**

```zig
// forApollo — Zig Build (Stage 2)
// Copyright The Fantastic Planet — By David Clabaugh
//
// Links prebuilt Fortran objects + upstream deps into final library.
// Usage:
//   zig build                          — default build
//   zig build -Duse-prebuilt=true      — use prebuilt Fortran objects
//   zig build -Dgenerate-prebuilt=true — regenerate prebuilt objects
//   zig build -Ddev=true               — source-build deps
//   zig build test                     — run Zig tests

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const use_prebuilt = b.option(bool, "use-prebuilt", "Use prebuilt Fortran objects") orelse false;
    const generate_prebuilt = b.option(bool, "generate-prebuilt", "Regenerate prebuilt objects") orelse false;
    const dev = b.option(bool, "dev", "Source-build deps (development only)") orelse false;
    _ = dev;
    _ = generate_prebuilt;

    // --- Static library ---
    const lib = b.addStaticLibrary(.{
        .name = "forapollo",
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Link prebuilt Fortran objects or compile from source
    if (use_prebuilt) {
        // Platform-specific prebuilt archive
        lib.addObjectFile(b.path("prebuilt/libforapollo_fortran.a"));
    } else {
        // Link Fortran objects compiled by Makefile (stage 1)
        lib.addObjectFile(b.path("build/libforapollo_fortran.a"));
    }

    // Link upstream dependency archives
    lib.addObjectFile(b.path("deps/formath/lib/libformath_linalg.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_ode.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_quaternion.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_liegroups.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_optimize.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_special.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_random.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_numdiff.a"));
    lib.addObjectFile(b.path("deps/formath/lib/libformath_fft.a"));
    lib.addObjectFile(b.path("deps/forfft/lib/libforfft.a"));
    lib.addObjectFile(b.path("deps/foropt/lib/libforopt.a"));
    lib.addObjectFile(b.path("deps/forternary/lib/libforternary.a"));
    lib.addObjectFile(b.path("deps/forgraph/lib/libforgraph.a"));

    // System libraries
    lib.linkSystemLibrary("gfortran");
    lib.linkSystemLibrary("gomp"); // OpenMP
    lib.linkLibC();

    b.installArtifact(lib);

    // --- Shared library (for Python ctypes) ---
    const shared = b.addSharedLibrary(.{
        .name = "forapollo",
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Same linking as static
    if (use_prebuilt) {
        shared.addObjectFile(b.path("prebuilt/libforapollo_fortran.a"));
    } else {
        shared.addObjectFile(b.path("build/libforapollo_fortran.a"));
    }
    shared.addObjectFile(b.path("deps/formath/lib/libformath_linalg.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_ode.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_quaternion.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_liegroups.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_optimize.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_special.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_random.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_numdiff.a"));
    shared.addObjectFile(b.path("deps/formath/lib/libformath_fft.a"));
    shared.addObjectFile(b.path("deps/forfft/lib/libforfft.a"));
    shared.addObjectFile(b.path("deps/foropt/lib/libforopt.a"));
    shared.addObjectFile(b.path("deps/forternary/lib/libforternary.a"));
    shared.addObjectFile(b.path("deps/forgraph/lib/libforgraph.a"));
    shared.linkSystemLibrary("gfortran");
    shared.linkSystemLibrary("gomp");
    shared.linkLibC();

    b.installArtifact(shared);

    // --- Tests ---
    const tests = b.addTest(.{
        .root_source_file = b.path("src/zig/exports.zig"),
        .target = target,
        .optimize = optimize,
    });
    const run_tests = b.addRunArtifact(tests);
    const test_step = b.step("test", "Run Zig tests");
    test_step.dependOn(&run_tests.step);
}
```

- [ ] **Step 3: Create stub exports.zig so build.zig parses**

```zig
// forApollo — Zig Exports (stub)
// Copyright The Fantastic Planet — By David Clabaugh
//
// C ABI exports for Python/Rust/C++ consumers.
// All functions prefixed forapollo_* and call fa_* Fortran kernels
// through the safety layer.

pub export fn forapollo_version() callconv(.c) u32 {
    return 0x000100; // 0.1.0
}
```

Write to `src/zig/exports.zig`.

- [ ] **Step 4: Verify zig build parses**

```bash
zig build --help 2>&1 | head -5
```

Expected: shows build options (will fail to link since deps don't exist yet — that's expected).

- [ ] **Step 5: Commit**

```bash
git add build.zig build.zig.zon src/zig/exports.zig && git commit -m "build: add Zig build system (stage 2) with dep linking"
```

### Task 1.4: Populate Dependency Stubs

**Files:**
- Create: `deps/formath/zig/formath.zig` (minimal extern stubs — only symbols forApollo calls)
- Create: `deps/forfft/zig/forfft.zig`
- Create: `deps/foropt/zig/foropt.zig`
- Create: `deps/forternary/zig/forternary.zig`
- Create: `deps/forgraph/zig/forgraph.zig`

These are placeholder Zig extern declarations. They'll be replaced with real prebuilt Zig modules when we copy them from the upstream repos.

- [ ] **Step 1: Write formath.zig stub**

```zig
// forMath — Minimal extern declarations for symbols forApollo uses.
// Copyright The Fantastic Planet — By David Clabaugh
// This file is auto-generated from forMath's prebuilt Zig modules.
// Only the subset of fm_* symbols that forApollo calls are declared here.

// Linear algebra
pub extern "c" fn fm_dgemm(transa: u8, transb: u8, m: c_int, n: c_int, k: c_int, alpha: f64, a: [*]const f64, lda: c_int, b: [*]const f64, ldb: c_int, beta: f64, c_ptr: [*]f64, ldc: c_int) void;
pub extern "c" fn fm_dgemv(trans: u8, m: c_int, n: c_int, alpha: f64, a: [*]const f64, lda: c_int, x: [*]const f64, incx: c_int, beta: f64, y: [*]f64, incy: c_int) void;

// ODE solvers
pub extern "c" fn fm_ode_rk4_step(n: c_int, x: [*]f64, t: f64, dt: f64, f_ptr: ?*const anyopaque, params: [*]const f64, np: c_int, info: *c_int) void;

// Quaternion
pub extern "c" fn fm_quat_multiply(a: [*]const f64, b: [*]const f64, c_ptr: [*]f64) void;
pub extern "c" fn fm_quat_to_dcm(q: [*]const f64, dcm: [*]f64) void;

// Numerical differentiation
pub extern "c" fn fm_numdiff_jacobian(n: c_int, m: c_int, x: [*]const f64, f_ptr: ?*const anyopaque, J: [*]f64, info: *c_int) void;

// Random
pub extern "c" fn fm_random_normal(n: c_int, mean: f64, std: f64, out: [*]f64) void;

const c_int = i32;
```

NOTE: This is a STUB. Real prebuilt Zig modules from forMath will replace this. The stub lets us compile and test the Zig layer before deps are fully wired.

- [ ] **Step 2: Write similar stubs for forfft, foropt, forternary, forgraph**

Each with only the extern fns forApollo actually calls. Keep minimal.

- [ ] **Step 3: Commit**

```bash
git add deps/ && git commit -m "deps: add Zig extern stubs for upstream dependencies"
```

---

## Phase 2: Models Layer (Dynamics + Measurements)

Models come first because the engine (estimators) dispatches TO them. Build bottom-up.

### Task 2.1: Implement forapollo_dynamics.f90 — Tracking Models

Start with the simplest dynamics: constant-velocity, constant-acceleration, constant-turn. These are trivially testable and exercise the full dispatch pattern.

**Files:**
- Create: `src/fortran/forapollo_dynamics.f90`
- Create: `tests/fortran/test_dynamics.f90`

- [ ] **Step 1: Write test for constant-velocity dynamics**

In `tests/fortran/test_dynamics.f90`:

```fortran
! test_dynamics.f90 — Unit tests for forapollo_dynamics
! Copyright The Fantastic Planet — By David Clabaugh
program test_dynamics
    use iso_c_binding
    implicit none

    ! CRITICAL: Interface block required for bind(C) subroutines with value args.
    ! Without this, gfortran passes integers by reference (Fortran default) instead
    ! of by value, causing segfaults or wrong results.
    interface
        subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
            bind(C, name="fa_dynamics_dispatch")
            use iso_c_binding
            integer(c_int), value, intent(in) :: model_id, n, nu, np
            real(c_double), intent(in) :: x(n), u(nu), t, params(np)
            real(c_double), intent(out) :: x_dot(n)
            integer(c_int), intent(out) :: info
        end subroutine
        subroutine fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info) &
            bind(C, name="fa_dynamics_jacobian")
            use iso_c_binding
            integer(c_int), value, intent(in) :: model_id, n, nu, np
            real(c_double), intent(in) :: x(n), u(nu), t, params(np)
            real(c_double), intent(out) :: F(n*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    integer(c_int) :: info
    real(c_double) :: x(4), x_dot(4), params(1), u(1)
    real(c_double) :: F(16), F_fd(16)  ! Jacobian + finite-diff check
    integer(c_int), parameter :: FA_DYN_CONST_VEL = 40
    logical :: pass

    pass = .true.

    ! Test 1: constant-velocity in 2D
    ! State: [x, y, vx, vy] = [0, 0, 1, 2]
    ! Expected x_dot: [vx, vy, 0, 0] = [1, 2, 0, 0]
    x = [0.0d0, 0.0d0, 1.0d0, 2.0d0]
    u = [0.0d0]
    params = [0.0d0]

    call fa_dynamics_dispatch(FA_DYN_CONST_VEL, 4, x, u, 0, 0.0d0, params, 0, x_dot, info)

    if (info /= 0) then
        print *, "FAIL: const_vel info =", info
        pass = .false.
    end if
    if (abs(x_dot(1) - 1.0d0) > 1.0d-12) then
        print *, "FAIL: x_dot(1) =", x_dot(1), "expected 1.0"
        pass = .false.
    end if
    if (abs(x_dot(2) - 2.0d0) > 1.0d-12) then
        print *, "FAIL: x_dot(2) =", x_dot(2), "expected 2.0"
        pass = .false.
    end if
    if (abs(x_dot(3)) > 1.0d-12 .or. abs(x_dot(4)) > 1.0d-12) then
        print *, "FAIL: acceleration should be zero"
        pass = .false.
    end if

    ! Test 2: Jacobian verification — analytic vs finite-difference
    ! For const-vel, F = [0 I; 0 0], constant. Verify against finite-diff.
    call fa_dynamics_jacobian(FA_DYN_CONST_VEL, 4, x, u, 0, 0.0d0, params, 0, F, info)
    if (info /= 0) then
        print *, "FAIL: const_vel Jacobian info =", info
        pass = .false.
    end if
    ! F(1,3)=1, F(2,4)=1, rest zero (row-major flat)
    if (abs(F(3) - 1.0d0) > 1.0d-12) then  ! row 1, col 3
        print *, "FAIL: F(1,3) =", F(3), "expected 1.0"
        pass = .false.
    end if

    ! Test 3: Edge case — invalid model_id should return info=3
    call fa_dynamics_dispatch(999, 4, x, u, 0, 0.0d0, params, 0, x_dot, info)
    if (info /= 3) then
        print *, "FAIL: invalid model_id should return info=3, got", info
        pass = .false.
    end if

    ! Test 4: Edge case — odd n for const_vel should return info=3
    call fa_dynamics_dispatch(FA_DYN_CONST_VEL, 3, x(1:3), u, 0, 0.0d0, params, 0, x_dot(1:3), info)
    if (info /= 3) then
        print *, "FAIL: odd n should return info=3, got", info
        pass = .false.
    end if

    if (pass) then
        print *, "PASS: all dynamics tests passed"
    else
        error stop "FAIL: dynamics tests failed"
    end if
end program test_dynamics
```

**NOTE: Every dynamics and observation model MUST have a Jacobian-vs-finite-difference verification test.** Use forMath's `fm_numdiff_jacobian` (once deps are wired) or a simple central-difference loop:
```fortran
! Central difference: dF/dx_j ≈ (f(x+h*e_j) - f(x-h*e_j)) / (2*h)
do j = 1, n
    x_plus = x; x_plus(j) = x_plus(j) + h
    x_minus = x; x_minus(j) = x_minus(j) - h
    call fa_dynamics_dispatch(model_id, n, x_plus, ..., f_plus, info)
    call fa_dynamics_dispatch(model_id, n, x_minus, ..., f_minus, info)
    F_fd((i-1)*n + j) = (f_plus(i) - f_minus(i)) / (2.0d0 * h)
end do
! Then: assert max|F_analytic - F_fd| < 1e-6
```

- [ ] **Step 2: Verify test doesn't compile (no implementation yet)**

```bash
gfortran -c tests/fortran/test_dynamics.f90 -o /dev/null 2>&1
```

Expected: error about missing `fa_dynamics_dispatch`.

- [ ] **Step 3: Implement forapollo_dynamics.f90 — dispatch + tracking models**

In `src/fortran/forapollo_dynamics.f90`:

```fortran
! forapollo_dynamics.f90 — Built-in dynamics model catalog with dispatch
! Copyright The Fantastic Planet — By David Clabaugh
!
! Each model computes x_dot = f(x, u, t; params) for a specific system.
! Models are selected by integer ID. Analytic Jacobians provided for each.
!
! Architecture: Fortran kernel (fa_*) → Zig safety → Python ctypes
! All matrices flat 1D row-major at bind(C) boundary.

! --- Model ID constants ---
! Orbital (1-9)
! integer(c_int), parameter :: FA_DYN_KEPLER = 1          (Task 2.2)
! Tracking (40-49)
! integer(c_int), parameter :: FA_DYN_CONST_VEL = 40
! integer(c_int), parameter :: FA_DYN_CONST_ACCEL = 41
! integer(c_int), parameter :: FA_DYN_CONST_TURN = 42

! ==========================================================================
! DISPATCH — routes model_id to specific dynamics subroutine
! ==========================================================================
subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
    bind(C, name="fa_dynamics_dispatch")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: model_id, n, nu, np
    real(c_double), intent(in) :: x(n), u(nu), t, params(np)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info

    info = 0
    x_dot = 0.0d0

    select case (model_id)
        ! --- Tracking models (40-49) ---
        case (40)  ! FA_DYN_CONST_VEL
            call fa_dyn_const_vel(n, x, x_dot, info)
        case (41)  ! FA_DYN_CONST_ACCEL
            call fa_dyn_const_accel(n, x, x_dot, info)
        case (42)  ! FA_DYN_CONST_TURN
            call fa_dyn_const_turn(n, x, x_dot, info)
        case default
            info = 3  ! invalid model_id
    end select
end subroutine fa_dynamics_dispatch

! ==========================================================================
! JACOBIAN DISPATCH — analytic df/dx for built-in models
! ==========================================================================
subroutine fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info) &
    bind(C, name="fa_dynamics_jacobian")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: model_id, n, nu, np
    real(c_double), intent(in) :: x(n), u(nu), t, params(np)
    real(c_double), intent(out) :: F(n*n)
    integer(c_int), intent(out) :: info

    info = 0
    F = 0.0d0

    select case (model_id)
        case (40)
            call fa_dyn_const_vel_jac(n, F, info)
        case (41)
            call fa_dyn_const_accel_jac(n, F, info)
        case (42)
            call fa_dyn_const_turn_jac(n, x, F, info)
        case default
            info = 3
    end select
end subroutine fa_dynamics_jacobian

! ==========================================================================
! CONSTANT VELOCITY: x_dot = [v; 0]
! State: [pos(d), vel(d)], n = 2*d, d inferred from n/2
! ==========================================================================
! NOTE: Internal subroutines do NOT have bind(C), so they MUST NOT use
! the 'value' attribute (only valid with bind(C) in Fortran 2008).
! Pass n by reference (Fortran default).
subroutine fa_dyn_const_vel(n, x, x_dot, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info
    integer :: d

    info = 0
    d = n / 2
    if (2*d /= n) then
        info = 3  ! n must be even
        return
    end if

    ! Position derivative = velocity
    x_dot(1:d) = x(d+1:n)
    ! Velocity derivative = zero (constant velocity assumption)
    x_dot(d+1:n) = 0.0d0
end subroutine fa_dyn_const_vel

subroutine fa_dyn_const_vel_jac(n, F, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(out) :: F(n*n)
    integer(c_int), intent(out) :: info
    integer :: d, i

    info = 0
    d = n / 2
    F = 0.0d0

    ! F = [0 I; 0 0] in row-major flat layout
    ! Row i (0-indexed), col d+i = 1.0 for i = 0..d-1
    do i = 1, d
        F((i-1)*n + d + i) = 1.0d0
    end do
end subroutine fa_dyn_const_vel_jac

! ==========================================================================
! CONSTANT ACCELERATION: x_dot = [v; a; 0]
! State: [pos(d), vel(d), acc(d)], n = 3*d, d inferred from n/3
! ==========================================================================
subroutine fa_dyn_const_accel(n, x, x_dot, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info
    integer :: d

    info = 0
    d = n / 3
    if (3*d /= n) then
        info = 3  ! n must be divisible by 3
        return
    end if

    ! pos_dot = vel
    x_dot(1:d) = x(d+1:2*d)
    ! vel_dot = acc
    x_dot(d+1:2*d) = x(2*d+1:n)
    ! acc_dot = 0 (constant acceleration)
    x_dot(2*d+1:n) = 0.0d0
end subroutine fa_dyn_const_accel

subroutine fa_dyn_const_accel_jac(n, F, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(out) :: F(n*n)
    integer(c_int), intent(out) :: info
    integer :: d, i

    info = 0
    d = n / 3
    F = 0.0d0

    ! F = [0 I 0; 0 0 I; 0 0 0] in row-major flat
    ! Row i, col d+i = 1 for i=0..d-1 (pos→vel)
    do i = 1, d
        F((i-1)*n + d + i) = 1.0d0
    end do
    ! Row d+i, col 2d+i = 1 for i=0..d-1 (vel→acc)
    do i = 1, d
        F((d+i-1)*n + 2*d + i) = 1.0d0
    end do
end subroutine fa_dyn_const_accel_jac

! ==========================================================================
! CONSTANT TURN RATE: 2D coordinated turn
! State: [x, y, vx, vy, omega], n = 5
! x_dot = [vx, vy, -omega*vy, omega*vx, 0]
! ==========================================================================
subroutine fa_dyn_const_turn(n, x, x_dot, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info
    real(c_double) :: vx, vy, omega

    info = 0
    if (n /= 5) then
        info = 3  ! constant turn requires n=5
        return
    end if

    vx = x(3)
    vy = x(4)
    omega = x(5)

    x_dot(1) = vx              ! dx/dt = vx
    x_dot(2) = vy              ! dy/dt = vy
    x_dot(3) = -omega * vy     ! dvx/dt = -omega * vy
    x_dot(4) = omega * vx      ! dvy/dt = omega * vx
    x_dot(5) = 0.0d0           ! domega/dt = 0 (constant turn rate)
end subroutine fa_dyn_const_turn

subroutine fa_dyn_const_turn_jac(n, x, F, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: F(n*n)
    integer(c_int), intent(out) :: info
    real(c_double) :: vx, vy, omega

    info = 0
    F = 0.0d0
    vx = x(3); vy = x(4); omega = x(5)

    ! Row-major flat: F(row, col) at index row*n + col (1-based: (row-1)*n + col)
    ! Row 1: dx/dt = vx → dF/dvx = 1
    F(3) = 1.0d0                      ! (1,3)
    ! Row 2: dy/dt = vy → dF/dvy = 1
    F(5+4) = 1.0d0                    ! (2,4)
    ! Row 3: dvx/dt = -omega*vy → dF/dvy=-omega, dF/domega=-vy
    F(2*5+4) = -omega                 ! (3,4)
    F(2*5+5) = -vy                    ! (3,5)
    ! Row 4: dvy/dt = omega*vx → dF/dvx=omega, dF/domega=vx
    F(3*5+3) = omega                  ! (4,3)
    F(3*5+5) = vx                     ! (4,5)
    ! Row 5: all zeros (constant omega)
end subroutine fa_dyn_const_turn_jac
```

- [ ] **Step 4: Compile and run test**

```bash
gfortran -O3 -std=f2008 -fall-intrinsics -c src/fortran/forapollo_dynamics.f90 -o build/obj/forapollo_dynamics.o -Jbuild/obj
gfortran -O3 -std=f2008 -fall-intrinsics tests/fortran/test_dynamics.f90 build/obj/forapollo_dynamics.o -o build/test_dynamics
./build/test_dynamics
```

Expected: `PASS: all dynamics tests passed`

- [ ] **Step 5: Commit**

```bash
git add src/fortran/forapollo_dynamics.f90 tests/fortran/test_dynamics.f90
git commit -m "feat: add tracking dynamics models (const-vel, const-accel, const-turn) with Jacobians"
```

### Task 2.2: Add Orbital Dynamics (Kepler, J2, CR3BP, Drag)

**Files:**
- Modify: `src/fortran/forapollo_dynamics.f90` — add orbital models to dispatch + implementations
- Modify: `tests/fortran/test_dynamics.f90` — add orbital tests

- [ ] **Step 1: Add Kepler test to test_dynamics.f90**

Test: circular orbit at 400km altitude. State: [r, 0, 0, 0, v_circ, 0]. After propagating x_dot, radial acceleration should be -mu/r^2, tangential velocity derivative should be zero (circular orbit → purely radial gravity).

```fortran
! Test Kepler: circular orbit, x_dot should give gravitational acceleration
x6 = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0]
params_k(1) = 3.986004418d14  ! mu_earth
call fa_dynamics_dispatch(1, 6, x6, u, 0, 0.0d0, params_k, 1, x_dot6, info)
! x_dot(1:3) = v = [0, 7669, 0]
! x_dot(4:6) = -mu/r^3 * r = [-mu/r^2, 0, 0]
expected_accel = -3.986004418d14 / (6778.0d3)**2
if (abs(x_dot6(4) - expected_accel) / abs(expected_accel) > 1.0d-6) then
    print *, "FAIL: Kepler accel"
    pass = .false.
end if
```

- [ ] **Step 2: Implement Kepler dynamics**

Add to `forapollo_dynamics.f90`:

```fortran
! ==========================================================================
! KEPLER: Two-body gravitational dynamics
! State: [rx, ry, rz, vx, vy, vz], n = 6
! Params: [mu] (gravitational parameter, m^3/s^2)
! x_dot = [v; -mu/|r|^3 * r]
! ==========================================================================
subroutine fa_dyn_kepler(n, x, params, np, x_dot, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n, np
    real(c_double), intent(in) :: x(n), params(np)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info
    real(c_double) :: mu, r_mag, r_mag3

    info = 0
    if (n /= 6 .or. np < 1) then
        info = 3
        return
    end if

    mu = params(1)
    r_mag = sqrt(x(1)**2 + x(2)**2 + x(3)**2)

    if (r_mag < 1.0d-10) then
        info = 1  ! degenerate — collision
        return
    end if

    r_mag3 = r_mag**3

    ! Velocity derivatives = velocity components
    x_dot(1) = x(4)
    x_dot(2) = x(5)
    x_dot(3) = x(6)

    ! Acceleration = -mu/r^3 * r (gravitational)
    x_dot(4) = -mu * x(1) / r_mag3
    x_dot(5) = -mu * x(2) / r_mag3
    x_dot(6) = -mu * x(3) / r_mag3
end subroutine fa_dyn_kepler
```

Add `case (1)` to dispatch. Add Kepler Jacobian. Then implement J2, CR3BP, drag similarly.

- [ ] **Step 3: Add J2 perturbation dynamics**

J2 adds the oblateness perturbation to Kepler:
```
a_J2 = -3/2 * J2 * mu * R_eq^2 / r^5 * [x*(1 - 5*z^2/r^2), y*(1 - 5*z^2/r^2), z*(3 - 5*z^2/r^2)]
```
Params: [mu, J2, R_eq]

- [ ] **Step 4: Add CR3BP dynamics**

Circular restricted three-body problem in rotating frame:
```
x_dot = [vx, vy, vz, 2*vy + dOmega/dx, -2*vx + dOmega/dy, dOmega/dz]
```
where Omega is the pseudo-potential. Params: [mu_ratio] (mass ratio m2/(m1+m2))

- [ ] **Step 5: Add atmospheric drag dynamics**

Kepler + drag: `a_drag = -0.5 * rho * Cd*A/m * |v_rel| * v_rel`
State: [r(3), v(3), beta] where beta = Cd*A/m (ballistic coefficient).
Params: [mu, rho0, h_scale] for exponential atmosphere.

- [ ] **Step 6: Run all tests, commit**

```bash
make all && gfortran tests/fortran/test_dynamics.f90 build/obj/*.o -o build/test_dynamics && ./build/test_dynamics
git add src/fortran/forapollo_dynamics.f90 tests/fortran/test_dynamics.f90
git commit -m "feat: add orbital dynamics (Kepler, J2, CR3BP, drag) with Jacobians"
```

### Task 2.3a: Add Ground Vehicle Dynamics (Bicycle, Ackermann, Diff Drive)

**Files:**
- Modify: `src/fortran/forapollo_dynamics.f90`
- Modify: `tests/fortran/test_dynamics.f90`

- [ ] **Step 1: Write bicycle model test**

Straight line at 10 m/s, zero steer. x_dot = [10*cos(0), 10*sin(0), 0, 0] = [10, 0, 0, 0].

- [ ] **Step 2: Implement bicycle, Ackermann, differential drive (cases 20, 21, 22)**

Bicycle: `x_dot = [v*cos(theta), v*sin(theta), v*tan(steer)/L, accel]`
Ackermann: steer is state, steer_rate is control.
Diff drive: `v = (v_l + v_r)/2, omega = (v_r - v_l)/sep`
Each with analytic Jacobian.

- [ ] **Step 3: Verify analytic Jacobians against finite-difference for all 3 models**
- [ ] **Step 4: Run tests, commit**

```bash
git commit -m "feat: add ground vehicle dynamics (bicycle, Ackermann, diff drive)"
```

### Task 2.3b: Add Rigid Body 6-DOF Dynamics

**Files:**
- Modify: `src/fortran/forapollo_dynamics.f90`
- Modify: `tests/fortran/test_dynamics.f90`

- [ ] **Step 1: Write rigid body test — torque-free spinning body**

13-state: [r(3), v(3), q(4), omega(3)]. No forces/torques. Quaternion norm should be preserved. Angular momentum should be conserved.

- [ ] **Step 2: Implement rigid body 6-DOF (case 10)**

Quaternion kinematics: `q_dot = 0.5 * q ⊗ [0, omega]`
Euler's rotation equations: `I*omega_dot = T - omega × (I*omega)`
Params: [mass, Ixx, Iyy, Izz]. Control: [Fx,Fy,Fz, Tx,Ty,Tz].

- [ ] **Step 3: Verify Jacobian (13x13) against finite-difference**
- [ ] **Step 4: Run tests, commit**

```bash
git commit -m "feat: add rigid body 6-DOF dynamics with quaternion kinematics"
```

### Task 2.3c: Add Aerial Vehicle Dynamics (Quadrotor, Fixed-Wing)

**Files:**
- Modify: `src/fortran/forapollo_dynamics.f90`
- Modify: `tests/fortran/test_dynamics.f90`

- [ ] **Step 1: Write quadrotor hover test**

12-state, hover at constant altitude. With thrust = mg, all derivatives except gravity cancellation should be zero.

- [ ] **Step 2: Implement quadrotor (case 30) and fixed-wing (case 31)**

Control: [thrust, tau_x, tau_y, tau_z] (nu=4).
Builds on rigid body rotation equations but uses Euler angles.

- [ ] **Step 3: Verify Jacobians against finite-difference**
- [ ] **Step 4: Run tests, commit**

```bash
git commit -m "feat: add aerial dynamics (quadrotor, fixed-wing)"
```

### Task 2.3d: Add Stochastic + Scalar Dynamics

**Files:**
- Modify: `src/fortran/forapollo_dynamics.f90`
- Modify: `tests/fortran/test_dynamics.f90`

- [ ] **Step 1: Write GBM and O-U tests**

GBM: `x_dot = mu*S` (deterministic drift). O-U: `x_dot = theta*(mu-X)` (mean reversion).

- [ ] **Step 2: Implement GBM (50), O-U (51), double integrator (60), spring-mass (61)**

All trivial — 5-10 lines each. Each with analytic Jacobian.

- [ ] **Step 3: Verify Jacobians, run all tests, commit**

```bash
git commit -m "feat: add stochastic and scalar dynamics (GBM, O-U, spring-mass)"
```

### Task 2.4: Implement forapollo_observe.f90 — Measurement Models

**Files:**
- Create: `src/fortran/forapollo_observe.f90`
- Create: `tests/fortran/test_observe.f90`

- [ ] **Step 1: Write test for position and range observations**

Position: h(x) = x(1:3) → trivial.
Range: h(x) = ||x(1:3) - station|| → single scalar.

- [ ] **Step 2: Implement dispatch + position, range, bearing, range+bearing**

Same dispatch pattern as dynamics. Each measurement model has an analytic Jacobian.

Range Jacobian: `H = (x - station)^T / ||x - station||`
Bearing: `H = [-sin(az), cos(az)] / range` (2D)

- [ ] **Step 3: Implement velocity, Doppler, attitude sensors**

Doppler: range-rate = `(x-s)·(v-vs) / ||x-s||`
Magnetometer/star tracker/sun sensor: body-frame observation of known inertial vector via quaternion rotation.

- [ ] **Step 4: Implement remaining models (accel, gyro, radar, scalar, pinhole, relative)**

Each with analytic Jacobian.

- [ ] **Step 5: Run tests, commit**

```bash
make all && gfortran tests/fortran/test_observe.f90 build/obj/*.o -o build/test_observe && ./build/test_observe
git add src/fortran/forapollo_observe.f90 tests/fortran/test_observe.f90
git commit -m "feat: add 16 measurement models with analytic Jacobians"
```

---

## Phase 3: Engine Core (Estimators + Propagation)

### Task 3.1: Implement forapollo_propagate.f90

**Files:**
- Create: `src/fortran/forapollo_propagate.f90`
- Create: `tests/fortran/test_propagate.f90`

- [ ] **Step 1: Write test — propagate constant-velocity state**

Start: [0, 0, 1, 0], dt=10. Expected: [10, 0, 1, 0].

- [ ] **Step 2: Implement fa_propagate using RK4**

Simple fixed-step RK4 integration. Calls dynamics dispatch internally. Later can dispatch to forMath ODE solvers for adaptive stepping.

```fortran
subroutine fa_propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, integrator_id, info) &
    bind(C, name="fa_propagate")
    use iso_c_binding
    integer(c_int), value :: n, nu, model_id, np, integrator_id
    real(c_double), intent(inout) :: x(n)
    real(c_double), intent(in) :: u(nu), params(np), dt
    type(c_funptr), value :: f_ptr   ! null → use model_id
    integer(c_int), intent(out) :: info
    ! integrator_id: 0=RK4(internal), 1=DP5(forMath), 2=DP8(forMath), 3=DOP853(forMath)
    ! When integrator_id > 0, dispatches to forMath ODE solvers via @import("formath")
end subroutine
```

**IMPORTANT:** `integrator_id=0` uses an internal RK4 (no forMath dependency). `integrator_id > 0` dispatches to forMath's adaptive solvers through Zig. This lets Phase 3 work with internal RK4 before deps are fully wired, then users get adaptive stepping by passing a different `integrator_id`.

- [ ] **Step 3: Implement fa_propagate_stm (state + state transition matrix)**

Propagate [x(n), Phi(n*n)] together. Phi_dot = F * Phi where F is the dynamics Jacobian.

- [ ] **Step 4: Implement fa_propagate_batch (OpenMP)**

```fortran
subroutine fa_propagate_batch(n, N_states, x_batch, ...) bind(C, name="fa_propagate_batch")
    !$omp parallel do private(i)
    do i = 1, N_states
        call fa_propagate(n, x_batch((i-1)*n+1:i*n), ...)
    end do
end subroutine
```

- [ ] **Step 5: Run tests, commit**

```bash
git commit -m "feat: add state propagation (RK4, STM, OpenMP batch)"
```

### Task 3.2: Implement forapollo_estimate.f90 — KF and EKF

**Files:**
- Create: `src/fortran/forapollo_estimate.f90`
- Create: `tests/fortran/test_estimate.f90`

Start with the two most fundamental estimators: linear KF and EKF.

- [ ] **Step 1: Write KF test — constant velocity tracking with position observations**

Known analytical solution: predict with constant-velocity model, update with position measurement. After several iterations, covariance should converge to steady state.

- [ ] **Step 2: Implement fa_kf_predict and fa_kf_update**

```fortran
! KF predict: x = F*x, P = F*P*F^T + Q
! KF update: K = P*H^T*(H*P*H^T + R)^-1, x = x + K*(z - H*x), P = (I-K*H)*P
```

Uses forMath BLAS for matrix operations via Zig @import.

- [ ] **Step 3: Write EKF test — Kepler orbit with range measurements**

Propagate a circular orbit, observe range from a ground station. EKF should converge.

- [ ] **Step 4: Implement fa_ekf_predict and fa_ekf_update**

EKF predict: propagate x with fa_propagate, compute STM for P propagation.
EKF update: evaluate h(x), compute H Jacobian, standard Kalman update.
Includes ternary gating on innovation: chi-squared test → forTernary gate.

- [ ] **Step 5: Run tests, commit**

```bash
git commit -m "feat: add linear KF and EKF with ternary measurement gating"
```

### Task 3.3: Add IEKF, UKF, ESKF, Square-Root Variants

**Files:**
- Modify: `src/fortran/forapollo_estimate.f90`
- Modify: `tests/fortran/test_estimate.f90`

- [ ] **Step 1: Implement IEKF (iterated EKF)**

Re-linearizes at updated state estimate up to `max_iter` times. Converges tighter than standard EKF for highly nonlinear measurement models.

- [ ] **Step 2: Implement UKF (sigma-point)**

Generate 2n+1 sigma points, propagate each through dynamics, reconstruct mean and covariance. No Jacobians needed.

- [ ] **Step 3: Implement ESKF (error-state — Apollo heritage)**

Predict: propagate nominal state, propagate error-state covariance.
Update: update error state, inject correction into nominal state.
Three-step: predict → update → inject.

- [ ] **Step 4: Implement SR-EKF and SR-UKF**

Work with Cholesky factor S where P = S*S^T. QR-based updates to maintain triangular form.

- [ ] **Step 5: Add ternary gating test**

Verify: large-innovation measurement (10x Mahalanobis distance) returns validity=-1, normal measurement returns validity=+1, borderline returns validity=0.

- [ ] **Step 6: Add NEES (Normalized Estimation Error Squared) consistency test**

Run Monte Carlo: 50 runs of EKF on known trajectory. Compute NEES = (x_true-x_est)^T * P^-1 * (x_true-x_est). Average NEES should be near n (state dimension) if covariance is consistent.

- [ ] **Step 7: Tests and commit**

```bash
git commit -m "feat: add IEKF, UKF, ESKF, SR-EKF, SR-UKF with ternary gating and NEES tests"
```

### Task 3.4a: Add Information Filters (IF, EIF)

**Files:**
- Modify: `src/fortran/forapollo_estimate.f90`
- Modify: `tests/fortran/test_estimate.f90`

- [ ] **Step 1: Implement IF predict and update**

Information form: Y = P^-1, y = P^-1 * x. Updates are additive: Y_new = Y + H^T*R^-1*H. Natural for multi-sensor fusion and decentralized estimation.

- [ ] **Step 2: Implement EIF (extended information filter)**

Nonlinear variant — linearize like EKF but in information form.

- [ ] **Step 3: Test: multi-sensor fusion — verify IF gives same result as KF**
- [ ] **Step 4: Commit**

```bash
git commit -m "feat: add information filters (IF, EIF)"
```

### Task 3.4b: Add SIR Particle Filter with OpenMP

**Files:**
- Modify: `src/fortran/forapollo_estimate.f90`
- Modify: `tests/fortran/test_estimate.f90`

- [ ] **Step 1: Implement fa_pf_sir_predict**

Propagate each particle through dynamics + process noise. Uses forMath `fm_random_normal` for noise sampling.

- [ ] **Step 2: Implement fa_pf_sir_update**

Weight by likelihood: `w_i = p(z | x_i)`. OpenMP parallel for weight computation.

- [ ] **Step 3: Implement fa_pf_sir_resample (systematic resampling)**

Inherently sequential. Generate one uniform U(0,1/N), stride by 1/N.

- [ ] **Step 4: Test: bimodal distribution — PF should handle, EKF should fail**
- [ ] **Step 5: Commit**

```bash
git commit -m "feat: add SIR particle filter with OpenMP parallel weights"
```

### Task 3.4c: Add Auxiliary PF and Rao-Blackwellized PF

**Files:**
- Modify: `src/fortran/forapollo_estimate.f90`
- Modify: `tests/fortran/test_estimate.f90`

- [ ] **Step 1: Implement auxiliary particle filter**

First-stage weights using predicted likelihood for better proposal distribution.

- [ ] **Step 2: Implement Rao-Blackwellized particle filter**

Hybrid: particles for nonlinear states, Kalman filter for linear substates. Reduces variance.

- [ ] **Step 3: Tests and commit**

```bash
git commit -m "feat: add auxiliary and Rao-Blackwellized particle filters"
```

### Task 3.4d: Add Smoothers (RTS, URTSS) and Batch (WLS, MAP)

**Files:**
- Modify: `src/fortran/forapollo_estimate.f90`
- Modify: `tests/fortran/test_estimate.f90`

- [ ] **Step 1: Implement RTS smoother**

Backward pass over stored forward KF results. Takes arrays of x_pred, x_upd, P_pred, P_upd, F for N timesteps.

- [ ] **Step 2: Implement URTSS (unscented smoother)**

Same backward pass but sigma-point based.

- [ ] **Step 3: Implement batch WLS and MAP**

WLS: minimize `sum ||z_i - h(x)||^2_R`. MAP: add prior `||x - x0||^2_P0`. Uses forOpt for optimization.

- [ ] **Step 4: Test: smoother should give tighter estimates than filter alone**
- [ ] **Step 5: Commit**

```bash
git commit -m "feat: add RTS/URTSS smoothers and batch WLS/MAP estimators"
```

### Task 3.5: Create Zig Safety Layer

**Files:**
- Create: `src/zig/fortran.zig` — extern fn declarations for all fa_* symbols
- Create: `src/zig/safety.zig` — bounds checking, error mapping
- Create: `src/zig/dispatch.zig` — size-based Zig/Fortran routing
- Create: `src/zig/kernels.zig` — Zig-native small-state implementations

- [ ] **Step 1: Write fortran.zig — extern declarations**

Declare every `fa_*` symbol as `extern "c" fn`. One-to-one with Fortran bind(C) exports.

- [ ] **Step 2: Write safety.zig — validation layer**

```zig
pub fn validateEstimatorArgs(n: i32, m: i32, x: ?[*]f64, P: ?[*]f64) !void {
    if (n <= 0) return error.InvalidDimension;
    if (x == null) return error.NullPointer;
    if (P == null) return error.NullPointer;
    // bounds check: P must be n*n elements
}
```

- [ ] **Step 3: Write dispatch.zig — size-based routing**

Small state (n < threshold) → Zig-native implementation.
Large state (n >= threshold) → Fortran via extern fn.
Threshold configurable per platform (Desktop: 12, Jetson: 6).

- [ ] **Step 4: Write kernels.zig — Zig-native small KF**

For n <= 6 (covers most tracking problems), implement KF/EKF directly in Zig for zero-FFI-overhead.

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: add Zig safety layer (extern, dispatch, kernels, safety)"
```

---

## Phase 4: Engine Extended (Guidance + Coords)

### Task 4.1a: Implement ZEM/ZEV and Proportional Navigation

**Files:**
- Create: `src/fortran/forapollo_guidance.f90`
- Create: `tests/fortran/test_guidance.f90`

- [ ] **Step 1: Write ZEM/ZEV test**

Test: given known position error and tgo, verify commanded acceleration direction/magnitude.

- [ ] **Step 2: Implement ZEM/ZEV and E-guidance**

ZEM/ZEV: `a_cmd = -C1*(r - r_tgt)/tgo^2 - C2*(v - v_tgt)/tgo` (Apollo P63/P64 core)

- [ ] **Step 3: Implement proportional navigation (pure, augmented, true)**

PN: `a_cmd = N * V_c * omega` (closing velocity × LOS rate)

- [ ] **Step 4: Tests and commit**

```bash
git commit -m "feat: add ZEM/ZEV, E-guidance, and proportional navigation"
```

### Task 4.1b: Implement Polynomial Guidance (Gravity Turn, PEG, Linear Tangent)

**Files:**
- Modify: `src/fortran/forapollo_guidance.f90`
- Modify: `tests/fortran/test_guidance.f90`

- [ ] **Step 1: Implement gravity turn, linear tangent steering, PEG**
- [ ] **Step 2: Tests and commit**

```bash
git commit -m "feat: add polynomial guidance (gravity turn, linear tangent, PEG)"
```

### Task 4.1c: Implement Optimal Control (LQR, iLQR, DDP)

**Files:**
- Modify: `src/fortran/forapollo_guidance.f90`
- Modify: `tests/fortran/test_guidance.f90`

- [ ] **Step 1: Implement LQR**

Solve Riccati equation `A^T*P + P*A - P*B*R^-1*B^T*P + Q = 0`. Uses forMath linalg.

- [ ] **Step 2: Implement iLQR (iterative) and DDP**

Forward-backward iterative linearization. Uses forOpt for line search.

- [ ] **Step 3: Test: double integrator LQR — verify gain matrix matches analytical solution**
- [ ] **Step 4: Commit**

```bash
git commit -m "feat: add optimal control (LQR, iLQR, DDP)"
```

### Task 4.1d: Implement MPC (Shooting + Collocation)

**Files:**
- Modify: `src/fortran/forapollo_guidance.f90`
- Modify: `tests/fortran/test_guidance.f90`

- [ ] **Step 1: Implement MPC via shooting method**

Forward simulate with candidate controls, optimize via forOpt.

- [ ] **Step 2: Implement MPC via direct collocation**

Discretize trajectory into NLP, solve via forOpt SLSQP.

- [ ] **Step 3: Tests and commit**

```bash
git commit -m "feat: add model predictive control (shooting, collocation)"
```

### Task 4.1e: Implement Targeting and Path Following

**Files:**
- Modify: `src/fortran/forapollo_guidance.f90`
- Modify: `tests/fortran/test_guidance.f90`

- [ ] **Step 1: Implement targeting (single/multi shooting, Lambert, differential correction)**
- [ ] **Step 2: Implement path following (pure pursuit, Stanley, trajectory tracking)**
- [ ] **Step 3: Tests and commit**

```bash
git commit -m "feat: add targeting methods and path following controllers"
```

### Task 4.1f: Implement Energy-Optimal Transfers + forGraph Integration

**Files:**
- Modify: `src/fortran/forapollo_guidance.f90`
- Modify: `tests/fortran/test_guidance.f90`

- [ ] **Step 1: Implement min-fuel, min-time, min-energy transfers**
- [ ] **Step 2: Integrate forGraph for waypoint sequencing (A*/Dijkstra)**

Waypoint graph → `fgr_astar` for optimal ordering → trajectory tracking.

- [ ] **Step 3: Tests and commit**

```bash
git commit -m "feat: add energy-optimal transfers and forGraph waypoint integration"
```

### Task 4.2: Implement forapollo_coords.f90

**Files:**
- Create: `src/fortran/forapollo_coords.f90`
- Create: `tests/fortran/test_coords.f90`

- [ ] **Step 1: Implement frame transform dispatch + ECI/ECEF/geodetic**

```fortran
subroutine fa_frame_transform(from_id, to_id, n, r_in, r_out, ref_params, nrp, t, info)
    ! Dispatch on (from_id, to_id) pair
end subroutine
```

ECI→ECEF requires Earth rotation angle (from time). Calls forMath quaternions.

- [ ] **Step 2: Implement NED, ENU, LVLH, orbital frames**

- [ ] **Step 3: Implement geodetic conversions (WGS84)**

```fortran
subroutine fa_geodetic_from_ecef(r_ecef, a, f, lat, lon, alt, info)
    ! Iterative Bowring method for lat/lon/alt from ECEF
end subroutine
```

- [ ] **Step 4: Tests and commit**

```bash
git commit -m "feat: add coordinate frame catalog and transforms"
```

### Task 4.3: Complete Zig Exports

**Files:**
- Modify: `src/zig/exports.zig` — add all forapollo_* C ABI exports

- [ ] **Step 1: Export all estimator functions through safety layer**

Each `forapollo_ekf_predict` calls `safety.validate()` then `dispatch.ekfPredict()`.

- [ ] **Step 2: Export guidance, coords, dynamics, observe**

- [ ] **Step 3: Verify zig build compiles**

```bash
make all && zig build
```

- [ ] **Step 4: Commit**

```bash
git commit -m "feat: add complete Zig C ABI exports for all fa_* functions"
```

---

## Phase 5: Domain Layer (Astro + Environment + Time)

### Task 5.1: Implement forapollo_astro.f90

**Files:**
- Create: `src/fortran/forapollo_astro.f90`
- Create: `tests/fortran/test_astro.f90`

- [ ] **Step 1: Implement orbital element conversions**

`fa_elements_from_rv`: Cartesian → (a, e, i, Omega, omega, nu)
`fa_rv_from_elements`: inverse
Test: round-trip conversion should recover original state.

- [ ] **Step 2: Implement vis-viva, period, SOI, transfer sizing**

- [ ] **Step 3: Implement planetary constants and eclipse geometry**

- [ ] **Step 4: Tests and commit**

```bash
git commit -m "feat: add astrodynamics utilities (elements, transfers, eclipse)"
```

### Task 5.2: Implement forapollo_environ.f90

**Files:**
- Create: `src/fortran/forapollo_environ.f90`
- Create: `tests/fortran/test_environ.f90`

- [ ] **Step 1: Implement US Standard Atmosphere 1976**

Altitude → density, pressure, temperature. Seven-layer piecewise model.

- [ ] **Step 2: Implement J2/J4 gravity and spherical harmonics**

J2 perturbation acceleration (uses forFFT for high-degree harmonics).

- [ ] **Step 3: Implement SRP and geodesics**

- [ ] **Step 4: Tests and commit**

```bash
git commit -m "feat: add environment models (atmosphere, gravity, SRP, geodesy)"
```

### Task 5.3: Implement forapollo_time.f90

**Files:**
- Create: `src/fortran/forapollo_time.f90`
- Create: `tests/fortran/test_time.f90`

- [ ] **Step 1: Implement UTC/TAI/TT/TDB/GPS conversions**

TAI = UTC + leap_seconds. TT = TAI + 32.184. TDB ≈ TT + periodic correction.

- [ ] **Step 2: Implement Julian date, MJD, Unix epoch**

- [ ] **Step 3: Implement leap second table**

Hardcoded table through current date. Comment noting staleness risk.

- [ ] **Step 4: Tests and commit**

```bash
git commit -m "feat: add time system conversions (future: migrates to forTime)"
```

---

## Phase 6: Python Bindings

### Task 6.1: Setup Python Package

**Files:**
- Create: `python/forapollo/__init__.py`
- Create: `python/forapollo/_lib.py`
- Create: `setup.py`
- Create: `requirements.txt`

- [ ] **Step 1: Create requirements.txt**

```
numpy>=1.21
```

- [ ] **Step 2: Create setup.py**

```python
# Copyright The Fantastic Planet — By David Clabaugh
from setuptools import setup, find_packages

setup(
    name="forapollo",
    version="0.1.0",
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    install_requires=["numpy>=1.21"],
    python_requires=">=3.10",
)
```

- [ ] **Step 3: Create _lib.py — ctypes loader**

Loads `libforapollo.so`/`.dylib`, declares all `forapollo_*` function signatures.

- [ ] **Step 4: Create __init__.py**

```python
# forApollo — Universal State Estimation & Guidance
# Copyright The Fantastic Planet — By David Clabaugh
from forapollo.estimators import EKF, UKF, ESKF, ParticleFilter
from forapollo.dynamics import KEPLER, CONST_VEL, BICYCLE, QUADROTOR
from forapollo.sensors import POSITION, RANGE, BEARING, RADAR, SCALAR
```

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: add Python package scaffold with ctypes loader"
```

### Task 6.2: Implement Python Estimator Classes

**Files:**
- Create: `python/forapollo/estimators.py`
- Create: `python/forapollo/dynamics.py`
- Create: `python/forapollo/sensors.py`
- Create: `tests/python/test_estimators.py`

- [ ] **Step 1: Implement EKF class**

```python
class EKF:
    def __init__(self, state_dim, meas_dim):
        self.n = state_dim
        self.m = meas_dim
        self.x = np.zeros(state_dim)
        self.P = np.eye(state_dim).flatten()  # flat row-major
        # ...

    def predict(self, dt, model=None, f_ptr=None, params=None):
        # calls _lib.forapollo_ekf_predict
        pass

    def update(self, z, sensor=None, h_ptr=None, params=None):
        # calls _lib.forapollo_ekf_update
        # returns validity array (ternary: -1/0/+1)
        pass
```

- [ ] **Step 2: Implement dynamics.py with model constants and @custom decorator**

- [ ] **Step 3: Implement sensors.py with observation constants**

- [ ] **Step 4: Write Python test — EKF tracking with const-vel + position obs**

```python
def test_ekf_const_vel_tracking():
    rng = np.random.default_rng(42)  # fixed seed for reproducibility
    ekf = EKF(state_dim=4, meas_dim=2)
    ekf.set_dynamics(CONST_VEL)
    ekf.set_sensor(POSITION)
    ekf.init(x0=np.array([0, 0, 1, 0]), P0=np.eye(4)*100)
    ekf.set_process_noise(np.diag([0.01, 0.01, 0.1, 0.1]))
    ekf.set_meas_noise(np.diag([1.0, 1.0]))

    # Simulate: true trajectory is constant velocity [1, 0]
    for i in range(20):
        ekf.predict(dt=1.0)
        z = np.array([float(i+1), 0.0]) + rng.standard_normal(2) * 1.0
        ekf.update(z)

    # After 20 steps, x[0] should be near 20, x[2] should be near 1
    assert abs(ekf.x[0] - 20.0) < 5.0, f"Position diverged: {ekf.x[0]}"
    assert abs(ekf.x[2] - 1.0) < 1.0, f"Velocity diverged: {ekf.x[2]}"
```

- [ ] **Step 5: Run test**

```bash
cd python && python -m pytest tests/python/test_estimators.py -v
```

- [ ] **Step 6: Commit**

```bash
git commit -m "feat: add Python EKF class with ctypes bindings"
```

### Task 6.3: Implement Remaining Python Wrappers

**Files:**
- Create: `python/forapollo/guidance.py`
- Create: `python/forapollo/coords.py`
- Create: `python/forapollo/astro.py`
- Create: `python/forapollo/environ.py`
- Create: `python/forapollo/time.py`
- Create: `tests/python/test_integration.py`

- [ ] **Step 1: Implement guidance.py, coords.py**

- [ ] **Step 2: Implement astro.py, environ.py, time.py**

- [ ] **Step 3: Write integration test — full orbit determination**

Test: simulate a satellite orbit, generate noisy range measurements, run EKF, verify convergence.

- [ ] **Step 4: Write integration test — robot tracking**

Test: bicycle model, range+bearing observations, UKF. Verify convergence.

- [ ] **Step 5: Run all tests, commit**

```bash
git commit -m "feat: complete Python bindings for all forApollo modules"
```

---

## Post-Implementation

### Task 7.1: Push to forKernels/forApollo

- [ ] **Step 1: Verify all tests pass**

```bash
make all && make test   # Fortran tests
zig build test          # Zig tests
cd python && python -m pytest tests/ -v  # Python tests
```

- [ ] **Step 2: Push**

```bash
git push -u origin main
```

---

## Dependency Checklist

Before Phase 3 (Engine) can fully compile and link, you need real prebuilt archives from upstream repos. The Zig stubs (Phase 1, Task 1.4) allow compilation and testing of the Fortran kernels independently, but the Zig dispatch layer needs real `.a` files.

| Dependency | Needed by Phase | Symbols used |
|------------|----------------|-------------|
| forMath linalg | Phase 3 | `fm_dgemm`, `fm_dgemv`, matrix inversion |
| forMath ODE | Phase 3 | `fm_ode_rk4_step`, `fm_ode_dp5_step` |
| forMath quaternion | Phase 4 | `fm_quat_multiply`, `fm_quat_to_dcm` |
| forMath numdiff | Phase 3 | `fm_numdiff_jacobian` |
| forMath random | Phase 3 | `fm_random_normal` (particle filter) |
| forFFT | Phase 5 | Spherical harmonic expansion |
| forOpt | Phase 4 | MPC solver, Riccati equation |
| forTernary | Phase 3 | `fk_ternary_gate` (measurement gating) |
| forGraph | Phase 4 | `fgr_astar`, `fgr_dijkstra` (guidance path planning) |

Copyright The Fantastic Planet — By David Clabaugh
