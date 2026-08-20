# =============================================================================
# forApollo — Fortran Compilation Makefile (Stage 1 of two-stage build)
# Copyright The Fantastic Planet — By David Clabaugh
#
# This Makefile compiles Fortran kernels to a static archive.
# Stage 2 (Zig safety layer + C ABI exports) is handled by build.zig.
#
# Usage:
#   make                        # compile all Fortran sources → libforapollo_fortran.a
#   make lib                    # same as all
#   make test                   # compile and run Fortran unit tests
#   make prebuilt               # copy archive to prebuilt/<platform>/
#   make clean                  # remove build/ artifacts
#
# Cross-compilation:
#   make MARCH=armv8-a          # override -march for Jetson / RPi
# =============================================================================

# --- Configurable flags -------------------------------------------------------

# NOTE: plain `=`, not `?=` — FC is a make BUILTIN (defaults to f77), so `?=`
# never fires. Command-line `make FC=...` still overrides.
FC       = gfortran
MARCH    ?= native

# =============================================================================
# -fwrapv: DEFINED WRAPPING FOR SIGNED OVERFLOW.
#
# src/fortran/forapollo_estimate.f90 runs a glibc-style LCG on 4-byte integers
# (integer(c_int) :: seed, and a local `integer :: local_seed`):
#     seed = mod(seed * 1103515245 + 12345, 2147483647)
# The multiply overflows signed 32-bit on essentially every call. The mod()
# bounds the RESULT, not the intermediate product, so it does not prevent it.
# Signed overflow is UNDEFINED in Fortran, so the optimizer may assume it never
# happens. Site enumerated by the multiplier constant present in the source, not
# from a site list (the 2026-08-19 list named a forMath file that does not exist
# and missed four forQuant files outright); this is forApollo's only one.
#
# NOT demonstrably miscompiled today: for this class of generator the isolated
# stream hash is identical across -O0/-O1/-O2/-O3 native/-Ofast and with
# -fwrapv, while -ftrapv aborts at the multiply. The overflow executes; gfortran
# is not currently exploiting it here. Do not report this as a repair.
#
# SO WHY ADD IT. Because that isolation test has a DEMONSTRATED BLIND SPOT.
# for3D's tests/test_density_to_mesh.f90 carries the same class of constant and
# IS live miscompiled on gfortran 15.2.0 -- the real fixture HANGS at -O2 and
# -O3 without -fwrapv -- while a standalone replication of that exact LCG,
# measured the same way, hashes IDENTICALLY at every level. The construct that
# provably hangs in situ looks clean in isolation, because the miscompile
# depends on surrounding context rather than on the LCG alone. A negative from
# the isolated test therefore licenses nothing.
#
# An estimator that silently reuses a frozen random draw degrades quietly --
# particle filters and UKF sigma-point jitter would still return numbers. The
# failure mode is silence: for3D presented as a HANG, not as bad numbers.
#
# REBUILDING: make keys off mtime, never off flags, so this triggers no
# recompile by itself. `make clean` first or it is inert.
# =============================================================================
FFLAGS_OVERFLOW := -fwrapv

FFLAGS   := $(FFLAGS_OVERFLOW) -O3 -ftree-vectorize -fPIC -cpp -std=f2008 -fall-intrinsics \
            -Wall -Wextra -Wimplicit-interface -fopenmp -march=$(MARCH)

# --- Target detection — short canonical names --------------------------------
#
# forKernels canon: thor | linX86 | winX86 | macos (NO long names).
# build.zig getTargetName() emits the SAME short names; the obj dir below MUST
# match build.zig's Stage-2 read path (prebuilt/<target>/obj). See
# ../forMath/docs/DELIVERY.md.
#
# Host auto-detection covers thor/linX86/macos. winX86 is a cross target:
#   make TARGET=winX86
# Stage 2 (build.zig) invokes `make lib TARGET=<target>` so Stage 1 writes to
# the same per-target obj dir Stage 2 reads — all 4 targets build, no clobber.
#
UNAME_S  := $(shell uname -s)
UNAME_M  := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
    DETECTED_TARGET := macos
else ifeq ($(UNAME_M),x86_64)
    DETECTED_TARGET := linX86
else
    DETECTED_TARGET := thor
endif

TARGET ?= $(DETECTED_TARGET)

# --- Directories --------------------------------------------------------------

SRC_DIR      := src/fortran
BUILD_DIR    := build
TEST_DIR     := tests/fortran
PREBUILT_DIR := prebuilt/$(TARGET)

# Fortran objects → prebuilt/<target>/obj — the SAME path build.zig (Stage 2)
# reads via `prebuilt/{target}/obj`. Per-target so `make TARGET=t` never
# clobbers another target's objects. (build intermediates: prebuilt/**/obj/
# is gitignored; only the final prebuilt/<target>/libforapollo.a is committed.)
OBJ_DIR      := $(PREBUILT_DIR)/obj

# --- Output archive -----------------------------------------------------------
#
# Intermediate Fortran archive goes into build/{target}/lib/ to match Zig
# build.zig's shared/test resolution (build/<target>/lib/libforapollo_fortran.a).
#
PLAT_LIB_DIR := $(BUILD_DIR)/$(TARGET)/lib
LIB          := $(PLAT_LIB_DIR)/libforapollo_fortran.a

# --- Source compilation order -------------------------------------------------
#
# Module dependency order (must compile lower layers before upper layers):
#
#   Layer 0 — Models (no intra-Apollo deps): dynamics, observe
#   Layer 1 — Engine core: estimate, propagate
#   Layer 2 — Engine aux: guidance, coords
#   Layer 3 — Domain utilities: astro, environ, time
#
SRCS := \
    $(SRC_DIR)/forapollo_dynamics.f90   \
    $(SRC_DIR)/forapollo_observe.f90    \
    $(SRC_DIR)/forapollo_estimate.f90   \
    $(SRC_DIR)/forapollo_propagate.f90  \
    $(SRC_DIR)/forapollo_guidance.f90   \
    $(SRC_DIR)/forapollo_coords.f90     \
    $(SRC_DIR)/forapollo_astro.f90      \
    $(SRC_DIR)/forapollo_environ.f90    \
    $(SRC_DIR)/forapollo_time.f90

# Derive object file names: src/fortran/forapollo_X.f90 → build/obj/forapollo_X.o
OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

# --- Test sources (tests/fortran/test_*.f90) ----------------------------------

TEST_SRCS := $(wildcard $(TEST_DIR)/test_*.f90)
TEST_BINS := $(patsubst $(TEST_DIR)/%.f90,$(BUILD_DIR)/%,$(TEST_SRCS))

# =============================================================================
# Targets
# =============================================================================

.PHONY: all lib prebuilt clean test

# Default target
all: lib

# --- lib: compile Fortran kernels → static archive ---------------------------

lib: $(LIB)

$(LIB): $(OBJS) | $(PLAT_LIB_DIR)
	@echo "[AR] $@"
	ar rcs $@ $^

$(PLAT_LIB_DIR):
	mkdir -p $(PLAT_LIB_DIR)

# Pattern rule: compile each .f90 → .o
# Modules (.mod) are emitted into OBJ_DIR so they are found by dependent units.
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	@echo "[FC] $<"
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@

# Ensure build directories exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# --- prebuilt: copy archive to prebuilt/<platform>/ --------------------------
#
# Called when generating a new committed prebuilt snapshot.
# Destination directory must already exist (created as part of repo layout).

prebuilt: $(LIB)
	@echo "[PREBUILT] copying to $(PREBUILT_DIR)/lib/"
	@mkdir -p $(PREBUILT_DIR)/lib
	cp $(LIB) $(PREBUILT_DIR)/lib/libforapollo_fortran.a
	@echo "[PREBUILT] done — $(PREBUILT_DIR)/lib/libforapollo_fortran.a"

# --- test: compile and run Fortran unit tests --------------------------------
#
# Each tests/fortran/test_*.f90 is compiled into its own executable and linked
# against the objects (not the archive, so there is no separate lib dep).
# Tests are run in sequence; any non-zero exit fails the target.

test: $(OBJS) $(TEST_BINS)
	@echo "[TEST] running $(words $(TEST_BINS)) test(s)"
	@failed=0; \
	for t in $(TEST_BINS); do \
	    echo "[RUN] $$t"; \
	    $$t; \
	    if [ $$? -ne 0 ]; then \
	        echo "[FAIL] $$t"; \
	        failed=$$((failed + 1)); \
	    else \
	        echo "[PASS] $$t"; \
	    fi; \
	done; \
	if [ $$failed -ne 0 ]; then \
	    echo "[TEST] $$failed test(s) failed"; \
	    exit 1; \
	else \
	    echo "[TEST] all tests passed"; \
	fi

# Pattern rule: compile test binary
# Links against the object files so module .mod files are available.
$(BUILD_DIR)/%: $(TEST_DIR)/%.f90 $(OBJS) | $(BUILD_DIR)
	@echo "[FC-TEST] $<"
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -I$(OBJ_DIR) $< $(OBJS) -o $@

# --- clean: remove all build artifacts ----------------------------------------

clean:
	@echo "[CLEAN] removing $(BUILD_DIR)/"
	rm -rf $(BUILD_DIR)
