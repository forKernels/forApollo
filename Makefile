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

FC       ?= gfortran
MARCH    ?= native

FFLAGS   := -O3 -ftree-vectorize -fPIC -cpp -std=f2008 -fall-intrinsics \
            -Wall -Wextra -Wimplicit-interface -fopenmp -march=$(MARCH)

# --- Directories --------------------------------------------------------------

SRC_DIR   := src/fortran
BUILD_DIR := build
OBJ_DIR   := $(BUILD_DIR)/obj
TEST_DIR  := tests/fortran

# --- Platform detection for prebuilt target -----------------------------------
#
# Darwin (macOS) with Apple Silicon → macos-arm64
# Linux aarch64 (Jetson, RPi)      → linux-arm64
# Linux x86_64 (cloud, CI)         → linux-x86_64
#
UNAME_S  := $(shell uname -s)
UNAME_M  := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
    PLATFORM := macos-arm64
else ifeq ($(UNAME_M),x86_64)
    PLATFORM := linux-x86_64
else
    PLATFORM := linux-arm64
endif

PREBUILT_DIR := prebuilt/$(PLATFORM)

# --- Output archive -----------------------------------------------------------
#
# Archive goes into build/{platform}/lib/ to match Zig build.zig resolution.
#
PLAT_LIB_DIR := $(BUILD_DIR)/$(PLATFORM)/lib
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
