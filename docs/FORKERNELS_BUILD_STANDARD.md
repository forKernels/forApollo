# forKernels Build Standard

## Target Platforms

Two targets. That's it.

| Target ID | Triple | Hardware |
|-----------|--------|----------|
| `macos-arm64` | `aarch64-macos` | Apple M1/M2/M3/M4 |
| `linux-arm64` | `aarch64-linux-gnu` | Graviton, Jetson, RPi5, new ARM servers |
| `windows-x86_64` | `x86_64-windows-gnu` | Windows x86_64 |

ARM only. No x86 targets.

## Directory Structure (every repo)

```
repo/
  build.zig              # Always in repo root (not zig/ or src/zig/)
  Makefile               # Fortran build (if repo has Fortran)
  src/
    fortran/             # Fortran sources (*.f90)
    zig/                 # Zig sources (*.zig)
  prebuilt/
    macos-arm64/
      lib/
        lib{name}.a      # Deliverable static archive
    linux-arm64/
      lib/
        lib{name}.a      # Deliverable static archive
  zig-out/               # Dev build output (gitignored)
    macos-arm64/
      lib/
        lib{name}.a
    linux-arm64/
      lib/
        lib{name}.a
```

### Rules

1. **`build.zig` lives at repo root.** Not in `zig/`, not in `src/zig/`. If it's elsewhere, move it.

2. **`prebuilt/` uses `{target}/lib/` layout.** No `darwin-arm64`, no `aarch64-macos`, no `macos-arm64-libs`. Exactly `macos-arm64` and `linux-arm64`.

3. **`zig-out/` uses the same `{target}/lib/` layout.** Pass `--prefix zig-out/{target}` to `zig build`. Never build into flat `zig-out/lib/`.

4. **Only deliverables go in `prebuilt/`.** The final `.a` archive(s) that downstream repos link. Not intermediate `.o` files, not debug artifacts, not test binaries.

5. **Naming convention:** `lib{reponame}.a` in lowercase. Examples:
   - forMath → `libformath_linalg.a`, `libformath_optimize.a`, etc.
   - forIO → `libforio.a`, `libforio_async.a`, etc.
   - forTorch → `libfortorch.a`
   - forML → `libforml.a`
   - forGraph → `libforgraph.a`

6. **Fortran objects stay internal.** If a repo has Fortran, the Zig build compiles Fortran via `make` or `zig cc`, then links into the final `.a`. Downstream repos never see `.o` files.

7. **No `_standalone`, `_combined`, `_pack` suffixes** unless the repo genuinely ships multiple libraries with different dependency sets. One repo = one library (or a small set of named modules like forMath).

## build.zig Standard

Every `build.zig` must:

```zig
// 1. Use --prefix for target-specific output
// (zig build handles this via b.install_prefix)

// 2. Resolve sibling deps from target-specific paths
const target_str = if (is_macos)
    if (is_arm64) "macos-arm64" else "macos-x86_64"
else if (is_linux)
    if (is_arm64) "linux-arm64" else "linux-x86_64"
else "unknown";

// 3. Search order for each dependency:
//    a) prebuilt/{target}/lib/     (checked-in prebuilts)
//    b) ../sibling/prebuilt/{target}/lib/  (sibling prebuilts)
//    c) ../sibling/zig-out/{target}/lib/   (sibling dev builds)
```

## Build Commands

```bash
# Native build (auto-detects target)
zig build -Doptimize=ReleaseFast --prefix zig-out/macos-arm64

# Cross-compile for Linux ARM (from macOS)
zig build -Doptimize=ReleaseFast -Dtarget=aarch64-linux-gnu --prefix zig-out/linux-arm64

# Docker build for Linux ARM (native in container)
docker run --platform linux/arm64 -v $(pwd):/src -w /src \
  debian:bookworm-slim bash -c "zig build -Doptimize=ReleaseFast --prefix zig-out/linux-arm64"

# Copy dev build to prebuilt (for commit)
cp -r zig-out/macos-arm64/lib/* prebuilt/macos-arm64/lib/
cp -r zig-out/linux-arm64/lib/* prebuilt/linux-arm64/lib/
```

## Makefile Standard (for Fortran repos)

```makefile
# Build Fortran objects for current platform
# Output: build/{target}/obj/*.o and build/{target}/lib/lib{name}.a
PLATFORM := $(shell uname -s | tr '[:upper:]' '[:lower:]')
ARCH := $(shell uname -m)
ifeq ($(ARCH),arm64)
  ARCH := arm64
endif
ifeq ($(ARCH),aarch64)
  ARCH := arm64
endif
TARGET := $(PLATFORM)-$(ARCH)
BUILD_DIR := build/$(TARGET)
```

## Downstream Resolution

When repo A depends on repo B, A's `build.zig` searches:

```
1. A/prebuilt/{target}/lib/          (A vendors B's prebuilts — for CI)
2. ../B/prebuilt/{target}/lib/       (B's checked-in prebuilts — for local dev)
3. ../B/zig-out/{target}/lib/        (B's dev build output — for active dev)
```

This means:
- CI/release: each repo vendors its deps in `prebuilt/`
- Local dev: sibling repos resolve automatically
- No `deps/` directory needed — removed in favor of direct sibling resolution

## Migration Checklist

For each repo:
- [ ] `build.zig` at repo root
- [ ] `prebuilt/macos-arm64/lib/` and `prebuilt/linux-arm64/lib/`
- [ ] Remove old prebuilt dirs (`darwin-arm64`, `aarch64-macos`, `*-libs`, etc.)
- [ ] `zig build --prefix zig-out/{target}` works
- [ ] No flat `zig-out/lib/` output
- [ ] Library names follow `lib{name}.a` convention
- [ ] Fortran objects are internal (not exposed to downstream)
- [ ] Doc dropped in `docs/FORKERNELS_BUILD_STANDARD.md`
