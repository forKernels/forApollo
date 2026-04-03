# forKernels Build System Fix — Paste Into Each Repo's Agent

## What you need to fix

This repo's `build.zig` needs two changes:

### 1. Per-target output directories

Every `b.installArtifact(lib)` call must route to `zig-out/{target}/lib/` instead of flat `zig-out/lib/`. This prevents macOS builds from overwriting Linux builds.

**Add this helper function** (if not already present):

```zig
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
            .aarch64 => "windows-x86_64",
            else => "windows-unknown",
        },
        else => "unknown",
    };
}
```

**Replace every** `b.installArtifact(lib);` with:

```zig
const target_name = getTargetName(target.result);  // compute once near top of build()

// Then for each library:
{
    const install = b.addInstallArtifact(lib, .{
        .dest_dir = .{ .override = .{ .custom = b.fmt("{s}/lib", .{target_name}) } },
    });
    b.getInstallStep().dependOn(&install.step);
}
```

This makes `zig build` output to `zig-out/macos-arm64/lib/` or `zig-out/linux-arm64/lib/` depending on target. Multiple targets coexist without overwriting.

### 2. Sibling dep resolution (no more `deps/` dir)

All sibling library dependencies must resolve from sibling repo prebuilt directories, NOT from a local `deps/` folder. The resolution order is:

```
1. ../sibling/prebuilt/{target}/lib/     (sibling's checked-in prebuilts)
2. ../sibling/zig-out/{target}/lib/      (sibling's dev build output)
3. ../sibling/prebuilt/lib/              (legacy flat — fallback only)
4. ../sibling/zig-out/lib/               (legacy flat — fallback only)
```

**Pattern for build.zig:**

```zig
const target_name = getTargetName(target.result);

// For each sibling dependency, add all resolution paths:
inline for (.{
    .{ "forMath", "../forMath" },
    .{ "forOpt", "../forOpt" },
    // ... add whichever siblings this repo depends on
}) |dep| {
    const sibling = dep[1];
    // Target-specific paths (preferred)
    lib.addLibraryPath(.{ .cwd_relative = sibling ++ "/prebuilt/" ++ target_name ++ "/lib" });
    lib.addLibraryPath(.{ .cwd_relative = sibling ++ "/zig-out/" ++ target_name ++ "/lib" });
    // Legacy flat paths (fallback)
    lib.addLibraryPath(.{ .cwd_relative = sibling ++ "/prebuilt/lib" });
    lib.addLibraryPath(.{ .cwd_relative = sibling ++ "/zig-out/lib" });
}

// Then link by name as usual:
lib.linkSystemLibrary("formath_linalg");
lib.linkSystemLibrary("foropt");
// etc.
```

**Remove** any `deps/` references. We don't use a `deps/` directory anymore. All repos are siblings in the same parent directory and resolve each other directly.

### 3. Fortran object resolution

If this repo has Fortran sources, the prebuilt `.o` files also use per-target paths:

```
prebuilt/obj/{target}/       ← prebuilt Fortran objects
build/{target}/obj/          ← locally compiled Fortran objects
```

Build.zig should check both locations:

```zig
const prebuilt_obj_dir = b.fmt("prebuilt/obj/{s}", .{target_name});
const build_obj_dir = b.fmt("build/{s}/obj", .{target_name});
```

## Target names (standardized across all repos)

| Target ID | Zig triple | Use for |
|-----------|-----------|---------|
| `macos-arm64` | `aarch64-macos` | Apple M1/M2/M3/M4 |
| `linux-arm64` | `aarch64-linux-gnu` | Google Cloud Axion, Graviton, Jetson, RPi5 |
| `windows-x86_64` | `x86_64-windows-gnu` | Windows x86_64 |

**Do NOT use:** `darwin-arm64`, `aarch64-macos`, `linux-aarch64`, `linux-x64`, `macos-arm64-libs`, or any other variant. Use exactly the names above.

## Build commands

```bash
# Native build (macOS arm64 on M4)
zig build -Doptimize=ReleaseFast

# Cross-compile for Linux ARM
zig build -Doptimize=ReleaseFast -Dtarget=aarch64-linux-gnu

# Per-target output (using --prefix)
zig build -Doptimize=ReleaseFast --prefix zig-out/macos-arm64
zig build -Doptimize=ReleaseFast -Dtarget=aarch64-linux-gnu --prefix zig-out/linux-arm64
```

## Directory structure (standard for all repos)

```
repo/
  build.zig              # Always at repo root
  src/
    fortran/             # Fortran sources
    zig/                 # Zig sources
  prebuilt/
    macos-arm64/lib/     # Checked-in deliverables
    linux-arm64/lib/
    windows-x86_64/lib/
    obj/                 # Prebuilt Fortran objects (if applicable)
      macos-arm64/
      linux-arm64/
  zig-out/               # Dev build output (gitignored)
    macos-arm64/lib/
    linux-arm64/lib/
```

## Common issues

**"file not found: deps/formath/lib/..."** → Remove `deps/` references. Use sibling resolution pattern above.

**"libformath.a not found"** → forMath doesn't produce a monolithic `libformath.a`. It produces modular archives: `libformath_linalg.a`, `libformath_optimize.a`, `libformath_sparse.a`, etc. Link the specific ones you need.

**"prebuilt/obj/macos-arm64/foo.o not found"** → Fortran objects need building first. Run `make MODE=release` or add a Fortran compile step to build.zig.

**"unable to find dynamic system library"** → Library search path is wrong. Add the correct `addLibraryPath` for sibling prebuilt dirs.

**Builds overwrite each other** → Use `--prefix zig-out/{target}` or add `dest_dir` override to `addInstallArtifact`.
