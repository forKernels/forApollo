## DEPENDENCY RULES — READ BEFORE WRITING ANY CODE

### THE ZIG LAYER IS THE API — NOT FORTRAN

When this library needs a function from forMath/forFFT/forBayes/forIO/forSim/etc:

  CORRECT:  const formath = @import("formath");
            const result = formath.linalg.dgemm(...);

  WRONG:    extern fn fm_dgemm(...) callconv(.C) void;
            (DO NOT redeclare extern fn — use the upstream Zig wrapper)

  WRONG:    Link formath_linalg.f90 into this build
            (DO NOT compile upstream Fortran source)

  WRONG:    Look in prebuilt/ for .o files and call symbols directly
            (DO NOT bypass the Zig safety layer)

The dependency chain is:
  upstream .f90 → upstream prebuilt .a → upstream Zig wrapper → THIS library imports Zig wrapper

Never skip layers. Never go direct to Fortran.

### LINKING RULES

1. NEVER add another library's .f90 files to this build.
   Link their prebuilt .a archive from their zig-out/lib/ directory.

2. NEVER rewrite Zig wrappers that exist in upstream libraries.
   Import them: const formath = @import("formath");

3. If you need a symbol from forMath/forFFT/forBayes/etc:
   - Check if it exists in their Zig module first (this is the API)
   - Import their Zig module for the safe wrapper
   - The prebuilt .a archive provides the Fortran symbols at link time
   - Do NOT copy their Fortran source into this repo
   - Do NOT write a new extern fn declaration — use theirs
   - Do NOT write a new Zig safe wrapper — use theirs

4. The only .f90 files in this repo are THIS library's kernels.
   Everything else comes from prebuilt archives via Zig imports.

5. Default build links prebuilt deps from zig-out/lib/. Source build is opt-in: -Ddev=true

### WHY THIS MATTERS

- Duplicated Fortran source across repos causes symbol collisions at link time
- Duplicated Zig wrappers create API inconsistencies
- Building from source cascades recompilation across the entire dependency tree
- Prebuilt archives are fast (milliseconds to link) and deterministic
