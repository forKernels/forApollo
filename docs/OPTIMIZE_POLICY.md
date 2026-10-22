# Optimize Policy - forApollo

**Classification: ZORTRAN -> delivery default `ReleaseSafe`**

Decided 2026-08-04. Supersedes the 2026-07-08 blanket-ReleaseFast decree.

## The rule

| Class | Default | Why |
|---|---|---|
| **Zortran** - has Fortran sources, or the build links Fortran | `ReleaseSafe` | Hot math is gfortran-compiled and unaffected by Zig optimize mode; the Zig layer is thin glue, so checks are cheap and convert host mistakes into clean panics instead of UB. |
| **Pure Zig** - neither | `ReleaseFast` | The compute IS the Zig, so safety checks are on the hot path. |

Ambiguous repos default to `ReleaseSafe`: being wrong toward safety is cheap, being wrong toward UB is not.

## This repo

- Fortran sources: **13**
- build.zig references to gfortran / Fortran objects: **5**
- Therefore: **ZORTRAN**, default **`ReleaseSafe`**

The hot math is gfortran-compiled at -O2/-O3 and is UNAFFECTED by Zig optimize mode. The Zig layer is thin validation and dispatch, so leaving bounds/overflow checks on costs little -- and it turns a bad pointer or length from the host into a clean panic at the C ABI boundary instead of undefined behaviour inside the consumer (Unreal).

## Overriding

`-Doptimize=` is still honoured. Use it for a measured hot path, not by habit:

```
zig build -Doptimize=ReleaseFast   # opt out of safety checks
zig build -Doptimize=Debug         # never for delivery
```

**Debug is never a delivery mode.** `b.standardOptimizeOption(.{})` with no baked default silently
yields Debug on a bare `zig build`; a fleet audit on 2026-08-04 found 15 of 52 modules doing exactly
that. Debug materializes `undefined` as real bytes and once shipped a 68 MB forCV archive.