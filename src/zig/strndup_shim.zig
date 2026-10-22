// =============================================================================
// strndup for Windows test executables
// =============================================================================
// libgfortran's string.o calls strndup, which Zig's bundled mingw-w64 lacks.
// Linking MSYS2's whole libmingwex for it drags in a conflicting _assert.
// build.zig compiles this as its own object and links it into the Windows test
// exe only; no archive links it.

const std = @import("std");

export fn strndup(s: [*:0]const u8, n: usize) callconv(.c) ?[*:0]u8 {
    const len = std.mem.indexOfScalar(u8, s[0..n], 0) orelse n;
    const out: [*]u8 = @ptrCast(std.c.malloc(len + 1) orelse return null);
    @memcpy(out[0..len], s[0..len]);
    out[len] = 0;
    return out[0..len :0];
}
