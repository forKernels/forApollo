// forApollo wrap tool — the final delivery step (run by build.zig).
// Pure Zig, no shell scripts. Produces a wrapped libforapollo.a where ONLY the
// Zig-origin fq_ exports are public.
//
// forQuant's Fortran kernels historically used bind(C, name="fq_*") — the
// PUBLIC prefix — so a '*fq_*' glob cannot separate the Zig boundary (426
// Fortran-origin publics in the 2026-07-02 audit) from the internals. The
// keep list is therefore an EXPLICIT MANIFEST of the Zig export names
// (tools/keep_public.txt, regenerate with scripts/gen_keep_manifest.sh)
// passed as '@<path>'. Fail-closed. This also internalizes the bundled-copy
// sibling symbols (forla_/gmc_/trit_/ternary_) the audit flagged.
//
//   1. copy input objects into workdir             (inputs never mutated)
//   2. localize *_MOD_* + duplicate globals        (so ld -r can combine)
//   3. ld -r  [zig export object + Fortran .o] -> one relocatable
//   4. objcopy --keep-global-symbol <each entry>   (Fortran goes internal)
//   5. ar rcs libforapollo.a combined.o
//   6. GATE: nm the archive — any public defined symbol outside the keep
//      list FAILS the build. The contract is enforced, not assumed.
//
// Same tool pattern as forMath/forGeo/forApollo/forBayes/forOpt/forTernary/
// forNet/forIO tools/wrap.zig; the keep list is the only per-library variable.
//
// argv: wrap <out_lib> <workdir> <keep_globs(comma-sep)|@manifest-file> <obj>...
//   @file form: one glob per line, '#' comments and blank lines ignored.
const std = @import("std");
const builtin = @import("builtin");

// objcopy is `llvm-objcopy` on macOS (no GNU objcopy); GNU `objcopy` elsewhere.
// Override with the OBJCOPY env var.
fn objcopyBin(a: std.mem.Allocator) []const u8 {
    if (std.process.getEnvVarOwned(a, "OBJCOPY")) |v| return v else |_| {}
    return if (builtin.os.tag == .macos) "/opt/homebrew/bin/llvm-objcopy" else "objcopy";
}

fn run(a: std.mem.Allocator, argv: []const []const u8) !void {
    var child = std.process.Child.init(argv, a);
    child.stdout_behavior = .Ignore;
    child.stderr_behavior = .Inherit;
    const term = try child.spawnAndWait();
    if (term != .Exited or term.Exited != 0) {
        std.debug.print("wrap: command failed: {s}\n", .{argv[0]});
        return error.CommandFailed;
    }
}

fn capture(a: std.mem.Allocator, argv: []const []const u8) ![]u8 {
    const res = try std.process.Child.run(.{ .allocator = a, .argv = argv, .max_output_bytes = 64 * 1024 * 1024 });
    a.free(res.stderr);
    return res.stdout;
}

pub fn main() !void {
    var arena_state = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena_state.deinit();
    const a = arena_state.allocator();

    const args = try std.process.argsAlloc(a);
    if (args.len < 5) {
        std.debug.print("usage: wrap <out_lib> <workdir> <keep_globs(comma-sep)> <obj>...\n", .{});
        return error.Usage;
    }
    const out_lib = args[1];
    const workdir = args[2];
    const inputs = args[4..];
    const objcopy = objcopyBin(a);

    // Parse the keep list into globs (for objcopy) and their bare tokens (for
    // the gate check — '*' stripped). Either inline comma-separated globs, or
    // '@<path>' — a manifest file with one glob per line (# comments OK).
    var globs = std.ArrayList([]const u8){};
    var tokens = std.ArrayList([]const u8){};
    {
        var keep_text: []const u8 = args[3];
        var delim: u8 = ',';
        if (std.mem.startsWith(u8, keep_text, "@")) {
            keep_text = try std.fs.cwd().readFileAlloc(a, keep_text[1..], 16 * 1024 * 1024);
            delim = '\n';
        }
        var it = std.mem.tokenizeScalar(u8, keep_text, delim);
        while (it.next()) |raw| {
            const g = std.mem.trim(u8, raw, " \t\r");
            if (g.len == 0 or g[0] == '#') continue;
            try globs.append(a, g);
            try tokens.append(a, std.mem.trim(u8, g, "*"));
        }
    }
    if (globs.items.len == 0) {
        std.debug.print("wrap: empty keep list\n", .{});
        return error.Usage;
    }

    // ---- copy every input object into workdir (localize-symbol mutates the
    // object, and Stage-1 outputs under prebuilt/obj must never be touched) ----
    std.fs.cwd().makePath(workdir) catch {};
    var objs = std.ArrayList([]const u8){};
    for (inputs, 0..) |in, i| {
        const base = std.fs.path.basename(in);
        const dst = try std.fs.path.join(a, &.{ workdir, try std.fmt.allocPrint(a, "{d}__{s}", .{ i, base }) });
        try std.fs.cwd().copyFile(in, std.fs.cwd(), dst, .{});
        try objs.append(a, dst);
    }

    // ---- build symbol -> [owning objects] map (nm -g), find duplicate globals ----
    var sym_owners = std.StringHashMap(std.ArrayListUnmanaged(usize)).init(a);
    for (objs.items, 0..) |o, idx| {
        const out = try capture(a, &.{ "nm", "-g", o });
        var it = std.mem.tokenizeScalar(u8, out, '\n');
        while (it.next()) |line| {
            var toks = std.mem.tokenizeAny(u8, line, " \t\r");
            var f: [3][]const u8 = undefined;
            var n: usize = 0;
            while (toks.next()) |t| : (n += 1) {
                if (n < 3) f[n] = t;
            }
            if (n < 2) continue;
            // "addr TYPE name" (defined) or "TYPE name" (undefined)
            const typ = if (n >= 3) f[1] else f[0];
            const name = if (n >= 3) f[2] else f[1];
            if (typ.len != 1) continue;
            const c = typ[0];
            if (!std.ascii.isUpper(c) or c == 'U' or c == 'C') continue;
            const gop = try sym_owners.getOrPut(name);
            if (!gop.found_existing) gop.value_ptr.* = .{};
            try gop.value_ptr.append(a, idx);
        }
    }

    // ---- localize duplicated globals in all owners but the first, so ld -r
    // can combine. NOTE (org rule, 2026-07-02): unlike sibling repos,
    // do NOT pre-localize every *_MOD_* symbol — Zig boundaries in several repos
    // extern the gfortran-mangled ___fordsp_*_MOD_* names
    // DIRECTLY, so their definitions must stay global until ld -r has bound
    // the references. The post-combine `--keep-global-symbol` pass hides them
    // afterwards (they don't match the keep globs), which is safe because by then
    // every reference is already resolved inside the combined object. ----
    var loc_lists = try a.alloc(std.ArrayListUnmanaged([]const u8), objs.items.len);
    for (loc_lists) |*l| l.* = .{};
    var sit = sym_owners.iterator();
    while (sit.next()) |e| {
        const owners = e.value_ptr.items;
        if (owners.len < 2) continue;
        for (owners[1..]) |oi| try loc_lists[oi].append(a, e.key_ptr.*);
    }
    var localized: usize = 0;
    for (objs.items, 0..) |o, idx| {
        if (loc_lists[idx].items.len == 0) continue;
        var argv = std.ArrayList([]const u8){};
        try argv.append(a, objcopy);
        for (loc_lists[idx].items) |s| {
            try argv.append(a, "--localize-symbol");
            try argv.append(a, s);
        }
        try argv.append(a, o);
        try run(a, argv.items);
        localized += 1;
    }

    // ---- ld -r : all objects -> combined.o ----
    const combined = try std.fs.path.join(a, &.{ workdir, "forapollo_combined.o" });
    {
        var argv = std.ArrayList([]const u8){};
        try argv.appendSlice(a, &.{ "ld", "-r", "-o", combined });
        for (objs.items) |o| try argv.append(a, o);
        try run(a, argv.items);
    }

    // ---- keep ONLY the listed globs global; all Fortran goes internal ----
    // Leading '*' in each glob matches both `_<name>` (Mach-O) and `<name>` (ELF).
    {
        var argv = std.ArrayList([]const u8){};
        try argv.appendSlice(a, &.{ objcopy, "--wildcard" });
        for (globs.items) |g| {
            try argv.append(a, "--keep-global-symbol");
            try argv.append(a, g);
        }
        try argv.append(a, combined);
        try run(a, argv.items);
    }

    // ---- post-keep MOD sweep: suffix-anchored keep entries (e.g. the forGame
    // compat list's `*sph_init`) can accidentally match module-mangled names
    // ending in the same routine name (___x_MOD_sph_init). References are
    // already bound by ld -r, so localizing survivors here is always safe. ----
    {
        const pub_syms = try capture(a, &.{ "nm", "-g", combined });
        var argv = std.ArrayList([]const u8){};
        try argv.append(a, objcopy);
        var n_mod: usize = 0;
        var it = std.mem.tokenizeScalar(u8, pub_syms, '\n');
        while (it.next()) |line| {
            if (std.mem.indexOf(u8, line, "_MOD_") == null) continue;
            if (std.mem.indexOf(u8, line, " U ") != null) continue;
            var toks = std.mem.tokenizeAny(u8, line, " \t\r");
            var name: []const u8 = "";
            while (toks.next()) |t| name = t;
            try argv.append(a, "--localize-symbol");
            try argv.append(a, name);
            n_mod += 1;
        }
        if (n_mod > 0) {
            try argv.append(a, combined);
            try run(a, argv.items);
        }
    }

    // ---- ar rcs out_lib combined.o ----
    if (std.fs.path.dirname(out_lib)) |d| std.fs.cwd().makePath(d) catch {};
    std.fs.cwd().deleteFile(out_lib) catch {};
    try run(a, &.{ "ar", "rcs", out_lib, combined });

    // ---- GATE: every public defined symbol must match a keep token ----
    var leaks: usize = 0;
    var kept: usize = 0;
    {
        const out = try capture(a, &.{ "nm", "-g", out_lib });
        var it = std.mem.tokenizeScalar(u8, out, '\n');
        while (it.next()) |line| {
            var toks = std.mem.tokenizeAny(u8, line, " \t\r");
            var f: [3][]const u8 = undefined;
            var n: usize = 0;
            while (toks.next()) |t| : (n += 1) {
                if (n < 3) f[n] = t;
            }
            if (n < 3) continue; // undefined (no addr) or noise
            const typ = f[1];
            const name = f[2];
            if (typ.len != 1) continue;
            const c = typ[0];
            if (!std.ascii.isUpper(c) or c == 'U' or c == 'C') continue;
            // _MOD_ is never a real C-ABI name — flag regardless of token match.
            const is_mod = std.mem.indexOf(u8, name, "_MOD_") != null;
            var matched = false;
            if (!is_mod) for (tokens.items) |t| {
                if (std.mem.indexOf(u8, name, t) != null) {
                    matched = true;
                    break;
                }
            };
            if (matched) {
                kept += 1;
            } else {
                std.debug.print("wrap: GATE LEAK: {s} {s}\n", .{ typ, name });
                leaks += 1;
            }
        }
    }
    if (leaks != 0) {
        std.debug.print("wrap: GATE FAILED — {d} public symbols outside keep list\n", .{leaks});
        return error.GateFailed;
    }

    // ---- GATE 2: no SEVERED *_MOD_* bindings. An undefined module-mangled
    // symbol whose definition exists LOCALLY (as a localized `t`) in this same
    // archive means the definition was hidden before ld -r could bind the
    // reference — fails at every consumer link (found live 2026-07-02).
    // Undefined MOD refs with NO local definition are cross-repo mangled
    // bindings (pre-existing debt — siblings now hide their MOD symbols, so
    // these can never resolve): WARN loudly, don't fail. ----
    var severed: usize = 0;
    {
        const all_syms = try capture(a, &.{ "nm", out_lib });
        const pub_syms = try capture(a, &.{ "nm", "-g", out_lib });
        var it = std.mem.tokenizeScalar(u8, pub_syms, '\n');
        while (it.next()) |line| {
            if (std.mem.indexOf(u8, line, "_MOD_") == null) continue;
            if (std.mem.indexOf(u8, line, " U ") == null) continue;
            var toks = std.mem.tokenizeAny(u8, line, " \t\r");
            var name: []const u8 = "";
            while (toks.next()) |t| name = t;
            const defined_locally = blk: {
                var lit = std.mem.tokenizeScalar(u8, all_syms, '\n');
                while (lit.next()) |l| {
                    if (std.mem.indexOf(u8, l, " U ") != null) continue;
                    if (std.mem.endsWith(u8, l, name)) break :blk true;
                }
                break :blk false;
            };
            if (defined_locally) {
                std.debug.print("wrap: GATE SEVERED: {s}\n", .{name});
                severed += 1;
            } else {
                std.debug.print("wrap: WARNING cross-repo mangled binding (cannot resolve against wrapped siblings): {s}\n", .{name});
            }
        }
    }
    if (severed != 0) {
        std.debug.print("wrap: GATE FAILED — {d} severed *_MOD_* bindings (definition localized pre-link)\n", .{severed});
        return error.GateFailed;
    }

    std.debug.print("wrap: {s} <- {d} objects ({d} dup/MOD-localized); {d} public Zig exports, all Fortran internal, gate PASS\n", .{ out_lib, objs.items.len, localized, kept });
}
