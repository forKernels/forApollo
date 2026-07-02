// forApollo wrap tool — the final Stage-2 delivery step (run by build.zig).
// Pure Zig, no shell scripts. Produces a fully Zortran-wrapped libforapollo.a:
// every Fortran symbol (fa_* bind(C), mat_*/vec_* helpers, module-mangled
// *_MOD_*) goes INTERNAL; ONLY the Zig forapollo_* C-ABI exports stay public.
//
//   1. copy input objects into workdir             (inputs never mutated)
//   2. localize *_MOD_* + duplicate globals        (so ld -r can combine)
//   3. ld -r  [zig export object + Fortran .o] -> one relocatable
//   4. objcopy --keep-global-symbol '*forapollo_*' (Fortran goes internal)
//   5. ar rcs libforapollo.a combined.o
//   6. GATE: nm the archive — any public defined symbol outside the keep
//      glob FAILS the build. The contract is enforced, not assumed.
//
// Why not zig's addLibrary alone: it bundles raw .o with all globals exposed
// (the fa_* leak), and Zig 0.15.2's MachO linker segfaults on duplicate
// globals in other forKernels libs — ld -r (system linker) is the canon path
// (same tool pattern as forMath/forGeo tools/wrap.zig; keep_glob is the only
// per-library variable).
//
// argv: wrap <out_lib> <workdir> <keep_glob> <obj>...
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

// keep_glob is '*<token>*' — reduce to the inner token for the gate check.
fn globToken(glob: []const u8) []const u8 {
    return std.mem.trim(u8, glob, "*");
}

pub fn main() !void {
    var arena_state = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena_state.deinit();
    const a = arena_state.allocator();

    const args = try std.process.argsAlloc(a);
    if (args.len < 5) {
        std.debug.print("usage: wrap <out_lib> <workdir> <keep_glob> <obj>...\n", .{});
        return error.Usage;
    }
    const out_lib = args[1];
    const workdir = args[2];
    const keep_glob = args[3];
    const inputs = args[4..];
    const objcopy = objcopyBin(a);

    // ---- copy every input object into workdir (localize-symbol mutates the
    // object, and Stage-1 outputs under build/obj must never be touched) ----
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

    // ---- localize: (a) every gfortran module-mangled symbol (`*_MOD_*` — never
    // a real C-ABI name, and a module named with the lib prefix would otherwise
    // survive the keep-glob and leak), in ALL owners; (b) each duplicated global
    // in all owners but the first. ld -r then resolves every now-local symbol
    // inside the single combined object. ----
    var loc_lists = try a.alloc(std.ArrayListUnmanaged([]const u8), objs.items.len);
    for (loc_lists) |*l| l.* = .{};
    var sit = sym_owners.iterator();
    while (sit.next()) |e| {
        const owners = e.value_ptr.items;
        const is_mod = std.mem.indexOf(u8, e.key_ptr.*, "_MOD_") != null;
        if (owners.len < 2 and !is_mod) continue;
        const start: usize = if (is_mod) 0 else 1; // _MOD: localize in ALL owners
        for (owners[start..]) |oi| try loc_lists[oi].append(a, e.key_ptr.*);
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

    // ---- keep ONLY the lib's prefix global; all Fortran goes internal ----
    // keep_glob (e.g. '*forapollo_*') matches both `_<prefix>` (Mach-O leading _)
    // and `<prefix>` (ELF).
    try run(a, &.{ objcopy, "--wildcard", "--keep-global-symbol", keep_glob, combined });

    // ---- ar rcs out_lib combined.o ----
    if (std.fs.path.dirname(out_lib)) |d| std.fs.cwd().makePath(d) catch {};
    std.fs.cwd().deleteFile(out_lib) catch {};
    try run(a, &.{ "ar", "rcs", out_lib, combined });

    // ---- GATE: every public defined symbol must match the keep glob ----
    const token = globToken(keep_glob);
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
            if (std.mem.indexOf(u8, name, token) != null) {
                kept += 1;
            } else {
                std.debug.print("wrap: GATE LEAK: {s} {s}\n", .{ typ, name });
                leaks += 1;
            }
        }
    }
    if (leaks != 0) {
        std.debug.print("wrap: GATE FAILED — {d} public symbols outside {s}\n", .{ leaks, keep_glob });
        return error.GateFailed;
    }

    std.debug.print("wrap: {s} <- {d} objects ({d} dup/MOD-localized); {d} public {s}, all Fortran internal, gate PASS\n", .{ out_lib, objs.items.len, localized, kept, keep_glob });
}
