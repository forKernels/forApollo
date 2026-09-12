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
// ===========================================================================
// Zig 0.15.2 / 0.16.0 compatibility layer.
//
// This tool SEALS the delivery archive, so the moment a target moves to 0.16
// and this file does not, that target stops shipping a sealed library at all.
// The build script configuring cleanly says nothing about it: `zig build -l`
// enumerates steps without ever compiling this file.
//
// Every divergence below is expressed by ALIASING one of two functions rather
// than branching inside one. An unreferenced function is never analysed, so
// 0.16 never sees argsAlloc and 0.15.2 never sees std.process.Init -- which a
// comptime branch inside a single body cannot achieve, because the two mains
// need different SIGNATURES.
//
// Taken from forNet/forIO/forBio tools/wrap.zig, which already carry it. Kept
// deliberately identical rather than improved, so one fix travels the fleet
// instead of twenty-five variants drifting. Duplicated rather than imported
// because a build tool must not acquire a cross-repo dependency to run.
// ===========================================================================
const zig16 = @import("builtin").zig_version.order(.{ .major = 0, .minor = 16, .patch = 0 }) != .lt;

fn objcopyBin(override: ?[]const u8) []const u8 {
    if (override) |v| return v;
    return if (builtin.os.tag == .macos) "/opt/homebrew/bin/llvm-objcopy" else "objcopy";
}

// 0.16 replaced std.process.Child with std.process.spawn.
const runRaw = if (zig16) run16 else run15;

// ---- command-line length guard -------------------------------------------
//
// Windows caps a process command line at 32767 characters. wrap passes one
// `--localize-symbol <name>` / `--keep-global-symbol <name>` PAIR per symbol,
// so a long `_MOD_` list or a large keep manifest overflows it and
// CreateProcessW fails with NameTooLong -- the delivery build then produces
// NO ARCHIVE AT ALL. Measured: forCV broke when its manifest grew 599 -> 921
// (argv ~43KB); forMath's largest per-object localize list is 25128 bytes
// against the 32767 cap. Invisible on Linux/macOS (ARG_MAX ~1MB).
//
// GNU binutils and llvm (ld / objcopy / ar) all read tail args from a
// `@response` file, so an over-long command spills args[1..] to one. Short
// commands stay inline, so an already-built repo is byte-for-byte unchanged.
//
// The response file goes in wrap's per-invocation workdir, not the CWD: a
// repo that builds many packs may run them in parallel and a fixed CWD name
// would have them clobber each other.
var g_wrap_workdir: []const u8 = "";

inline fn wrapWriteFile(sub_path: []const u8, data: []const u8) !void {
    return if (zig16)
        std.Io.Dir.cwd().writeFile(cwdIo(), .{ .sub_path = sub_path, .data = data })
    else
        std.fs.cwd().writeFile(.{ .sub_path = sub_path, .data = data });
}

fn run(a: std.mem.Allocator, argv: []const []const u8) !void {
    var total: usize = 0;
    for (argv) |x| total += x.len + 3;
    if (total <= 24000 or argv.len <= 2) return runRaw(a, argv);

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    for (argv[1..]) |x| {
        try buf.append(a, '"');
        try buf.appendSlice(a, x);
        try buf.appendSlice(a, "\"\n");
    }
    const dir = if (g_wrap_workdir.len != 0) g_wrap_workdir else ".";
    const rsp = try std.fs.path.join(a, &.{ dir, "wrap_args.rsp" });
    try wrapWriteFile(rsp, buf.items);
    const at = try std.fmt.allocPrint(a, "@{s}", .{rsp});
    try runRaw(a, &.{ argv[0], at });
}

fn run15(a: std.mem.Allocator, argv: []const []const u8) !void {
    var child = std.process.Child.init(argv, a);
    child.stdout_behavior = .Ignore;
    child.stderr_behavior = .Inherit;
    return checkTerm(try child.spawnAndWait(), argv);
}

fn run16(a: std.mem.Allocator, argv: []const []const u8) !void {
    var threaded: std.Io.Threaded = .init(a, .{});
    defer threaded.deinit();
    const tio = threaded.io();
    var child = try std.process.spawn(tio, .{ .argv = argv, .stdout = .ignore });
    return checkTerm(try child.wait(tio), argv);
}

/// 0.16 lowercased the Term union tags: .Exited -> .exited. `anytype` lets one
/// body serve both; the tag is resolved per toolchain at comptime.
fn checkTerm(term: anytype, argv: []const []const u8) !void {
    const ok = if (comptime @hasField(@TypeOf(term), "Exited"))
        term == .Exited and term.Exited == 0
    else
        term == .exited and term.exited == 0;
    if (!ok) {
        std.debug.print("wrap: command failed: {s}\n", .{argv[0]});
        return error.CommandFailed;
    }
}

const capture = if (zig16) capture16 else capture15;

fn capture15(a: std.mem.Allocator, argv: []const []const u8) ![]u8 {
    const res = try std.process.Child.run(.{ .allocator = a, .argv = argv, .max_output_bytes = 64 * 1024 * 1024 });
    a.free(res.stderr);
    return res.stdout;
}

/// 0.16 has no run-and-capture, so spawn with stdout redirected to a temp file
/// and read it back. The name is salted with the argv pointer so two concurrent
/// wrap processes cannot collide on it.
fn capture16(a: std.mem.Allocator, argv: []const []const u8) ![]u8 {
    const tmp = try std.fmt.allocPrint(a, ".forapollo-capture-{x}.tmp", .{@intFromPtr(argv.ptr)});
    defer deleteFile(tmp) catch {};
    {
        var threaded: std.Io.Threaded = .init(a, .{});
        defer threaded.deinit();
        const tio = threaded.io();
        const out_file = try std.Io.Dir.cwd().createFile(tio, tmp, .{ .read = true });
        defer out_file.close(tio);
        var child = try std.process.spawn(tio, .{ .argv = argv, .stdout = .{ .file = out_file } });
        _ = try child.wait(tio);
    }
    return readFileAlloc(a, tmp, 64 * 1024 * 1024);
}

// 0.16 removed std.fs.cwd(); directory handles live on std.Io.Dir and every
// file call takes an `io`.
inline fn cwdIo() if (zig16) std.Io else void {
    if (zig16) return std.Io.Threaded.global_single_threaded.io();
    return {};
}

/// 0.16 RENAMED makePath to createDirPath. The old name still compiles under
/// 0.15.2, so a plain rename would have broken the older toolchain silently.
inline fn makePath(path: []const u8) !void {
    return if (zig16)
        std.Io.Dir.cwd().createDirPath(cwdIo(), path)
    else
        std.fs.cwd().makePath(path);
}

inline fn deleteFile(path: []const u8) !void {
    return if (zig16)
        std.Io.Dir.cwd().deleteFile(cwdIo(), path)
    else
        std.fs.cwd().deleteFile(path);
}

/// 0.16 INSERTS an `io` parameter before `options`. The four dir/path
/// parameters keep their order on both toolchains:
///
///   0.15.2  copyFile(source_dir, source_path, dest_dir, dest_path, options)
///   0.16.0  copyFile(source_dir, source_path, dest_dir, dest_path, io, options)
///
/// An earlier version of this comment claimed the ORDER changed and that the
/// wrong order would still compile. Both halves were false, and the claim was
/// copied into twenty of these files before the macOS agent checked it against
/// the two signatures rather than trusting it. The code was always right; only
/// the comment lied.
inline fn copyFile(src: []const u8, dst: []const u8) !void {
    return if (zig16)
        std.Io.Dir.cwd().copyFile(src, std.Io.Dir.cwd(), dst, cwdIo(), .{})
    else
        std.fs.cwd().copyFile(src, std.fs.cwd(), dst, .{});
}

inline fn readFileAlloc(gpa: std.mem.Allocator, path: []const u8, max: usize) ![]u8 {
    return if (zig16)
        std.Io.Dir.cwd().readFileAlloc(cwdIo(), path, gpa, std.Io.Limit.limited(max))
    else
        std.fs.cwd().readFileAlloc(gpa, path, max);
}

pub const main = if (zig16) main16 else main15;

fn main15() !void {
    var arena_state = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena_state.deinit();
    const a = arena_state.allocator();
    const args = try std.process.argsAlloc(a);
    const objcopy_override = std.process.getEnvVarOwned(a, "OBJCOPY") catch null;
    return wrapMain(a, args, objcopy_override);
}

fn main16(init: std.process.Init.Minimal) !void {
    var arena_state = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena_state.deinit();
    const a = arena_state.allocator();
    var it = try std.process.Args.Iterator.initAllocator(init.args, a);
    defer it.deinit();
    var list: std.ArrayListUnmanaged([]const u8) = .empty;
    while (it.next()) |arg| try list.append(a, arg);
    // getAlloc, not getPosix: getPosix is POSIX-only and will not compile for
    // a Windows target.
    const objcopy_override = init.environ.getAlloc(a, "OBJCOPY") catch null;
    return wrapMain(a, list.items, objcopy_override);
}



fn wrapMain(a: std.mem.Allocator, args: []const []const u8, objcopy_override: ?[]const u8) !void {
    if (args.len < 5) {
        std.debug.print("usage: wrap <out_lib> <workdir> <keep_globs(comma-sep)> <obj>...\n", .{});
        return error.Usage;
    }
    const out_lib = args[1];
    const workdir = args[2];
    g_wrap_workdir = workdir;
    const inputs = args[4..];
    const objcopy = objcopyBin(objcopy_override);

    // Parse the keep list into globs (for objcopy) and their bare tokens (for
    // the gate check — '*' stripped). Either inline comma-separated globs, or
    // '@<path>' — a manifest file with one glob per line (# comments OK).
    var globs: std.ArrayListUnmanaged([]const u8) = .empty;
    var tokens: std.ArrayListUnmanaged([]const u8) = .empty;
    {
        var keep_text: []const u8 = args[3];
        var delim: u8 = ',';
        if (std.mem.startsWith(u8, keep_text, "@")) {
            keep_text = try readFileAlloc(a, keep_text[1..], 16 * 1024 * 1024);
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
    makePath(workdir) catch {};
    var objs: std.ArrayListUnmanaged([]const u8) = .empty;
    for (inputs, 0..) |in, i| {
        const base = std.fs.path.basename(in);
        const dst = try std.fs.path.join(a, &.{ workdir, try std.fmt.allocPrint(a, "{d}__{s}", .{ i, base }) });
        try copyFile(in, dst);
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
            if (!gop.found_existing) gop.value_ptr.* = .empty;
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
    for (loc_lists) |*l| l.* = .empty;
    var sit = sym_owners.iterator();
    while (sit.next()) |e| {
        const owners = e.value_ptr.items;
        if (owners.len < 2) continue;
        for (owners[1..]) |oi| try loc_lists[oi].append(a, e.key_ptr.*);
    }
    var localized: usize = 0;
    for (objs.items, 0..) |o, idx| {
        if (loc_lists[idx].items.len == 0) continue;
        var argv: std.ArrayListUnmanaged([]const u8) = .empty;
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
        var argv: std.ArrayListUnmanaged([]const u8) = .empty;
        try argv.appendSlice(a, &.{ "ld", "-r", "-o", combined });
        for (objs.items) |o| try argv.append(a, o);
        try run(a, argv.items);
    }

    // ---- keep ONLY the listed globs global; all Fortran goes internal ----
    // Leading '*' in each glob matches both `_<name>` (Mach-O) and `<name>` (ELF).
    {
        var argv: std.ArrayListUnmanaged([]const u8) = .empty;
        try argv.appendSlice(a, &.{ objcopy, "--wildcard" });
        // LINKAGE MACHINERY — must stay GLOBAL, and this is not optional.
        // On MinGW/PE, GCC and Zig reach an imported global through an
        // indirection stub `.refptr.<sym>`; for a stack-protected function that
        // is `.refptr.__stack_chk_guard`, loaded RIP-relative in the PROLOGUE.
        // It is emitted per-object and must stay global so it COMDAT-merges to
        // one copy. Localized, the linker leaves the displacement at ZERO and
        // the function dereferences its own opcode bytes -> 0xC0000005, in any
        // gfortran-linked consumer. Toolchain stubs, not API names.
        //   nm <lib> | grep -w '.refptr.__stack_chk_guard'   R = fine, r = BROKEN
        try argv.appendSlice(a, &.{ "--keep-global-symbol", ".refptr.*" });
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
        var argv: std.ArrayListUnmanaged([]const u8) = .empty;
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
    if (std.fs.path.dirname(out_lib)) |d| makePath(d) catch {};
    deleteFile(out_lib) catch {};
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
            // COFF `.refptr.<sym>` stubs are toolchain linkage machinery, not
            // API. The keep pass deliberately holds them global so the stack
            // protector's RIP-relative prologue load resolves; localizing one
            // strands the displacement at 0 and the function faults on its own
            // opcode bytes. A stub carries no callable surface, so this opens
            // no bind-Fortran path past the seal.
            var matched = std.mem.indexOf(u8, name, ".refptr.") != null;
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
                    // Last token, not endsWith on the raw line: nm pipes CRLF on
                    // winX86, so every line ends "<name>\r" and endsWith(name)
                    // is false for EVERY symbol. defined_locally was therefore
                    // always false, silently DOWNGRADING this gate from a failure
                    // to the cross-repo WARNING -- a real severed binding would
                    // have shipped with a warning instead of stopping the build.
                    // Measured in forCV: 38488 bytes of nm output, 923 CR.
                    var ltoks = std.mem.tokenizeAny(u8, l, " \t\r");
                    var llast: []const u8 = "";
                    while (ltoks.next()) |t| llast = t;
                    if (std.mem.eql(u8, llast, name)) break :blk true;
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
