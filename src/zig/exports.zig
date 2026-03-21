// forApollo — Zig C ABI Exports
// Copyright The Fantastic Planet — By David Clabaugh
//
// All functions prefixed forapollo_* and call fa_* Fortran kernels
// through the safety layer. Python/Rust/C++ call these.

pub export fn forapollo_version() callconv(.c) u32 {
    return 0x000100; // 0.1.0
}
