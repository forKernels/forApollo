# linX86 — x86_64-linux-gnu target branch

This branch carries patches specific to the x86_64-linux-gnu target
(Ubuntu 24.04 noble, Zig 0.15.2, gfortran 14.3.0).

## Reference machine
- 2013 MacBook Pro, Intel x86_64, Ubuntu 24.04
- Zig 0.15.2 at ~/.local/bin/zig
- gfortran 14.3.0 at /usr/bin/gfortran-14
- libgfortran.so at /usr/lib/gcc/x86_64-linux-gnu/14/

## Branch relationships
- main        — target-agnostic truth
- jetson-thor — aarch64 Linux + Jetson AGX Thor GPU (nvfortran -gpu=cc110)
- linX86      — this branch
- winX86      — Windows x86_64 (MSYS2 UCRT64)

## linX86 changes typical of this module
- FC env var plumbing where scripts hardcoded bare gfortran
- /usr/lib/gcc/x86_64-linux-gnu/14 in link search paths
- x86_64-linux-gnu in per-module TargetId enums where missing
- forCUDA sibling path resolution extended to linux-x86_64

## Delivery build (2026-11-05)
The delivered `prebuilt/linX86/` archive is built in Docker
(`forkernels-linx86:zig016`: Ubuntu 24.04, Zig 0.16.0, gfortran 14.2.0) with
`zig build -Dtarget=x86_64-linux-gnu.2.34 -Dcpu=x86_64_v3`: glibc floor 2.34,
no AVX-512, no libmvec (`_ZGV*`), PIC. Tests run at native.
