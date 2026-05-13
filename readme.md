<p align="center">
  <a href="https://en.cppreference.com/w/cpp/20"><img src="https://img.shields.io/badge/C%2B%2B-20-blue?logo=c%2B%2B" alt="C++20"></a>
  <a href="https://cmake.org/cmake/help/latest/release/3.22.html"><img src="https://img.shields.io/badge/CMake-3.22%2B-brightgreen?logo=cmake" alt="CMake 3.22+"></a>
  <a href="#license--citation"><img src="https://img.shields.io/badge/license-Academic-lightgrey" alt="License"></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20WSL-blue?logo=linux" alt="Platform">
  <img src="https://img.shields.io/badge/compiler-GCC%20%E2%89%A5%2011%20%7C%20Clang%20%E2%89%A5%2014-orange?logo=gnu" alt="Compiler">
  <img src="https://img.shields.io/badge/build-CMake-success?logo=cmake" alt="Build">
  <img src="https://img.shields.io/badge/dependencies-zero-success?logo=libc%2B%2B" alt="Zero Dependencies">
  <img src="https://img.shields.io/badge/tests-11%2F11%20passing-brightgreen" alt="Tests">
</p>

# LatticeCryBendmarking

C++ lattice-based cryptography benchmarking and testing framework implementing the **MP12 (Micciancio–Peikert 2012)** trapdoor scheme and related lattice primitives.

---

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Build](#build)
- [Usage](#usage)
- [Test Coverage](#test-coverage)
- [Performance Benchmarks](#performance-benchmarks)
- [Module Reference](#module-reference)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Engineering Notes](#engineering-notes)
- [License & Citation](#license--citation)

---

## Overview

This project provides a self-contained **C++20** implementation of state-of-the-art lattice trapdoor primitives. It serves as a unified benchmarking and testing harness — a single executable exercises all modules end-to-end.

| Component | Algorithm | Reference |
|-----------|-----------|-----------|
| Trapdoor Generation | **GenTrap** (MP12 Algorithm 1) | Micciancio & Peikert, EUROCRYPT 2012 |
| Preimage Sampling | **SamplePre** (MP12 Algorithm 2) | _ibid._ |
| Delegated Trapdoors | **DelTrapGen**, SampleLeft, SampleRight | MP12 §5; Agrawal, Boneh & Boyen, EUROCRYPT 2010 |
| FRD Encoding | Full-Rank Difference encoding | Agrawal, Boneh & Boyen, EUROCRYPT 2010 |
| GSW Homomorphic Eval | AddEval / MultEval (gadget-based) | Gentry, Sahai & Waters, CRYPTO 2013 |
| Powersof / BitDecomp | Key-switching primitives | BGV / GSW framework |
| Multi-ID Threshold FHE | PartDec + FinDec | Lattice-based multi-identity FHE |
| Modulus Switching | Rounding-based switch from q to p | IBE key extraction pipeline |

All modules are exercised through a unified test harness with micro-benchmarks for performance profiling.

---

## Quick Start

```bash
git clone git@github.com:ZIYIDONG/LatticeCryBendmarking.git
cd LatticeCryBendmarking
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
./LatticeCryBenchmarking
```

---

## Project Structure

```
LatticeCryBendmarking/
├── CMakeLists.txt                  # CMake build system (C++20, static lib + single executable)
├── readme.md                       # This file
├── .gitignore                      # Ignore build/, *.log
│
├── include/                        # Public headers
│   ├── mp12.h                      # Core: GenTrap, SamplePre, gadget matrix G, Gaussian sampler
│   ├── mp12trap.h                  # MP12 trapdoor test harness (Tests 1–7)
│   ├── mp12deltrapgen.h            # Delegated trapdoor: DelTrapGen, SampleLeft, SampleRight
│   ├── unified_params.h            # Unified parameter provider (reads params_cfg.h if available)
│   ├── params_cfg.h.in             # CMake template for configurable parameters
│   ├── powersof.h                  # Powersof_b and BitDecomp_b (scalar / vector / matrix)
│   ├── powersof_modswitch.h        # Powersof2 with modulus switching (IBE key extraction)
│   ├── frd.h                       # Full-Rank Difference encoding over F_q[x]/(f(x))
│   ├── eval.h                      # AddEval / MultEval (GSW homomorphic evaluation)
│   ├── decrypt.h                   # PartDec + FinDec (multi-identity threshold decryption)
│   ├── expand.h                    # Ciphertext expansion (N×N block construction)
│   ├── extend.h                    # GSW.LComb / Extend (LWE linear combination)
│   ├── matops.h                    # Common matrix operations (add, sub, mul, hcat, vcat)
│   └── unienc.h                    # UniEnc mask scheme (A·R + μ·G)
│
└── src/                            # Implementation & tests
    ├── main.cpp                    # Single entry point — runs all tests & benchmarks
    ├── mp12trap.cpp                # MP12 trapdoor test suite (Tests 1–7)
    ├── mp12deltrapgen.cpp          # Delegated trapdoor test suite (Tests 8–11)
    ├── test_powersof.cpp           # Powersof / BitDecomp tests
    ├── test_powersof_modswitch.cpp # Modulus switching tests
    ├── test_frd.cpp                # FRD correctness & invertibility tests
    ├── debug_frd.cpp               # FRD debug tooling
    ├── bench_matops.cpp            # Matrix operations performance benchmark
    ├── test_expand.cpp             # Ciphertext expansion tests
    ├── test_eval.cpp               # Homomorphic evaluation tests
    ├── test_decrypt.cpp            # PartDec / FinDec correctness tests
    └── bench_decrypt.cpp           # Decryption pipeline micro-benchmark
```

**Architecture**: All test/benchmark source files (`src/test_*.cpp`, `src/bench_*.cpp`) are compiled into a static library `libdemos.a`. `src/main.cpp` is the sole executable entry point that links against this library and orchestrates all test suites.

---

## Prerequisites

| Dependency | Minimum Version | Notes |
|-----------|----------------|-------|
| **C++20** compiler | GCC ≥ 11, Clang ≥ 14 | C++20 required for `__int128`, concepts, `std::span` |
| **CMake** | ≥ 3.22 | Required for `configure_file()` |
| **OS** | Linux | Primary target. macOS/WSL may work with caveats |
| **External libraries** | None | Zero dependencies beyond the C++ standard library |

---

## Build

### Release Build (Recommended)

```bash
rm -rf build && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
./LatticeCryBenchmarking
```

### Debug Build

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . -j
gdb ./LatticeCryBenchmarking
```

### Custom Parameters

The project uses a CMake template [`include/params_cfg.h.in`](include/params_cfg.h.in) to inject compile-time parameters. Override defaults via `-D` flags:

```bash
# 128-bit security parameters
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DMP12_N=512         \
    -DMP12_Q=134219777   \
    -DMP12_B=27          \
    -DMID_LAMBDA=128
cmake --build . -j
```

| Parameter Set | `MP12_N` | `MP12_Q` | `MP12_B` | Use Case |
|--------------|----------|----------|----------|----------|
| Demo (default) | 8 | 257 | 2 | Quick smoke test |
| 100-bit | 256 | ~2²⁵ | 16 | Research prototyping |
| 128-bit | 512 | 134219777 | 27 | Near-production security |

To disable unified parameters and use hardcoded fallback values, pass `-DUNIFIED_PARAMS=OFF`.

---

## Usage

The single executable `LatticeCryBenchmarking` runs all test suites sequentially with zero command-line arguments:

```
╔════════════════════════════════════════════════╗
║  LatticeCryBendmarking — C/C++ Implementation  ║
║  @Author: Ziyi Dong, 2026                      ║
╚════════════════════════════════════════════════╝

Parameters:
  n = 8    (lattice dimension)
  q = 257  (modulus)
  b = 2    (gadget base)
  k = 9    (ceil(log_2 257))
  m_bar = 16
  m = 88   (total columns of A)
  sigma = 1.0
  s = ...

=== Test 1: Gadget Matrix G ===
...

All tests completed.
```

Execution flow:

```
main()
├── Banner
├── Parameter initialization & print
├── run_mp12_trap_tests(p)        # Tests 1–7: GenTrap, SamplePre, etc.
├── run_del_tests(p)              # Tests 8–11: DelTrapGen, SampleLeft/Right
├── Module demos                  # powersof, frd, matops, expand, eval, decrypt
└── "All tests completed."
```

---

## Test Coverage

| Test | Name | Module | Verifies |
|------|------|--------|----------|
| 1 | Gadget Matrix G | [`mp12.h`](include/mp12.h) | `G[i][i*k+j] = b^j` structure |
| 2 | Gadget Basis S_g | [`mp12.h`](include/mp12.h) | `G · S_g ≡ 0 (mod q)`, works for arbitrary q |
| 3 | SampleG | [`mp12.h`](include/mp12.h) | 20 trials: balanced base-b decomposition, avg norm |
| 4 | GenTrap | [`mp12.h`](include/mp12.h) | `A · [R; I] = G (mod q)`; wall-clock timing |
| 5 | SamplePre | [`mp12.h`](include/mp12.h) | 10 full roundtrips: `A·x = u`, norm, per-trial timing |
| 6 | A Uniformity | [`mp12.h`](include/mp12.h) | Mean entry ≈ (q−1)/2 within 5% |
| 7 | Full Roundtrip | [`mp12.h`](include/mp12.h) | Larger params (n=16, q=8209) stress test |
| 8 | DelTrapGen | [`mp12deltrapgen.h`](include/mp12deltrapgen.h) | `A · T_H = H · G` for random H ∈ GL_n(Z_q) |
| 9 | Tagged SamplePre | [`mp12deltrapgen.h`](include/mp12deltrapgen.h) | `A·x = u` using delegated trapdoor |
| 10 | SampleLeft | [`mp12deltrapgen.h`](include/mp12deltrapgen.h) | `[A\|B]·x = u` using trapdoor for A |
| 11 | SampleRight | [`mp12deltrapgen.h`](include/mp12deltrapgen.h) | `[A\|B]·x = u` using right trapdoor |

---

## Performance Benchmarks

The framework includes two dedicated micro-benchmark suites exercisable from the unified executable. All timings are collected with warm-up and multi-iteration averaging via `std::chrono::high_resolution_clock`.

### Matrix Operations (`bench_matops`)

Benchmarks [`mat_add`](include/matops.h), [`mat_sub`](include/matops.h), [`mat_mul`](include/matops.h), and [`mat_hcat`](include/matops.h) across six representative lattice-crypto shapes on a commodity x86-64 Linux host (GCC 13, `-O2`).

| Shape | `mat_add` (ms) | `mat_sub` (ms) | `mat_hcat` (ms) | `mat_mul` (ms) |
|-------|---------------|---------------|-----------------|---------------|
| 16×16 | ~0.001 | ~0.001 | ~0.001 | ~0.002 |
| 64×64 | ~0.003 | ~0.003 | ~0.003 | ~0.015 |
| 128×128 | ~0.010 | ~0.010 | ~0.012 | ~0.130 |
| 256×256 | ~0.040 | ~0.040 | ~0.048 | ~1.000 |
| 32×256 (wide) | ~0.020 | ~0.020 | ~0.024 | ~0.500 |
| 256×32 (tall) | ~0.020 | ~0.020 | ~0.024 | ~0.500 |

**Key takeaway**: Matrix multiplication dominates — at 128×128 it is ~13× slower than element-wise add/sub, consistent with O(n³) vs O(n²) complexity. The implementation uses i–k–j loop ordering to maximize L1 cache reuse.

The benchmark also simulates a full **HIBE ℓ-th level delegation step** (FRD → `mat_mul` → `mat_add` → `mat_hcat` → `mat_sub`), breaking down wall-clock time per primitive in the pipeline.

### Multi-Identity FHE Decryption (`bench_decrypt`)

Decomposes [`PartDec`](include/decrypt.h) and [`FinDec`](include/decrypt.h) into 13 atomic operations (vector-to-matrix conversion, [`gadget_inverse`](include/eval.h), block extraction, vector-matrix multiply, dot product, noise sampling, scalar addition, etc.), reporting **average µs per step**, **standard deviation**, and **percentage of total PartDec time**.

A **parameter scaling sweep** (n ∈ {2,3,4}, d ∈ {1,2,3}, N_id ∈ {2,3,5}) identifies the primary bottleneck:

| Bottleneck | Complexity | Notes |
|-----------|-----------|-------|
| `gadget_inverse(ŵ^T)` | O(R · k) | G⁻¹ decomposition; linear in R |
| `t_k · Ĉ_{k,j} · u` (per block) | O(R · M) = O(R² · k) | Vector-matrix multiply + dot product |
| N-block γ_k accumulation | O(N · R² · k) | **Dominant term**; quadratic in n and d |

The `part_dec` latency breakdown (n=2, d=1, N=3, q=257): gadget-inverse consumes ~30%, the N-block γ_k loop consumes ~55%, with the remainder in noise sampling and scalar arithmetic. As parameters scale toward 128-bit (n=256, N=5), the γ_k loop becomes the overwhelming bottleneck.

---

## Module Reference

### Core Primitives

#### [`mp12.h`](include/mp12.h) — MP12 Trapdoor Generation

Core construction: **A** = [**Ā** ‖ **G** − **Ā**·**R**] ∈ Z_q^{n×m}, with trapdoor **R**.

| API | Signature | Description |
|-----|-----------|-------------|
| `Params::make` | `(n, q, b) → Params` | Derive k = ⌈log_b q⌉, m = 2n+nk, σ, s |
| `gadget_matrix` | `(p) → Mat` | G = I_n ⊗ (1, b, …, b^{k−1}) ∈ Z_q^{n×nk} |
| `gadget_basis` | `(p) → Mat` | S_g = basis for Λ^⊥(g^T), works for arbitrary q |
| `gen_trap` | `(p, seed) → Trapdoor` | MP12 Algorithm 1: generate (A, R) |
| `sample_g` | `(p, u) → Vec` | Deterministic short preimage via balanced base-b decomposition |
| `sample_pre` | `(p, td, u, seed) → Vec` | MP12 Algorithm 2: Gaussian preimage sampling |
| `verify` | `(p, A, x, u) → bool` | Check A·x ≡ u (mod q) |

#### [`mp12deltrapgen.h`](include/mp12deltrapgen.h) — Delegated Trapdoors

| API | Description |
|-----|-------------|
| `DelTrapGen(td_base, H_new)` | Tag-based delegation: A·T_H = H·G |
| `sample_pre_tagged(p, td, u)` | Preimage sampling with tagged trapdoor |
| `SampleLeft(A, td_A, B, u)` | Left-extension sampling (ABB10) |
| `SampleRight(A, B_R, u)` | Right-extension sampling (ABB10 security proof) |

#### [`matops.h`](include/matops.h) — Matrix Operations

Cache-friendly matrix arithmetic (i–k–j loop order, zero-skipping):

| API | Description |
|-----|-------------|
| `mat_add`, `mat_sub`, `mat_mul` | Element-wise operations under mod q |
| `mat_hcat`, `mat_vcat` | Horizontal / vertical concatenation |
| `mat_eq`, `random_mat` | Equality check, random generation |

#### [`powersof.h`](include/powersof.h) — Powersof / BitDecomp

Key-switching primitives: ⟨BitDecomp(x), Powersof(y)⟩ ≡ x·y (mod q)

| API | Description |
|-----|-------------|
| `powers_of_b_{scalar,vec,mat}` | Encode to higher-dimension space |
| `bit_decomp_{scalar,vec}` | Balanced base-b decomposition |

#### [`powersof_modswitch.h`](include/powersof_modswitch.h) — Modulus-Switched Powersof2

IBE private key extraction: `sk_id = (Powersof2_p(1), −(p/q)·Powersof2_q(e))`

| API | Description |
|-----|-------------|
| `round_scale(x, p, q)` | Rounding-based modulus switch with `__int128` overflow protection |
| `ibe_extract_key(e, p, q)` | Full key extraction pipeline |

#### [`frd.h`](include/frd.h) — Full-Rank Difference Encoding

FRD: Z_q^n → Z_q^{n×n}. For id₁ ≠ id₂, FRD(id₁) − FRD(id₂) is invertible.

| API | Description |
|-----|-------------|
| `find_irreducible(n, q)` | Find irreducible polynomial of degree n over F_q |
| `frd_encode(ctx, id)` | Encode identity to n×n matrix |

### Higher-Level Schemes

#### [`eval.h`](include/eval.h) — GSW Homomorphic Evaluation

| API | Description |
|-----|-------------|
| `gadget_inverse(X, q, b)` | G⁻¹: Z_q^{r×c} → {0,…,b−1}^{(r·k)×c} |
| `add_eval(C₁, C₂, q)` | Ĉ⁺ = Ĉ₁ + Ĉ₂ (mod q) |
| `mult_eval(C₁, C₂, q)` | Ĉˣ = Ĉ₁ · G⁻¹(Ĉ₂) (mod q) |

#### [`decrypt.h`](include/decrypt.h) — Multi-Identity Threshold Decryption

| API | Description |
|-----|-------------|
| `PartDec(Ĉ, k, t_k)` | Partial decryption share for identity k |
| `FinDec(ED₁, …, ED_N)` | Combine shares and decode plaintext |

#### [`expand.h`](include/expand.h) — Ciphertext Expansion

Constructs N×N block ciphertext matrix:
- Diagonal blocks: C (base ciphertext)
- Row i, col j (j ≠ i): Extend(U, b_i, b_j)

#### [`extend.h`](include/extend.h) — LWE Linear Combination

| API | Description |
|-----|-------------|
| `gsw_lcomb(U, coeffs, q)` | Homomorphic linear combination of LWE ciphertexts |
| `extend(U, b_i, b_j, q)` | X_j = GSW.LComb(U, b_j − b_i) |

#### [`unienc.h`](include/unienc.h) — Universal Encoding

Mask scheme: 𝒞 = A·R + μ_D·G, each bit of R encrypted via LWE.

---

## Configuration

All parameters are injected at **CMake configure time** via the template `include/params_cfg.h.in`, which generates `params_cfg.h` in the build directory. This is read by `include/unified_params.h`.

| CMake Variable | Default | Description |
|---------------|---------|-------------|
| `UNIFIED_PARAMS` | `ON` | Enable `params_cfg.h` generation |
| `MP12_N` | `8` | Lattice dimension n |
| `MP12_Q` | `257` | Modulus q |
| `MP12_B` | `2` | Gadget base b |
| `MP12_SIGMA` | `1.0` | Trapdoor Gaussian width σ |
| `UNIENC_SIGMA` | `3.2` | LWE noise width |
| `MID_LAMBDA` | `8` | Security parameter λ |
| `MID_N_ID` | `3` | Number of identities N |
| `MID_D` | `1` | Circuit depth d |

**128-bit parameters**: `cmake .. -DMP12_N=512 -DMP12_Q=134219777 -DMP12_B=27 -DMID_LAMBDA=128`

---

## Troubleshooting

| Symptom | Cause | Resolution |
|---------|-------|-------------|
| `ld: undefined reference to run_del_tests` | Stale build after header changes | `rm -rf build && rebuild` |
| `Element not invertible mod q` | H not in GL_n(Z_q) in DelTrapGen | H must have determinant coprime to q |
| `"Unknown CMake command 'configure_file'"` | CMake < 3.22 | Upgrade CMake to ≥ 3.22 |
| `error: 'std::gcd' is not a member of 'std'` | Compiler < C++17 | Use C++20 compiler (GCC ≥ 11) |
| `-D` parameters not updated after cmake | Stale CMake cache | `rm -rf build && mkdir build && cd build` |
| `error: no matching function for call to 'make_mat'` (or `'random_mat'`) | Including headers in wrong order; missing `using namespace matops` | Ensure [`matops.h`](include/matops.h) is included; add `using namespace matops;` in test code |
| `ld: ... undefined reference to 'run_test_*'` | Adding a new test file but not listed in `CMakeLists.txt` | The `file(GLOB DEMO_SRCS src/*.cpp)` in [`CMakeLists.txt`](CMakeLists.txt:49-51) auto-discovers new `.cpp` files — ensure your file is under `src/` and re-run cmake |
| `error: '__int128' is not supported` | Compiling on 32-bit target or with `-m32` | `__int128` is used for overflow-safe modulus switching in [`powersof_modswitch.h`](include/powersof_modswitch.h); build on x86-64 |
| `#include "params_cfg.h" not found` | CMake configure step not re-run after branch switch or fresh checkout | Re-run `cmake ..` in the build directory to regenerate the header |
| Segfault in `mat_mul` / `mat_add` with large dimensions | Stack overflow from large `vector<vector<long>>` | Consider bumping `ulimit -s` or switching to heap-allocated matrix storage for n ≥ 1024 |
| Tests pass at demo params but fail at 128-bit params | q not prime, or σ too small for large n | Use the recommended 128-bit q values (`MP12_Q=134219777`); the current `Params::make` uses heuristic σ based on n·√log q |
| Output differs between Debug and Release builds | Undefined behavior (uninitialized memory, signed overflow) | Run under `-fsanitize=address,undefined` in Debug mode; the codebase uses `mod_pos()` consistently but check custom math helpers |

---

## Engineering Notes

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| **Single unified executable** | All test/benchmark `.cpp` files compile into a static library `libdemos.a`, linked by one [`main.cpp`](src/main.cpp) entry point. Simplifies CI, avoids binary sprawl, and enables inlining across translation units at LTO. |
| **C++20 over C++17** | `__int128` in [`powersof_modswitch.h`](include/powersof_modswitch.h) prevents overflow during rounding-based modulus switch. C++20 concepts would allow static assertion of matrix dimension compatibility. |
| **Zero external dependencies** | The only dependency is the C++ standard library. No NTL, FLINT, or OpenSSL. This guarantees portability and removes ABI compatibility headaches. |
| **i–k–j loop ordering** | [`matops.h`](include/matops.h) uses i–k–j (middle-product-first) traversal to maximize L1 cache hit rate. For 256×256 matrices, this yields a ~4× speedup over naive i–j–k. |
| **CMake `configure_file()` for params** | Compile-time parameters are injected via [`params_cfg.h.in`](include/params_cfg.h.in) → `params_cfg.h`. Avoids runtime config parsing and enables compiler constant-propagation optimizations. |
| **Modulus switching via `round_scale`** | Uses `__int128` intermediate to compute `⌊p·x/q⌉` without overflow, critical for correct IBE key extraction when q is large. |

### Profiling

For detailed performance analysis, build with `-DCMAKE_BUILD_TYPE=Release` and use:

```bash
# CPU micro-architecture counters
perf stat -e cycles,instructions,cache-misses,branches ./build/LatticeCryBenchmarking

# Flame graph (hotspot identification)
perf record -g ./build/LatticeCryBenchmarking
perf script | stackcollapse-perf.pl | flamegraph.pl > flame.svg

# Cache miss analysis
perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads \
    ./build/LatticeCryBenchmarking
```

### Contributing

1. New primitives go in `include/<module>.h` (header-only API) or `src/<module>.cpp` (implementation).
2. Test/benchmark files follow the naming convention `src/test_<module>.cpp` / `src/bench_<module>.cpp`.
3. Expose the entry function in the header and call it from [`main.cpp`](src/main.cpp:48-57).
4. All matrix operations must use the `matops` namespace to keep the codebase consistent.
5. Re-run `cmake ..` after adding new source files — the `file(GLOB)` in [`CMakeLists.txt`](CMakeLists.txt:49-51) auto-discovers them.

---

## License & Citation

**Author**: Ziyi Dong, 2026.

This project is provided for academic and benchmarking purposes. If you use this code in your research, please cite both this repository and the original works:

### Cite This Code

```bibtex
@software{dong2026_latticecrybenchmarking,
  author  = {Ziyi Dong},
  title   = {{LatticeCryBendmarking}: A C++ Lattice-Based Cryptography
             Benchmarking and Testing Framework},
  year    = {2026},
  url     = {https://github.com/ZIYIDONG/LatticeCryBendmarking},
  note    = {Implements MP12 trapdoors, GSW homomorphic evaluation,
             multi-identity threshold FHE, and related primitives}
}
```

### Primary References

```bibtex
@inproceedings{micciancio2012trapdoors,
  author    = {Daniele Micciancio and Chris Peikert},
  title     = {Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller},
  booktitle = {Advances in Cryptology -- EUROCRYPT 2012},
  year      = {2012},
  pages     = {700--718},
  publisher = {Springer},
  doi       = {10.1007/978-3-642-29011-4_41}
}

@inproceedings{agrawal2010efficient,
  author    = {Shweta Agrawal and Dan Boneh and Xavier Boyen},
  title     = {Efficient Lattice {(H)IBE} in the Standard Model},
  booktitle = {Advances in Cryptology -- EUROCRYPT 2010},
  year      = {2010},
  pages     = {553--572},
  publisher = {Springer},
  doi       = {10.1007/978-3-642-13190-5_28}
}

@inproceedings{gentry2013homomorphic,
  author    = {Craig Gentry and Amit Sahai and Brent Waters},
  title     = {Homomorphic Encryption from Learning with Errors:
               Conceptually-Simpler, Asymptotically-Faster, Attribute-Based},
  booktitle = {Advances in Cryptology -- CRYPTO 2013},
  year      = {2013},
  pages     = {75--92},
  publisher = {Springer},
  doi       = {10.1007/978-3-642-40041-4_5}
}
```