<p align="center">
  <img src="https://img.shields.io/badge/C%2B%2B-20-blue?logo=c%2B%2B" alt="C++20">
  <img src="https://img.shields.io/badge/CMake-3.22%2B-brightgreen?logo=cmake" alt="CMake 3.22+">
  <img src="https://img.shields.io/badge/license-Academic-lightgrey" alt="License">
  <img src="https://img.shields.io/badge/build-passing-success" alt="Build">
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
- [Module Reference](#module-reference)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
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
| Parameters not updated after `-D` | Stale CMake cache | `rm -rf build && mkdir build && cd build` |

---

## License & Citation

**Author**: Ziyi Dong, 2026.

This project is provided for academic and benchmarking purposes. If you use this code in your research, please cite the original works:

- Micciancio, D. and Peikert, C. *"Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller"*, EUROCRYPT 2012.
- Agrawal, S., Boneh, D., and Boyen, X. *"Efficient Lattice (H)IBE in the Standard Model"*, EUROCRYPT 2010.
- Gentry, C., Sahai, A., and Waters, B. *"Homomorphic Encryption from Learning with Errors"*, CRYPTO 2013.