#ifndef MP12TRAP_H
#define MP12TRAP_H

#include "mp12.h"
#include "mp12deltrapgen.h"
#include "unified_params.h"

#include <vector>

namespace mp12 {

// Diagnostics / demos extracted from main.cpp

void test_gadget(const Params& p);
void test_gadget_basis(const Params& p);
void test_sample_g(const Params& p);
void test_gen_trap(const Params& p);
void test_sample_pre(const Params& p);
void test_uniformity(const Params& p);
void test_full_roundtrip_large(const Params& p);

void run_mp12_trap_tests(const Params& p);

double vec_norm(const Vec& v);

} // namespace mp12

/* ── External test entry points from other translation units ── */
void run_test_powersof();
void run_test_powersof_modswitch();
void run_test_frd();
void run_debug_frd();
void run_bench_matops();
void run_test_expand();
void run_test_eval();
void run_test_decrypt();
void run_bench_decrypt();

#endif // MP12TRAP_H
