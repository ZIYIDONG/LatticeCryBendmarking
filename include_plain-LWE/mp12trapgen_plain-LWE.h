#ifndef MP12TRAPGEN_H
#define MP12TRAPGEN_H

#include "mp12_plain-LWE.h"
#include "mp12deltrapgen_plain-LWE.h"
#include "unified_params_plain-LWE.h"

#include <vector>

// Diagnostics / demos extracted from main.cpp (global scope)

void test_gadget(const mp12::Params& p);
void test_gadget_basis(const mp12::Params& p);
void test_sample_g(const mp12::Params& p);
void test_gen_trap(const mp12::Params& p);
void test_sample_pre(const mp12::Params& p);
void test_uniformity(const mp12::Params& p);
void test_full_roundtrip_large(const mp12::Params& p);

void bench_gen_trap(const mp12::Params& p);

void run_mp12_trap_tests(const mp12::Params& p);

double vec_norm(const mp12::Vec& v);

/* ── External test entry points from other translation units ── */
void run_test_encrypt();
void run_test_powersof();
void run_test_powersof_modswitch();
void run_test_frd();
void run_bench_matops();
void run_test_expand();
void run_test_eval();
void run_test_decrypt();
void run_test_security();
void run_bench_decrypt();

#endif // MP12TRAPGEN_H
