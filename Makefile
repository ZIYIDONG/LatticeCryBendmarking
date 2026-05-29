# LatticeCryBenchmarking — Top-level build & test dispatcher
#
# Submodules:
#   plain-LWE   Matrix-based LWE  (MP12 trapdoors, GSW eval, FHE decrypt, benchmarks)
#   ibags       Ring-based RLWE/NTRU  (IBAGS protocols, NTT, trapgen, unit tests)
#
# Usage:
#   make build-L1           Build both submodules with NIST Level 1 params
#   make build-demo         Build both submodules with Demo params
#   make build-L3           Build both submodules with NIST Level 3
#   make build-L5           Build both submodules with NIST Level 5
#
#   make run-plain-L1       Run Plain-LWE  (L1, --no-decrypt for speed)
#   make run-plain-L1-full  Run Plain-LWE  (L1, full including decrypt)
#   make run-ibags-L1       Run IBAGS test (L1, quiet mode)
#   make run-ibags-L1-v     Run IBAGS test (L1, verbose mode)
#
#   make all-L1             Build + Run both  (L1, skip decrypt)
#   make all-demo           Build + Run both  (Demo)
#
#   make log-L1             Build + Run both → logs/L1-<ts>/  (auto skip decrypt)
#
#   make clean              Remove all build directories
#   make clean-logs         Remove all log directories
#   make help

CMAKE    := cmake
BUILD    := cmake --build
JOBS     := $(shell nproc 2>/dev/null || echo 4)
LOGS_DIR := logs
TIMESTAMP := $(shell date +%Y%m%d-%H%M%S)

# ── Internal: cmake configure + build (uses single build dir per level) ──
define build_all
	@echo "============================================================"
	@echo "  Building LatticeCryBenchmarking (level=$(1))"
	@echo "  Targets: LatticeCryBenchmarking + test_ibags_all"
	@echo "============================================================"
	$(CMAKE) -B build_$(1) -DCMAKE_BUILD_TYPE=Release $(2)
	@echo ""
	$(BUILD) build_$(1) -j $(JOBS)
	@echo ""
	@echo "  Binaries:"
	@echo "    build_$(1)/LatticeCryBenchmarking  (plain-LWE)"
	@echo "    build_$(1)/test_ibags_all           (IBAGS)"
	@echo "============================================================"
endef

# ── Build targets ──────────────────────────────────────────────────
.PHONY: build-demo build-L1 build-L3 build-L5
build-demo: ; $(call build_all,demo,)
build-L1:   ; $(call build_all,L1,-DSECURITY_LEVEL=1)
build-L3:   ; $(call build_all,L3,-DSECURITY_LEVEL=3)
build-L5:   ; $(call build_all,L5,-DSECURITY_LEVEL=5)

# ── Run: Plain-LWE (default: skip decrypt for speed; use *-full for all) ─
.PHONY: run-plain-demo run-plain-L1 run-plain-L3 run-plain-L5
run-plain-demo: build-demo ; @echo ">>> Running Plain-LWE (Demo)"  && build_demo/LatticeCryBenchmarking
run-plain-L1:   build-L1   ; @echo ">>> Running Plain-LWE (L1, --no-decrypt)"  && build_L1/LatticeCryBenchmarking --no-decrypt
run-plain-L3:   build-L3   ; @echo ">>> Running Plain-LWE (L3, --no-decrypt)"  && build_L3/LatticeCryBenchmarking --no-decrypt
run-plain-L5:   build-L5   ; @echo ">>> Running Plain-LWE (L5, --no-decrypt)"  && build_L5/LatticeCryBenchmarking --no-decrypt

.PHONY: run-plain-demo-full run-plain-L1-full run-plain-L3-full run-plain-L5-full
run-plain-demo-full: build-demo ; @echo ">>> Running Plain-LWE (Demo, full)"  && build_demo/LatticeCryBenchmarking
run-plain-L1-full:   build-L1   ; @echo ">>> Running Plain-LWE (L1, full)"    && build_L1/LatticeCryBenchmarking
run-plain-L3-full:   build-L3   ; @echo ">>> Running Plain-LWE (L3, full)"    && build_L3/LatticeCryBenchmarking
run-plain-L5-full:   build-L5   ; @echo ">>> Running Plain-LWE (L5, full)"    && build_L5/LatticeCryBenchmarking

# ── Run: IBAGS (quiet by default, add -v for verbose) ────────────
.PHONY: run-ibags-demo run-ibags-L1 run-ibags-L3 run-ibags-L5
run-ibags-demo: build-demo ; @echo ">>> Running IBAGS  (Demo)"    && build_demo/test_ibags_all -q
run-ibags-L1:   build-L1   ; @echo ">>> Running IBAGS  (L1)"      && build_L1/test_ibags_all -q
run-ibags-L3:   build-L3   ; @echo ">>> Running IBAGS  (L3)"      && build_L3/test_ibags_all -q
run-ibags-L5:   build-L5   ; @echo ">>> Running IBAGS  (L5)"      && build_L5/test_ibags_all -q

.PHONY: run-ibags-demo-v run-ibags-L1-v run-ibags-L3-v run-ibags-L5-v
run-ibags-demo-v: build-demo ; @echo ">>> Running IBAGS  (Demo, verbose)"  && build_demo/test_ibags_all
run-ibags-L1-v:   build-L1   ; @echo ">>> Running IBAGS  (L1, verbose)"    && build_L1/test_ibags_all
run-ibags-L3-v:   build-L3   ; @echo ">>> Running IBAGS  (L3, verbose)"    && build_L3/test_ibags_all
run-ibags-L5-v:   build-L5   ; @echo ">>> Running IBAGS  (L5, verbose)"    && build_L5/test_ibags_all

# ── Bench: IBAGS ──────────────────────────────────────────────────
.PHONY: bench-ibags-demo bench-ibags-L1 bench-ibags-L3 bench-ibags-L5
bench-ibags-demo: build-demo ; @echo ">>> Benchmarking IBAGS (Demo)" && build_demo/bench_ibags
bench-ibags-L1:   build-L1   ; @echo ">>> Benchmarking IBAGS (L1)"   && build_L1/bench_ibags
bench-ibags-L3:   build-L3   ; @echo ">>> Benchmarking IBAGS (L3)"   && build_L3/bench_ibags
bench-ibags-L5:   build-L5   ; @echo ">>> Benchmarking IBAGS (L5)"   && build_L5/bench_ibags

# ── Log targets (终端 + 日志文件同时输出, auto --no-decrypt) ────
.PHONY: log-demo log-L1 log-L3 log-L5
log-demo: build-demo
	@mkdir -p $(LOGS_DIR)/demo-$(TIMESTAMP)
	@echo ">>> Running Plain-LWE (Demo)  ->  $(LOGS_DIR)/demo-$(TIMESTAMP)/plain-LWE.log"
	@build_demo/LatticeCryBenchmarking 2>&1 | tee $(LOGS_DIR)/demo-$(TIMESTAMP)/plain-LWE.log
	@echo ">>> Running IBAGS  (Demo)     ->  $(LOGS_DIR)/demo-$(TIMESTAMP)/ibags.log"
	@build_demo/test_ibags_all -q 2>&1 | tee $(LOGS_DIR)/demo-$(TIMESTAMP)/ibags.log
	@echo ""
	@echo ">>> Logs saved in $(LOGS_DIR)/demo-$(TIMESTAMP)/"

log-L1: build-L1
	@mkdir -p $(LOGS_DIR)/L1-$(TIMESTAMP)
	@echo ">>> Running Plain-LWE (L1, --no-decrypt)  ->  $(LOGS_DIR)/L1-$(TIMESTAMP)/plain-LWE.log"
	@build_L1/LatticeCryBenchmarking --no-decrypt 2>&1 | tee $(LOGS_DIR)/L1-$(TIMESTAMP)/plain-LWE.log
	@echo ">>> Running IBAGS  (L1)     ->  $(LOGS_DIR)/L1-$(TIMESTAMP)/ibags.log"
	@build_L1/test_ibags_all -q 2>&1 | tee $(LOGS_DIR)/L1-$(TIMESTAMP)/ibags.log
	@echo ""
	@echo ">>> Logs saved in $(LOGS_DIR)/L1-$(TIMESTAMP)/"

log-L3: build-L3
	@mkdir -p $(LOGS_DIR)/L3-$(TIMESTAMP)
	@echo ">>> Running Plain-LWE (L3, --no-decrypt)  ->  $(LOGS_DIR)/L3-$(TIMESTAMP)/plain-LWE.log"
	@build_L3/LatticeCryBenchmarking --no-decrypt 2>&1 | tee $(LOGS_DIR)/L3-$(TIMESTAMP)/plain-LWE.log
	@echo ">>> Running IBAGS  (L3)     ->  $(LOGS_DIR)/L3-$(TIMESTAMP)/ibags.log"
	@build_L3/test_ibags_all -q 2>&1 | tee $(LOGS_DIR)/L3-$(TIMESTAMP)/ibags.log
	@echo ""
	@echo ">>> Logs saved in $(LOGS_DIR)/L3-$(TIMESTAMP)/"

log-L5: build-L5
	@mkdir -p $(LOGS_DIR)/L5-$(TIMESTAMP)
	@echo ">>> Running Plain-LWE (L5, --no-decrypt)  ->  $(LOGS_DIR)/L5-$(TIMESTAMP)/plain-LWE.log"
	@build_L5/LatticeCryBenchmarking --no-decrypt 2>&1 | tee $(LOGS_DIR)/L5-$(TIMESTAMP)/plain-LWE.log
	@echo ">>> Running IBAGS  (L5)     ->  $(LOGS_DIR)/L5-$(TIMESTAMP)/ibags.log"
	@build_L5/test_ibags_all -q 2>&1 | tee $(LOGS_DIR)/L5-$(TIMESTAMP)/ibags.log
	@echo ""
	@echo ">>> Logs saved in $(LOGS_DIR)/L5-$(TIMESTAMP)/"

# ── Aggregate targets ──────────────────────────────────────────────
.PHONY: all-demo all-L1 all-L3 all-L5
all-demo: run-plain-demo run-ibags-demo
all-L1:   run-plain-L1   run-ibags-L1
all-L3:   run-plain-L3   run-ibags-L3
all-L5:   run-plain-L5   run-ibags-L5

# ── Clean ──────────────────────────────────────────────────────────
.PHONY: clean
clean:
	rm -rf build_demo build_L1 build_L3 build_L5 build_plain-* build_module-* build_ring-*
	@echo ">>> All build directories removed"

.PHONY: clean-logs
clean-logs:
	rm -rf $(LOGS_DIR)
	@echo ">>> All log directories removed"

# ── Help ───────────────────────────────────────────────────────────
.PHONY: help
help:
	@echo "LatticeCryBenchmarking — Two Submodules"
	@echo ""
	@echo "┌─ Build ───────────────────────────────────────────────┐"
	@echo "│ make build-demo     Demo params  (n=64,   ~40-bit)    │"
	@echo "│ make build-L1       NIST L1      (n=512,  ~128-bit)   │"
	@echo "│ make build-L3       NIST L3      (n=768,  ~192-bit)   │"
	@echo "│ make build-L5       NIST L5      (n=1024, ~256-bit)   │"
	@echo "├─ Run: Plain-LWE ──────────────────────────────────────┤"
	@echo "│ make run-plain-L1       Quick  (--no-decrypt)         │"
	@echo "│ make run-plain-L1-full  Full   (includes slow decrypt) │"
	@echo "│ make run-plain-demo     Demo params                   │"
	@echo "├─ Run: IBAGS ──────────────────────────────────────────┤"
	@echo "│ make run-ibags-L1       Quiet mode (~50 lines)        │"
	@echo "│ make run-ibags-L1-v     Verbose mode                  │"
	@echo "├─ All-in-one ──────────────────────────────────────────┤"
	@echo "│ make all-L1             Build + Run both (L1, quick)  │"
	@echo "│ make all-demo           Build + Run both (Demo)       │"
	@echo "├─ Log file output ─────────────────────────────────────┤"
	@echo "│ make log-L1             Run both -> logs/L1-<ts>/     │"
	@echo "│ make log-demo           Run both -> logs/demo-<ts>/   │"
	@echo "├─ Utilities ───────────────────────────────────────────┤"
	@echo "│ make clean              Remove all build directories  │"
	@echo "│ make clean-logs         Remove all log directories    │"
	@echo "│ make help               Show this help                │"
	@echo "└───────────────────────────────────────────────────────┘"
	@echo ""
	@echo "Binaries per level (single build directory):"
	@echo "  build_L1/LatticeCryBenchmarking   build_L1/test_ibags_all"
	@echo ""
	@echo "Log files per level (timestamped):"
	@echo "  logs/L1-YYYYMMDD-HHMMSS/plain-LWE.log"
	@echo "  logs/L1-YYYYMMDD-HHMMSS/ibags.log"
	@echo ""
	@echo "Decrypt skip: L1+ params auto-use --no-decrypt to avoid"
	@echo "  hours-long O(n^3) matrix ops. Use *-full targets for full run."
