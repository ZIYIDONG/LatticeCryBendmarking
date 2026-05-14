# LatticeCryBenchmarking — Top-level build dispatcher
#
# Usage:
#   make plain-demo      Plain-LWE  demo params  → build_plain-LWE/
#   make plain-128       Plain-LWE  128-bit      → build_plain-LWE/
#   make ring-demo       Ring-LWE   demo params  → build_ring-LWE/
#   make ring-128        Ring-LWE   128-bit      → build_ring-LWE/
#   make module-demo     Module-LWE demo params  → build_module-LWE/
#   make module-128      Module-LWE 128-bit      → build_module-LWE/
#
#   make all-demo        Build all three variants with demo params
#   make all-128         Build all three variants with 128-bit params
#   make clean           Remove all build directories
#
# Under the hood: cmake -B <dir>  (creates dir if it doesn't exist)

CMAKE  := cmake
BUILD  := cmake --build
JOBS   := $(shell nproc 2>/dev/null || echo 4)

# ── Variant-to-source mapping ─────────────────────────────
# Each LWE variant has its own src_<variant>-LWE/ and include_<variant>-LWE/
# directories.  By default the top-level CMakeLists.txt builds the
# "plain" variant.  Future ring/module variants will need their own
# CMakeLists.txt or a parametrised build.

define build_variant
	@echo ">>> Building $(1)-LWE ($(2) params) -> build_$(1)-LWE/"
	$(CMAKE) -B build_$(1)-LWE -DCMAKE_BUILD_TYPE=Release $(3)
	$(BUILD) build_$(1)-LWE -j $(JOBS)
	@echo ">>> Done: build_$(1)-LWE/LatticeCryBenchmarking"
endef

# ── Plain-LWE ────────────────────────────────────────────
.PHONY: plain-demo
plain-demo:
	$(call build_variant,plain,demo,)

.PHONY: plain-128
plain-128:
	$(call build_variant,plain,128-bit,-DUSE_128BIT_PARAMS=ON)

# ── Ring-LWE (placeholder) ───────────────────────────────
.PHONY: ring-demo
ring-demo:
	$(call build_variant,ring,demo,)

.PHONY: ring-128
ring-128:
	$(call build_variant,ring,128-bit,-DUSE_128BIT_PARAMS=ON)

# ── Module-LWE (placeholder) ─────────────────────────────
.PHONY: module-demo
module-demo:
	$(call build_variant,module,demo,)

.PHONY: module-128
module-128:
	$(call build_variant,module,128-bit,-DUSE_128BIT_PARAMS=ON)

# ── Aggregate targets ────────────────────────────────────
.PHONY: all-demo
all-demo: plain-demo ring-demo module-demo

.PHONY: all-128
all-128: plain-128 ring-128 module-128

# ── Clean ────────────────────────────────────────────────
.PHONY: clean
clean:
	rm -rf build_plain-LWE build_ring-LWE build_module-LWE
	@echo ">>> All build directories removed"

# ── Help ─────────────────────────────────────────────────
.PHONY: help
help:
	@echo "LatticeCryBenchmarking — Build Targets"
	@echo ""
	@echo "  make plain-demo       Plain-LWE  demo  (n=8,  q=257)"
	@echo "  make plain-128        Plain-LWE  128-bit (n=512, q=134219777)"
	@echo "  make ring-demo        Ring-LWE   demo  (future)"
	@echo "  make ring-128         Ring-LWE   128-bit (future)"
	@echo "  make module-demo      Module-LWE demo  (future)"
	@echo "  make module-128       Module-LWE 128-bit (future)"
	@echo "  make all-demo         Build all three with demo params"
	@echo "  make all-128          Build all three with 128-bit params"
	@echo "  make clean            Remove all build directories"
