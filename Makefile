# LatticeCryBenchmarking — Top-level build dispatcher
#
# LWE variants in structural-richness order:
#   plain  (no structure, matrix ops, slowest)  →  module  (block circulant, balanced)  →  ring  (full algebraic, fastest)
#
# Usage:
#   make plain-demo        Plain-LWE  demo  (n=8)    → build_plain-demo/
#   make plain-Level-1     Plain-LWE  L1   (n=512)   → build_plain-Level-1/
#   make plain-Level-3     Plain-LWE  L3   (n=768)   → build_plain-Level-3/
#   make plain-Level-5     Plain-LWE  L5   (n=1024)  → build_plain-Level-5/
#   make module-demo       Module-LWE demo           → build_module-demo/
#   make module-Level-1    Module-LWE L1             → build_module-Level-1/
#   make module-Level-3    Module-LWE L3             → build_module-Level-3/
#   make module-Level-5    Module-LWE L5             → build_module-Level-5/
#   make ring-demo         Ring-LWE   demo           → build_ring-demo/
#   make ring-Level-1      Ring-LWE   L1             → build_ring-Level-1/
#   make ring-Level-3      Ring-LWE   L3             → build_ring-Level-3/
#   make ring-Level-5      Ring-LWE   L5             → build_ring-Level-5/
#
#   make all-demo          Build all three with demo params
#   make all-Level-1       Build all three with NIST Level 1
#   make all-Level-3       Build all three with NIST Level 3
#   make all-Level-5       Build all three with NIST Level 5
#   make clean             Remove all build directories

CMAKE  := cmake
BUILD  := cmake --build
JOBS   := $(shell nproc 2>/dev/null || echo 4)

define build_variant
	@echo ">>> Building $(1)-LWE ($(2)) -> build_$(1)-$(2)/"
	$(CMAKE) -B build_$(1)-$(2) -DCMAKE_BUILD_TYPE=Release $(3)
	$(BUILD) build_$(1)-$(2) -j $(JOBS)
	@echo ">>> Done: build_$(1)-$(2)/LatticeCryBenchmarking"
endef

# ── Plain-LWE ────────────────────────────────────────────
.PHONY: plain-demo plain-Level-1 plain-Level-3 plain-Level-5
plain-demo:   ; $(call build_variant,plain,demo,)
plain-Level-1:; $(call build_variant,plain,Level-1,-DSECURITY_LEVEL=1)
plain-Level-3:; $(call build_variant,plain,Level-3,-DSECURITY_LEVEL=3)
plain-Level-5:; $(call build_variant,plain,Level-5,-DSECURITY_LEVEL=5)

# ── Module-LWE (placeholder) ─────────────────────────────
.PHONY: module-demo module-Level-1 module-Level-3 module-Level-5
module-demo:    ; $(call build_variant,module,demo,)
module-Level-1: ; $(call build_variant,module,Level-1,-DSECURITY_LEVEL=1)
module-Level-3: ; $(call build_variant,module,Level-3,-DSECURITY_LEVEL=3)
module-Level-5: ; $(call build_variant,module,Level-5,-DSECURITY_LEVEL=5)

# ── Ring-LWE (placeholder) ───────────────────────────────
.PHONY: ring-demo ring-Level-1 ring-Level-3 ring-Level-5
ring-demo:    ; $(call build_variant,ring,demo,)
ring-Level-1: ; $(call build_variant,ring,Level-1,-DSECURITY_LEVEL=1)
ring-Level-3: ; $(call build_variant,ring,Level-3,-DSECURITY_LEVEL=3)
ring-Level-5: ; $(call build_variant,ring,Level-5,-DSECURITY_LEVEL=5)

# ── Aggregate targets ────────────────────────────────────
.PHONY: all-demo all-Level-1 all-Level-3 all-Level-5
all-demo:    plain-demo    module-demo    ring-demo
all-Level-1: plain-Level-1 module-Level-1 ring-Level-1
all-Level-3: plain-Level-3 module-Level-3 ring-Level-3
all-Level-5: plain-Level-5 module-Level-5 ring-Level-5

# ── Clean ────────────────────────────────────────────────
.PHONY: clean
clean:
	rm -rf build_plain-* build_module-* build_ring-*
	@echo ">>> All build directories removed"

# ── Help ─────────────────────────────────────────────────
.PHONY: help
help:
	@echo "LatticeCryBenchmarking — Build Targets (plain → module → ring)"
	@echo ""
	@echo "  make plain-demo         Plain-LWE  Demo      (n=8,   q=257)"
	@echo "  make plain-Level-1      Plain-LWE  NIST L1   (n=512, q~2^27)"
	@echo "  make plain-Level-3      Plain-LWE  NIST L3   (n=768, q~2^32)"
	@echo "  make plain-Level-5      Plain-LWE  NIST L5   (n=1024,q~2^32)"
	@echo "  make module-*           Module-LWE (future)"
	@echo "  make ring-*             Ring-LWE   (future)"
	@echo "  make all-demo/Level-*   Build all three variants"
	@echo "  make clean              Remove all build directories"
