.PHONY: clean distclean
.PHONY: test_integration
ARTIFACT_PATTERNS := \
	"*.so" \
	"*.pcm" \
	"*.d" \
	"*_C.dll" \
	"*_C_ACLiC_dict*" \
	"*_C*so" \
	"*_C*pcm" \
	"*_dict*" \
	"test.root"
SEARCH_DIRS := $(shell find . -type d -not -path '*/.git*')

define do_clean
	@echo "  [DIR] $(1)" ; \
	for pattern in $(ARTIFACT_PATTERNS); do \
	  find $(1) -maxdepth 1 -type f -name $$pattern -print -delete; \
	done;
endef

clean:
	@echo "[CLEAN] Removing temporary build artifacts..."
	$(foreach dir,$(SEARCH_DIRS),$(call do_clean,$(dir)))
distclean: clean

# Integration tests (ROOT macros in tests/)
TEST_SOURCES := $(wildcard tests/test_*.C)

test_integration:
	@echo "[TEST] Running integration tests..."
	@set -e; \
	for t in $(TEST_SOURCES); do \
	  bn=$$(basename $$t); \
	  printf "%-30s" "$$bn ..."; \
	  if root -l -b -q "$$t+()" >/dev/null 2>&1; then \
	    printf "\033[0;32mPASSED\033[0m\n"; \
	  else \
	    printf "\033[0;31mFAILED\033[0m\n"; \
	    exit 1; \
	  fi; \
	done
	@echo "[TEST] All integration tests passed."

