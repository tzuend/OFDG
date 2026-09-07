MFEM_DIR ?= ../mfem
MFEM_BUILD_DIR ?= $(MFEM_DIR)
MFEM_INSTALL_DIR ?= $(MFEM_DIR)

CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
                    $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))
-include $(CONFIG_MK)

BUILD_DIR := build
LOCAL_QUICK_DIR := measurements/study/local-quick
DEBUG_DIR := $(BUILD_DIR)/debug
RELEASE_DIR := $(BUILD_DIR)/release
DEBUG_TEST_DIR := $(DEBUG_DIR)/tests
RELEASE_TEST_DIR := $(RELEASE_DIR)/tests
DEBUG_EXAMPLE_DIR := $(DEBUG_DIR)/examples
RELEASE_EXAMPLE_DIR := $(RELEASE_DIR)/examples

CPPFLAGS := $(MFEM_CPPFLAGS) $(MFEM_INCFLAGS) \
	-DOFDG_MFEM_DATA_DIR=\"$(abspath $(MFEM_DIR))/data\"
CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra -Wpedantic \
	-Wno-unused-parameter -Wno-unused-variable
DEBUG_FLAGS := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined
RELEASE_FLAGS := -O3 -DNDEBUG
SANITIZER_LIBS := -fsanitize=address,undefined

TEST_NAMES := test_curved_geometry test_ofdg_geometry test_mixed_geometry test_kxrcf \
	test_ofdg_regression test_oedg_2024 \
	test_euler_positivity test_face_physics test_rk4_cadence test_experiment_rk
DEBUG_TESTS := $(addprefix $(DEBUG_TEST_DIR)/,$(TEST_NAMES))
RELEASE_TESTS := $(addprefix $(RELEASE_TEST_DIR)/,$(TEST_NAMES))
EXAMPLE_NAMES := advection burgers euler
DEBUG_EXAMPLES := $(addprefix $(DEBUG_EXAMPLE_DIR)/,$(EXAMPLE_NAMES))
RELEASE_EXAMPLES := $(addprefix $(RELEASE_EXAMPLE_DIR)/,$(EXAMPLE_NAMES))
FILTER_HEADERS := src/ofdg.hpp src/kxrcf.hpp src/oedg_2024.hpp \
	src/face_physics.hpp src/study_filter.hpp src/study_method.hpp
FILTER_CORE_HEADERS := src/curved_geometry.hpp src/ofdg_core.hpp src/ofdg.hpp src/oedg_2024.hpp \
	src/kxrcf.hpp src/face_physics.hpp
DEBUG_FILTER_OBJECTS := $(DEBUG_DIR)/src/ofdg_core.o $(DEBUG_DIR)/src/ofdg.o \
	$(DEBUG_DIR)/src/oedg_2024.o $(DEBUG_DIR)/src/kxrcf.o
RELEASE_FILTER_OBJECTS := $(RELEASE_DIR)/src/ofdg_core.o $(RELEASE_DIR)/src/ofdg.o \
	$(RELEASE_DIR)/src/oedg_2024.o $(RELEASE_DIR)/src/kxrcf.o
DEBUG_FILTER_LIBRARY := $(DEBUG_DIR)/libofdg.a
RELEASE_FILTER_LIBRARY := $(RELEASE_DIR)/libofdg.a

.PHONY: all test test-parallel test-public-headers test-release test-study examples release \
	benchmark-kxrcf reference-solver study-plan study-quick study-full \
	figures figures-draft report report-draft report-final clean

all: test

test: $(DEBUG_TESTS)
	@for test in $(DEBUG_TESTS); do $$test || exit 1; done
	$(MAKE) test-parallel

test-parallel: $(DEBUG_TEST_DIR)/test_parallel_consistency $(RELEASE_EXAMPLES)
	mpirun -np 2 $<
	PYTHONDONTWRITEBYTECODE=1 python3 tests/test_parallel_examples.py \
		--example-dir $(RELEASE_EXAMPLE_DIR)

test-release: test-public-headers $(RELEASE_TESTS)
	@for test in $(RELEASE_TESTS); do $$test || exit 1; done

test-public-headers:
	@for header in face_physics.hpp ofdg.hpp kxrcf.hpp oedg_2024.hpp \
		study_method.hpp study_filter.hpp; do \
		$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) -fsyntax-only -x c++ \
			-include src/$$header /dev/null || exit 1; \
	done
	@! grep -R -n '^[[:space:]]*using namespace mfem' src
	@! grep -nE 'MFEM_(VERIFY|ASSERT)|(^|[^[:alnum:]_])assert[[:space:]]*\(' \
		src/ofdg.hpp src/ofdg.cpp src/ofdg_core.hpp src/ofdg_core.cpp \
		src/kxrcf.hpp src/kxrcf.cpp
	@! grep -nE 'chrono|InternalTiming|INTERNAL_TIMING|ResetInternalTimings|PrintInternalTimings' \
		src/ofdg.hpp src/ofdg.cpp src/ofdg_core.hpp src/ofdg_core.cpp \
		src/kxrcf.hpp src/kxrcf.cpp
	@! grep -nE '^[[:space:]]*using mfem::' src/*.hpp
	@! grep -nE '#include "(kxrcf|oedg_2024|study_filter|study_method)\.hpp"' \
		src/ofdg.hpp
	@test ! -e src/ofdg_serial_optimized.hpp

test-study: reference-solver
	PYTHONDONTWRITEBYTECODE=1 python3 tests/test_study_analysis.py
	PYTHONDONTWRITEBYTECODE=1 python3 tests/test_study_design.py
	PYTHONDONTWRITEBYTECODE=1 python3 tests/test_reference_solver.py $(RELEASE_DIR)/benchmarks/euler_reference_weno

examples: $(DEBUG_EXAMPLES)

release: $(RELEASE_EXAMPLES)

benchmark-kxrcf: $(RELEASE_DIR)/benchmarks/kxrcf
	$<

reference-solver: $(RELEASE_DIR)/benchmarks/euler_reference_weno

study-quick: release reference-solver
	python3 scripts/run_study.py --profile quick --skip-build --output $(LOCAL_QUICK_DIR)
	python3 scripts/analyze_study.py --profile quick \
		--input $(LOCAL_QUICK_DIR)/results.csv --report-mode draft

study-plan:
	python3 scripts/run_study.py --profile full --dry-run

# The runner checks approval before it builds or launches anything.
study-full:
	python3 scripts/run_study.py --profile full
	python3 scripts/analyze_study.py --profile full --report-mode final

figures: figures-draft

figures-draft:
	python3 scripts/analyze_study.py --profile quick --report-mode draft

$(RELEASE_DIR)/benchmarks/euler_reference_weno: benchmarks/euler_reference_weno.cpp | $(RELEASE_DIR)/benchmarks
	$(MFEM_CXX) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@

$(DEBUG_TEST_DIR) $(RELEASE_TEST_DIR) $(DEBUG_EXAMPLE_DIR) $(RELEASE_EXAMPLE_DIR):
	mkdir -p $@

$(RELEASE_DIR)/benchmarks:
	mkdir -p $@

$(DEBUG_DIR)/src $(RELEASE_DIR)/src:
	mkdir -p $@

$(DEBUG_DIR)/src/ofdg.o: src/ofdg.cpp $(FILTER_CORE_HEADERS) | $(DEBUG_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) -c $< -o $@

$(DEBUG_DIR)/src/ofdg_core.o: src/ofdg_core.cpp $(FILTER_CORE_HEADERS) | $(DEBUG_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) -c $< -o $@

$(DEBUG_DIR)/src/oedg_2024.o: src/oedg_2024.cpp $(FILTER_CORE_HEADERS) | $(DEBUG_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) -c $< -o $@

$(DEBUG_DIR)/src/kxrcf.o: src/kxrcf.cpp $(FILTER_CORE_HEADERS) | $(DEBUG_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) -c $< -o $@

$(RELEASE_DIR)/src/ofdg.o: src/ofdg.cpp $(FILTER_CORE_HEADERS) | $(RELEASE_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) -c $< -o $@

$(RELEASE_DIR)/src/ofdg_core.o: src/ofdg_core.cpp $(FILTER_CORE_HEADERS) | $(RELEASE_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) -c $< -o $@

$(RELEASE_DIR)/src/oedg_2024.o: src/oedg_2024.cpp $(FILTER_CORE_HEADERS) | $(RELEASE_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) -c $< -o $@

$(RELEASE_DIR)/src/kxrcf.o: src/kxrcf.cpp $(FILTER_CORE_HEADERS) | $(RELEASE_DIR)/src
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) -c $< -o $@

$(DEBUG_FILTER_LIBRARY): $(DEBUG_FILTER_OBJECTS)
	ar rcs $@ $^

$(RELEASE_FILTER_LIBRARY): $(RELEASE_FILTER_OBJECTS)
	ar rcs $@ $^

$(RELEASE_DIR)/benchmarks/kxrcf: benchmarks/kxrcf_benchmark.cpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_DIR)/benchmarks
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(DEBUG_TEST_DIR)/test_ofdg_geometry: tests/test_ofdg_geometry.cpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_kxrcf: tests/test_kxrcf.cpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_mixed_geometry: tests/test_mixed_geometry.cpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_ofdg_regression: tests/test_ofdg_regression.cpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_oedg_2024: tests/test_oedg_2024.cpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_euler_positivity: tests/test_euler_positivity.cpp src/euler_positivity.hpp src/curved_geometry.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_face_physics: tests/test_face_physics.cpp src/face_physics.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_rk4_cadence: tests/test_rk4_cadence.cpp src/rk4.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_experiment_rk: tests/test_experiment_rk.cpp src/experiment_rk.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_parallel_consistency: tests/test_parallel_consistency.cpp $(FILTER_HEADERS) src/euler_positivity.hpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(RELEASE_TEST_DIR)/test_ofdg_geometry: tests/test_ofdg_geometry.cpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_kxrcf: tests/test_kxrcf.cpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_mixed_geometry: tests/test_mixed_geometry.cpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_ofdg_regression: tests/test_ofdg_regression.cpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_oedg_2024: tests/test_oedg_2024.cpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_euler_positivity: tests/test_euler_positivity.cpp src/euler_positivity.hpp src/curved_geometry.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_face_physics: tests/test_face_physics.cpp src/face_physics.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_rk4_cadence: tests/test_rk4_cadence.cpp src/rk4.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_experiment_rk: tests/test_experiment_rk.cpp src/experiment_rk.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(DEBUG_EXAMPLE_DIR)/advection: examples/advection/example_advection.cpp src/conservation.hpp src/profile_output.hpp src/glvis_output.hpp $(FILTER_HEADERS) $(DEBUG_FILTER_LIBRARY) | $(DEBUG_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_EXAMPLE_DIR)/burgers: examples/burgers/burgers.cpp examples/euler/euler.hpp src/conservation.hpp src/profile_output.hpp src/glvis_output.hpp $(FILTER_HEADERS) $(DEBUG_FILTER_LIBRARY) | $(DEBUG_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_EXAMPLE_DIR)/euler: examples/euler/euler.cpp examples/euler/euler.hpp src/conservation.hpp src/profile_output.hpp src/glvis_output.hpp src/euler_positivity.hpp $(FILTER_HEADERS) $(DEBUG_FILTER_LIBRARY) | $(DEBUG_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(RELEASE_EXAMPLE_DIR)/advection: examples/advection/example_advection.cpp src/conservation.hpp src/profile_output.hpp src/glvis_output.hpp $(FILTER_HEADERS) $(RELEASE_FILTER_LIBRARY) | $(RELEASE_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_EXAMPLE_DIR)/burgers: examples/burgers/burgers.cpp examples/euler/euler.hpp src/conservation.hpp src/profile_output.hpp src/glvis_output.hpp $(FILTER_HEADERS) $(RELEASE_FILTER_LIBRARY) | $(RELEASE_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(RELEASE_EXAMPLE_DIR)/euler: examples/euler/euler.cpp examples/euler/euler.hpp src/conservation.hpp src/profile_output.hpp src/glvis_output.hpp src/euler_positivity.hpp $(FILTER_HEADERS) $(RELEASE_FILTER_LIBRARY) | $(RELEASE_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

report: report-draft

report-draft: report-preview-assets
	$(MAKE) -C report draft
	mkdir -p output/pdf
	cp report/thesis.pdf output/pdf/thesis-revised.pdf

.PHONY: report-preview-assets
report-preview-assets:
	python3 scripts/analyze_preview.py --output report/generated/preview --report
	python3 scripts/export_report_preview.py

report-final:
	python3 scripts/check_final_report.py
	$(MAKE) -C report final

clean:
	rm -rf $(BUILD_DIR)

$(RELEASE_TEST_DIR)/test_curved_geometry: tests/test_curved_geometry.cpp $(FILTER_CORE_HEADERS) src/euler_positivity.hpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

$(DEBUG_TEST_DIR)/test_curved_geometry: tests/test_curved_geometry.cpp $(FILTER_CORE_HEADERS) src/euler_positivity.hpp $(DEBUG_FILTER_LIBRARY) | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(DEBUG_FILTER_LIBRARY) $(MFEM_LIBS) $(SANITIZER_LIBS)

$(RELEASE_TEST_DIR)/probe_native_nurbs: tests/probe_native_nurbs.cpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_curved_euler: tests/test_curved_euler.cpp examples/euler/euler.hpp $(FILTER_CORE_HEADERS) src/euler_positivity.hpp $(RELEASE_FILTER_LIBRARY) | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)

.PHONY: test-curved
test-curved: $(RELEASE_TEST_DIR)/test_curved_geometry $(RELEASE_TEST_DIR)/test_curved_euler $(RELEASE_EXAMPLE_DIR)/advection
	python3 tests/run_bounded.py $(RELEASE_TEST_DIR)/test_curved_geometry
	python3 tests/run_bounded.py mpirun -np 2 $(RELEASE_TEST_DIR)/test_curved_geometry
	python3 tests/run_bounded.py $(RELEASE_TEST_DIR)/test_curved_euler
	python3 tests/run_bounded.py mpirun -np 2 $(RELEASE_TEST_DIR)/test_curved_euler
	python3 tests/run_bounded.py python3 tests/test_curved_advection.py $(RELEASE_EXAMPLE_DIR)/advection

$(RELEASE_TEST_DIR)/test_preview: tests/test_preview.cpp src/profile_output.hpp examples/euler/euler.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

.PHONY: study-preview
study-preview: release reference-solver
	python3 scripts/run_preview.py
	python3 scripts/analyze_preview.py

$(RELEASE_DIR)/benchmarks/filter_performance: benchmarks/filter_performance.cpp examples/euler/euler.hpp $(FILTER_HEADERS) $(RELEASE_FILTER_LIBRARY) | $(RELEASE_DIR)/benchmarks
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)
