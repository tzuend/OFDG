MFEM_DIR ?= ../mfem
MFEM_BUILD_DIR ?= $(MFEM_DIR)
MFEM_INSTALL_DIR ?= $(MFEM_DIR)

CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
                    $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))
-include $(CONFIG_MK)

BUILD_DIR := build
DEBUG_DIR := $(BUILD_DIR)/debug
RELEASE_DIR := $(BUILD_DIR)/release
DEBUG_TEST_DIR := $(DEBUG_DIR)/tests
RELEASE_TEST_DIR := $(RELEASE_DIR)/tests
DEBUG_EXAMPLE_DIR := $(DEBUG_DIR)/examples
RELEASE_EXAMPLE_DIR := $(RELEASE_DIR)/examples

CPPFLAGS := $(MFEM_CPPFLAGS) $(MFEM_INCFLAGS)
CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra -Wpedantic \
	-Wno-unused-parameter -Wno-unused-variable
DEBUG_FLAGS := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined
RELEASE_FLAGS := -O3 -DNDEBUG
SANITIZER_LIBS := -fsanitize=address,undefined

TEST_NAMES := test_ofdg_geometry test_kxrcf test_ofdg_regression test_oedg_2024 \
	test_euler_positivity test_face_physics test_rk4_cadence test_experiment_rk
DEBUG_TESTS := $(addprefix $(DEBUG_TEST_DIR)/,$(TEST_NAMES))
RELEASE_TESTS := $(addprefix $(RELEASE_TEST_DIR)/,$(TEST_NAMES))
EXAMPLE_NAMES := advection burgers euler
DEBUG_EXAMPLES := $(addprefix $(DEBUG_EXAMPLE_DIR)/,$(EXAMPLE_NAMES))
RELEASE_EXAMPLES := $(addprefix $(RELEASE_EXAMPLE_DIR)/,$(EXAMPLE_NAMES))

.PHONY: all test test-parallel test-release test-study examples release \
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

test-release: $(RELEASE_TESTS)
	@for test in $(RELEASE_TESTS); do $$test || exit 1; done

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
	python3 scripts/run_study.py --profile quick --skip-build
	python3 scripts/analyze_study.py --profile quick --report-mode draft

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

$(RELEASE_DIR)/benchmarks/kxrcf: benchmarks/kxrcf_benchmark.cpp src/kxrcf.hpp src/face_physics.hpp | $(RELEASE_DIR)/benchmarks
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) -DKXRCF_INTERNAL_TIMING $< -o $@ $(MFEM_LIBS)

$(DEBUG_TEST_DIR)/test_ofdg_geometry: tests/test_ofdg_geometry.cpp src/ofdg_serial_optimized.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_kxrcf: tests/test_kxrcf.cpp src/kxrcf.hpp src/face_physics.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_ofdg_regression: tests/test_ofdg_regression.cpp src/ofdg_serial_optimized.hpp src/face_physics.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_oedg_2024: tests/test_oedg_2024.cpp src/oedg_2024.hpp src/ofdg_serial_optimized.hpp src/face_physics.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_euler_positivity: tests/test_euler_positivity.cpp src/euler_positivity.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_face_physics: tests/test_face_physics.cpp src/face_physics.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_rk4_cadence: tests/test_rk4_cadence.cpp src/rk4.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_experiment_rk: tests/test_experiment_rk.cpp src/experiment_rk.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_TEST_DIR)/test_parallel_consistency: tests/test_parallel_consistency.cpp src/ofdg_serial_optimized.hpp src/oedg_2024.hpp src/euler_positivity.hpp src/kxrcf.hpp src/face_physics.hpp | $(DEBUG_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(RELEASE_TEST_DIR)/test_ofdg_geometry: tests/test_ofdg_geometry.cpp src/ofdg_serial_optimized.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_kxrcf: tests/test_kxrcf.cpp src/kxrcf.hpp src/face_physics.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_ofdg_regression: tests/test_ofdg_regression.cpp src/ofdg_serial_optimized.hpp src/face_physics.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_oedg_2024: tests/test_oedg_2024.cpp src/oedg_2024.hpp src/ofdg_serial_optimized.hpp src/face_physics.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_euler_positivity: tests/test_euler_positivity.cpp src/euler_positivity.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_face_physics: tests/test_face_physics.cpp src/face_physics.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_rk4_cadence: tests/test_rk4_cadence.cpp src/rk4.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_TEST_DIR)/test_experiment_rk: tests/test_experiment_rk.cpp src/experiment_rk.hpp | $(RELEASE_TEST_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(DEBUG_EXAMPLE_DIR)/advection: examples/advection/example_advection.cpp src/conservation.hpp src/glvis_output.hpp src/study_filter.hpp src/oedg_2024.hpp src/ofdg_serial_optimized.hpp src/kxrcf.hpp | $(DEBUG_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_EXAMPLE_DIR)/burgers: examples/burgers/burgers.cpp examples/euler/euler.hpp src/conservation.hpp src/glvis_output.hpp src/study_filter.hpp src/oedg_2024.hpp src/ofdg_serial_optimized.hpp src/kxrcf.hpp | $(DEBUG_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(DEBUG_EXAMPLE_DIR)/euler: examples/euler/euler.cpp examples/euler/euler.hpp src/conservation.hpp src/glvis_output.hpp src/study_filter.hpp src/oedg_2024.hpp src/euler_positivity.hpp src/ofdg_serial_optimized.hpp src/kxrcf.hpp | $(DEBUG_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(DEBUG_FLAGS) $< -o $@ $(MFEM_LIBS) $(SANITIZER_LIBS)

$(RELEASE_EXAMPLE_DIR)/advection: examples/advection/example_advection.cpp src/conservation.hpp src/glvis_output.hpp src/study_filter.hpp src/oedg_2024.hpp src/ofdg_serial_optimized.hpp src/kxrcf.hpp | $(RELEASE_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_EXAMPLE_DIR)/burgers: examples/burgers/burgers.cpp examples/euler/euler.hpp src/conservation.hpp src/glvis_output.hpp src/study_filter.hpp src/oedg_2024.hpp src/ofdg_serial_optimized.hpp src/kxrcf.hpp | $(RELEASE_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

$(RELEASE_EXAMPLE_DIR)/euler: examples/euler/euler.cpp examples/euler/euler.hpp src/conservation.hpp src/glvis_output.hpp src/study_filter.hpp src/oedg_2024.hpp src/euler_positivity.hpp src/ofdg_serial_optimized.hpp src/kxrcf.hpp | $(RELEASE_EXAMPLE_DIR)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(MFEM_LIBS)

report: report-draft

report-draft: figures-draft
	$(MAKE) -C report draft

report-final:
	python3 scripts/check_final_report.py
	$(MAKE) -C report final

clean:
	rm -rf $(BUILD_DIR)
