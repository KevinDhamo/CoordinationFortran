# Compiler and flags
FC = gfortran
FFLAGS_BASE = -cpp -fopenmp

# Release optimization
FFLAGS_OPT = $(FFLAGS_BASE) \
             -O3 \
             -march=native \
             -ffast-math \
             -funroll-loops \
             -ftree-vectorize \
             -flto

# Debug flags
FFLAGS_DEBUG = $(FFLAGS_BASE) \
               -g \
               -fcheck=all \
               -fbacktrace \
               -Wall \
               -Wextra \
               -Warray-temporaries \
               -fimplicit-none \
               -finit-real=nan \
               -finit-integer=-99999999

# Source files in dependency order
SRCS = types_mod.f90 \
       error_mod.f90 \
       config_mod.f90 \
       cell_list_mod.f90 \
       atom_selection_mod.f90 \
       coordination_mod.f90 \
       progress_mod.f90 \
       data_file_io_mod.f90 \
       trajectory_io_mod.f90 \
       output_cache_mod.f90 \
       output_io_mod.f90 \
       io_mod.f90 \
       main.f90

# Benchmark source
BENCHMARK_SRC = benchmark_main.f90

# Object files
OBJS = $(SRCS:.f90=.o)

# Module files
MODS = $(SRCS:.f90=.mod)

# Common objects (all except main.o)
COMMON_OBJS = $(filter-out main.o,$(OBJS))

# Executables
PROG = coordination_analysis
PROG_DEBUG = $(PROG)_debug
BENCHMARK = benchmark_coordination

# Default target
all: release debug benchmark

# Release build (optimized)
release: FFLAGS = $(FFLAGS_OPT)
release: $(PROG)

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(PROG_DEBUG)

# Benchmark build
benchmark: FFLAGS = $(FFLAGS_OPT)
benchmark: $(BENCHMARK)

# Build rules
$(PROG): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(PROG_DEBUG): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(BENCHMARK): $(COMMON_OBJS) benchmark_main.o
	$(FC) $(FFLAGS) -o $@ $^

# Dependencies
types_mod.o: types_mod.f90
error_mod.o: error_mod.f90
config_mod.o: config_mod.f90 types_mod.o error_mod.o
cell_list_mod.o: cell_list_mod.f90 types_mod.o config_mod.o error_mod.o
atom_selection_mod.o: atom_selection_mod.f90 types_mod.o config_mod.o error_mod.o
coordination_mod.o: coordination_mod.f90 types_mod.o config_mod.o cell_list_mod.o error_mod.o atom_selection_mod.o
progress_mod.o: progress_mod.f90 types_mod.o config_mod.o
data_file_io_mod.o: data_file_io_mod.f90 types_mod.o config_mod.o error_mod.o
trajectory_io_mod.o: trajectory_io_mod.f90 types_mod.o config_mod.o error_mod.o
output_cache_mod.o: output_cache_mod.f90 types_mod.o config_mod.o error_mod.o
output_io_mod.o: output_io_mod.f90 types_mod.o config_mod.o error_mod.o coordination_mod.o cell_list_mod.o atom_selection_mod.o
io_mod.o: io_mod.f90 types_mod.o trajectory_io_mod.o data_file_io_mod.o output_io_mod.o progress_mod.o
main.o: main.f90 types_mod.o config_mod.o cell_list_mod.o coordination_mod.o io_mod.o error_mod.o atom_selection_mod.o
benchmark_main.o: benchmark_main.f90 types_mod.o config_mod.o cell_list_mod.o coordination_mod.o io_mod.o error_mod.o atom_selection_mod.o

# Compilation rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean rules
clean:
	rm -f *.o *.mod $(PROG) $(PROG_DEBUG) $(BENCHMARK)

# Help target
help:
	@echo "Available targets:"
	@echo "  make        - Build all versions"
	@echo "  make release - Build optimized version"
	@echo "  make debug   - Build debug version"
	@echo "  make benchmark - Build benchmark version"
	@echo "  make clean   - Remove compiled files"
	@echo "  make help    - Show this help message"

.PHONY: all release debug benchmark clean help