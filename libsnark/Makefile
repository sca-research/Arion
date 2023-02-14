#### BEGIN USER OPTIONS ####
# Enables debug mode
DEBUG := 0
# Choose which curve to use: CURVE_{ALT_BN128, BLS12_381, BN128, EDWARDS, MNT4, MNT6}
ELLIPTIC_CURVE := CURVE_ALT_BN128
# Enable micro-benchmarking during tests, useful to check performance variation across builds: 0/1
MEASURE_PERFORMANCE := 0
# Enable OpenMP parallelization: O/1
MULTICORE := 0
# Override default compiler and flags: 0/1
OVERRIDE_DEFAULT_CXX := 1
# Enable assembly code (improves performance, might be buggy): 0/1
USE_ASM := 0
#### END USER OPTIONS ####


#### BEGIN TARGETS ####
# Targets which only have tests
TARGETS_ONLYTEST :=
TARGETS_ONLYTEST += arion
TARGETS_ONLYTEST +=	arion_gadget
TARGETS_ONLYTEST += arion_v2
TARGETS_ONLYTEST +=	arion_v2_gadget
#TARGETS_ONLYTEST += abr_gadget
TARGETS_ONLYTEST += fixed_abr
TARGETS_ONLYTEST += fixed_mtree
#TARGETS_ONLYTEST += fixed_mtree_gadget
TARGETS_ONLYTEST += griffin
TARGETS_ONLYTEST += griffin_gadget
TARGETS_ONLYTEST += mimc256
TARGETS_ONLYTEST += mimc256_gadget
TARGETS_ONLYTEST += mimc512f
TARGETS_ONLYTEST += mimc512f_gadget
TARGETS_ONLYTEST += mimc512f2k
TARGETS_ONLYTEST += mimc512f2k_gadget
TARGETS_ONLYTEST += mtree
TARGETS_ONLYTEST += mtree_gadget
TARGETS_ONLYTEST += poseidon5
TARGETS_ONLYTEST += poseidon5_gadget
TARGETS_ONLYTEST += pow_gadget
TARGETS_ONLYTEST += sha256
TARGETS_ONLYTEST += sha256_gadget
TARGETS_ONLYTEST += sha512

# Targets which have tests and an additional executable (e.g. benchmarks)
TARGETS_TEST :=

# Targets which do not have tests
TARGETS_NOTEST :=
TARGETS_NOTEST += benchmark_mtree
#TARGETS_NOTEST += benchmark_abr

# Name of the library to build
LIBNAME := libzkp

# Put all targets together
TARGETS := $(TARGETS_ONLYTEST) $(TARGETS_TEST) $(TARGETS_NOTEST)
#### END TARGETS ####


#### BEGIN FILE SETUP ####
# Directory for binaries
BINPATH := bin
# Directory for object files
BUILDPATH := build
# Directory for include files
INCPATH := include
# Directory for library files
LIBPATH := lib
# Directory for source files
SRCPATH := src
# Directory for test files
TESTPATH := test

# Prefix for test binaries
TEST_PRE = test_

# Include files
INC := $(wildcard $(INCPATH)/*.h $(INCPATH)/**/*.h $(INCPATH)/*.hpp $(INCPATH)/**/*.hpp)
# Source files
SRC := $(wildcard $(SRCPATH)/*.cpp $(SRCPATH)/**/*.cpp)
# Test files
TEST := $(wildcard $(TESTPATH)/*.cpp $(TESTPATH)/**/*.cpp)

LIBNAMES := $(TARGETS_TEST) $(TARGETS_NOTEST)
SRCNAMES := $(notdir $(SRC))
TESTNAMES := $(notdir $(TEST))

# Object files
OBJ := $(SRCNAMES:%.cpp=$(BUILDPATH)/%.$(OEXT)) $(TESTNAMES:%.cpp=$(BUILDPATH)/$(TEST_PRE)%.$(OEXT))
# Library files
LIB := $(LIBNAMES:%=$(BUILDPATH)/%.$(OEXT))
# Dependency files (allows recompiling when a header changes without `make clean`)
DEP := $(OBJ:%.$(OEXT)=%.d) 

# Command to run all the tests at once, will be printed to a bash script
TEST_COMMAND := $(TARGETS_TEST:%=$$SCRIPTPATH/$(TEST_PRE)%;\n)$(TARGETS_ONLYTEST:%=$$SCRIPTPATH/$(TEST_PRE)%;\n)

#### END FILE SETUP ####


#### BEGIN COMPILER/LINKER SETUP ####
# If there is no default compiler, or if we want to override it, set the compiler here
ifeq ($(CXX), )
    CXX := c++
    CXXFLAGS := -Ofast -march=native
endif
ifeq ($(OVERRIDE_DEFAULT_CXX), 1)
    CXX := c++
    CXXFLAGS := -Ofast -march=native -flto=auto
endif

# Override flags for debugging
ifeq ($(DEBUG), 1)
    CXXFLAGS := -g -flto=auto
endif

# Additional compiler flags 
CXXFLAGS += -std=c++17 -Wall -fopenmp -I$(INCPATH) -MMD -MP

#### BEGIN MACRO DEFINITION ARGUMENTS ####
# Enable performance measurament code in test files
ifeq ($(MEASURE_PERFORMANCE), 1)
    CXXFLAGS += -DMEASURE_PERFORMANCE
endif

# Enable parallel code
ifeq ($(MULTICORE), 1)
    CXXFLAGS += -DMULTICORE
endif

# Enable assembly code
ifeq ($(USE_ASM), 1)
    CXXFLAGS += -DUSE_ASM
endif

# Additional macro definitions
CXXFLAGS += -D$(ELLIPTIC_CURVE) -DPERFORMANCE -DUSE_LIBFF
#### END MACRO DEFINITION ARGUMENTS ####

# Flags for c++20 targets
CXX20FLAGS := $(CXXFLAGS) -std=c++20

# Linker flags
LDFLAGS := -lff -lgmp -lgmpxx -lsnark -lprocps

# Archive command
AR := ar
# Archive flags
ARFLAGS := r
# Archive extension
LIBEXT := a

# Object file extenstion (some compilers use .obj)
OEXT := o
# Debug file extension (for Windows only)
PDB :=
#### END COMPILER/LINKER SETUP ####


#### BEGIN COMMAND ALIASING ####
# Make useful commands OS independent
MKDIR = mkdir -p $(1)
RM = rm -rf $(1)/*
MV = mv $(1) $(2)/
CP = cp $(1) $(2)/
#### END COMMAND ALIASING ####

# Build all targets
all: dirs testsuite $(TARGETS) library

# Remove all build files
clean:
	$(call RM,$(BUILDPATH))
	$(call RM,$(BINPATH))
	$(call RM,$(LIBPATH))

# Create build directories
dirs:
	$(call MKDIR,$(BUILDPATH))
	$(call MKDIR,$(BINPATH))
	$(call MKDIR,$(LIBPATH))

# Generate "run all tests" script
testsuite:
	@echo -e '#!/bin/sh\nSCRIPTPATH=$$(dirname "$$(readlink -f "$$0")")\n$(TEST_COMMAND)' > $(BINPATH)/test.sh

# Build library
library: $(TARGETS)
	$(AR) $(ARFLAGS) $(LIBOUTFLAGS)$(LIBPATH)/$(LIBNAME).$(LIBEXT) $(LIBS)

# Build targets which only have test files
$(TARGETS_ONLYTEST): %: $(BUILDPATH)/$(TEST_PRE)%.$(OEXT)
	$(CXX) $(CXXFLAGS) $^ -o $(BINPATH)/$(TEST_PRE)$@ $(LDFLAGS)

# Build targets which do not have test files
$(TARGETS_NOTEST):  %: $(BUILDPATH)/%.$(OEXT)
	$(CXX) $(CXXFLAGS) $^ -o $(BINPATH)/$@ $(LDFLAGS)

# Include additional dependencies to detect changes in headers
-include $(DEP)

# Compile tests with c++20 support (main directory)
$(BUILDPATH)/$(TEST_PRE)%.$(OEXT): $(TESTPATH)/cpp20/%.cpp
	$(CXX) $(CXX20FLAGS) -c $< -o $@

# Compile tests with c++20 support (recursive)
$(BUILDPATH)/$(TEST_PRE)%.$(OEXT): $(TESTPATH)/cpp20/**/%.cpp
	$(CXX) $(CXX20FLAGS) -c $< -o $@

# Compile tests (main directory)
$(BUILDPATH)/$(TEST_PRE)%.$(OEXT): $(TESTPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile tests (recursive)
$(BUILDPATH)/$(TEST_PRE)%.$(OEXT): $(TESTPATH)/**/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile sources with c++20 support (main directory)
$(BUILDPATH)/%.$(OEXT): $(SRCPATH)/cpp20/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile sources with c++20 support (recursive)
$(BUILDPATH)/%.$(OEXT): $(SRCPATH)/cpp20/**/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile sources (main directory)
$(BUILDPATH)/%.$(OEXT): $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile sources (recursive)
$(BUILDPATH)/%.$(OEXT): $(SRCPATH)/**/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
