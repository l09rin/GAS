# -*- mode: makefile -*-

#==================================#
# Compiler Configurations
CXX = g++
DEFAULT_CXXFLAGS = -Wall -std=c++11 -lm
ADD_COMPILER_FLAGS ?=
CXXFLAGS := $(sort $(DEFAULT_CXXFLAGS) $(ADD_COMPILER_FLAGS))

NVCC = nvcc
ifndef CUDA_CC
	CUDA_CC := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
endif
NVCCFLAGS = -x cu -std=c++14 --diag-warn -Xcompiler "-Wall -Wextra" -gencode arch=compute_$(CUDA_CC),code=sm_$(CUDA_CC) -v
#==================================#

#==================================#
# Debug and Release Flags
DEBUG_FLAGS = -g -fno-inline
RELEASE_FLAGS = -O3
ifeq ($(DEBUG), YES)
	CXXFLAGS += $(DEBUG_FLAGS)
	NVCCFLAGS += -g -G
else
	CXXFLAGS += $(RELEASE_FLAGS)
	NVCCFLAGS += $(RELEASE_FLAGS)
endif
#==================================#

# Source File Path

_SRC = $(wildcard src/numerical_simulation/*.cpp src/numerical_simulation/*.ipp src/numerical_simulation/*.tpp src/numerical_simulation/Interactions/*.tpp src/numerical_simulation/Integrators/*.tpp)
SRC = $(filter-out %/polymer.tpp %/rings.tpp , $(_SRC))
SRC += $(wildcard src/numerical_simulation/*.cu src/numerical_simulation/Interactions/*.cu src/numerical_simulation/Integrators/*.cu)
HEADERS = $(filter-out %/polymer.h %/rings.h , $(wildcard src/numerical_simulation/*.h))
HEADERS += $(wildcard src/numerical_simulation/*.cuh src/numerical_simulation/Interactions/*.cuh)

SRC_MYLIBS = lib/matrix.cpp lib/record.cpp
OBJ_MYLIBS = $(patsubst lib/%.cpp,obj/%.o,$(SRC_MYLIBS))
DEPS_MYLIBS = lib/matrix.h lib/record.h lib/system_comm.h lib/system_comm.tpp
DEPS_MYLIBS += lib/gpu_comm.cuh

DEPS = $(SRC) $(HEADERS)
OBJ = obj/gas.o
EXE = gas.exe

all : test_gpu obj $(EXE)  ## Build GPU test, objects, and the main executable
	@echo "Built all the targets: "$(EXE)

$(OBJ_MYLIBS): obj/%.o: lib/%.cpp obj ${DEPS_MYLIBS}  # Compile library object files
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/gas.o : src/gas.cpp $(DEPS) $(SRC_MYLIBS) $(DEPS_MYLIBS) obj  # Compile CUDA gas object
	$(NVCC) $(NVCCFLAGS) -c src/gas.cpp -o $@

gas.exe : obj/gas.o $(OBJ_MYLIBS)  ## Link final executable
	$(NVCC) $(NVCCFLAGS) $^ -o $@

obj :  # Create object directory
	mkdir -p obj

clean-obj :  # Remove all object files
	rm -f $(OBJ) $(OBJ_MYLIBS)

clean : clean-obj  ## Remove objects and executables
	rm -f $(EXE)
	rm -f test/gpu/check_gpu.exe

test_gpu : lib/gpu_comm.cuh test/gpu/check_gpu.cu  ## Compile and run GPU test
	cd test/gpu/ ; \
	$(NVCC) $(NVCCFLAGS) check_gpu.cu -o check_gpu.exe ; \
	./check_gpu.exe
	@echo "GPU compilation test executed successfully!"

help:  ## Show this help
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+ *:.*?##' $(MAKEFILE_LIST) | \
	    awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

.PHONY : all clean clean-obj obj test_gpu help
