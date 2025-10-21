# Compiler and flags
CC ?= gcc
CXX ?= g++
CFLAGS = -Wall -Wextra -O3 -std=c11 -MMD -MP
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -MMD -MP

# --- Experimental Evaluation params ---
ALPHA ?= 4
BETA ?= 24
PREPROCESSOR_VARS = -DALPHA=$(ALPHA) -DBETA=$(BETA)

ifeq ($(FRONTIER_DEBUG), 1)
PREPROCESSOR_VARS += -DFRONTIER_DEBUG
endif

ifeq ($(USE_PAPI), 1)
PREPROCESSOR_VARS += -DUSE_PAPI -lpapi -I${PAPI_DIR}/include -L${PAPI_DIR}/lib
endif

# --- Library Configuration ---
# Set the path to the root of the distributed_mmio library.
DIST_MMIO_PATH = ../distributed_mmio
LIB_STATIC_FULL_PATH = $(DIST_MMIO_PATH)/build/libdistributed_mmio.a

# CPPFLAGS: Pre-processor flags, primarily for include paths (-I).
CPPFLAGS = -I$(DIST_MMIO_PATH)/include -Iinclude
# LDFLAGS: Linker flags, primarily for library search paths (-L).
LDFLAGS = -L$(DIST_MMIO_PATH)/build

# LDLIBS: The libraries to link against.
LDLIBS = -ldistributed_mmio -lstdc++ -fopenmp

# --- Project Directories ---
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# --- Source, Object, and Dependency Files ---
SRCS =  $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))
DEPS = $(OBJS:.o=.d)

# The final executable target.
TARGET = $(BIN_DIR)/bfs

# Explicitly define the full path to the static library we depend on.
LIB_STATIC_FULL_PATH = $(DIST_MMIO_PATH)/build/libdistributed_mmio.a

# Print the GCC version
$(info GCC version: $(shell $(CXX) -dumpversion))

# Rules
all: $(TARGET)

$(LIB_STATIC_FULL_PATH):
	@echo "==> Configuring and building distributed_mmio library..."
	@mkdir -p $(DIST_MMIO_PATH)/build
	@cd $(DIST_MMIO_PATH)/build && cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) ..
	@$(MAKE) -C $(DIST_MMIO_PATH)/build

# Rule to link the final executable.
# It depends on all object files AND the static library itself.
$(TARGET): $(OBJS) $(LIB_STATIC_FULL_PATH)
	@echo "==> Linking objects with the distributed_mmio library..."
	@mkdir -p $(BIN_DIR)
	$(CXX) -fopenmp -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(PREPROCESSOR_VARS)
	@echo "==> Build successful: $(TARGET)"

# Pattern rule to compile a .cpp file into a .o file.
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "==> Compiling: $<"
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(PREPROCESSOR_VARS) -fopenmp -c $< -o $@

# Include auto-generated dependency files if they exist
-include $(DEPS)

# Rule to clean up all generated files.
.PHONY: clean
clean:
	@echo "==> Cleaning up build artifacts..."
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	@echo "==> Cleaning up distributed_mmio library..."
	@rm -rf $(DIST_MMIO_PATH)/build
	@echo "==> Cleanup complete."

.PHONY: all clean

# --- Debugging ---
# To build for debugging: make debug
# This will clean the project and then rebuild with debugging symbols.
.PHONY: debug
debug:
	$(MAKE) clean
	$(MAKE) all CXXFLAGS="$(filter-out -O3,$(CXXFLAGS)) -g -O0"