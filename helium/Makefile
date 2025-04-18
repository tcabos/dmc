# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -static-libstdc++ -static-libgcc -O3 -Wall -Wextra -ffast-math -fopenmp

ifdef WRITE_STEP_ENERGY_FILE
	CFLAGS += -DWRITE_STEP_ENERGY_FILE
endif

ifdef WRITE_BLOCK_FILE
	CFLAGS += -DWRITE_BLOCK_FILE
endif

ifdef WRITE_CUMUL_BLOCK_FILES
	CFLAGS += -DWRITE_CUMUL_BLOCK_FILES
endif

ifdef HISTOGRAM
	CFLAGS += -DHISTOGRAM
endif

ifdef SNAPSHOT_WALKERS_START_END
	CFLAGS += -DSNAPSHOT_WALKERS_START_END
endif

ifdef PRINT_BLOCK_DATA_TO_CONSOLE
	CFLAGS += -DPRINT_BLOCK_DATA_TO_CONSOLE
endif

ifdef SNAPSHOT_WALKERS_EQUILIBRATION
	CFLAGS += -DSNAPSHOT_WALKERS_EQUILIBRATION
endif

ifdef SNAPSHOT_WALKERS_PRODUCTION
	CFLAGS += -DSNAPSHOT_WALKERS_PRODUCTION
endif

# Directories
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = $(BIN_DIR)/obj

# Target executable name
TARGET = $(BIN_DIR)/helium.app

# Source and object files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(OBJECTS) -o $(TARGET) $(CFLAGS)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create directories if they don't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean up
clean:
	rm -rf $(BIN_DIR)

# Phony targets
.PHONY: all clean
