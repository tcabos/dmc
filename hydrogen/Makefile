# Compiler and flags
CC = g++
CFLAGS = -std=c++20 -Wall -Wextra -fopenmp -O3 -march=native -ffast-math

# Directories
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = $(BIN_DIR)/obj

# Target executable name
TARGET = $(BIN_DIR)/hydrogen.app

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
