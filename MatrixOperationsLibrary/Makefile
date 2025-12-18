# Makefile for the Matrix Market linear solver (single-file version)

# Compiler
CXX = g++

# Compiler flags
# -std=c++20 : Required for proper template and filesystem usage
# -Wall -Wextra : Enable useful warnings
# -O2       : Optimize for performance
# -g        : Include debugging symbols (remove for pure release)
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -g

# Name of the final executable
TARGET = solver

# Source file (single file version)
SOURCES = main.cpp

# Object files (automatically derived)
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets (not actual files)
.PHONY: all clean

# Run the program (convenience target)
# Assumes A.dat and B.dat exist in the current directory
run: $(TARGET)
	./$(TARGET)

# Rebuild everything from scratch
rebuild: clean all
