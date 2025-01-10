# Compiler and flags
CC = g++
CXXFLAGS = -fopenmp -I./include -L./libs
LDFLAGS = -lcdhitlib -lz

# Debugging flags (optional)
ifeq ($(debug),yes)
    CXXFLAGS += -ggdb
else
    CXXFLAGS += -O2
endif

# Source files
SRCS = main.cpp

# Output binary
TARGET = main

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Rule to build the main executable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET)

# Rule to compile .cpp to .o
%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJS) $(TARGET)

# Install target (optional)
install:
	install -m 0755 $(TARGET) $(PREFIX)/bin
