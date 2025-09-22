# 编译器配置
CXX := g++
CC := gcc
CXXFLAGS := -std=c++17 -O3 -g -I. -I./aahash
CFLAGS := -O3 -g -I. -I./aahash
AR := ar
ARFLAGS := rcs
LDFLAGS := -lz -pthread -lstdc++fs

# 源文件搜索配置
SRC_DIR := .
AAHASH_DIR := aahash

# 递归查找所有源文件（.cpp 和 .c）
ALL_SRC := $(shell find $(SRC_DIR) $(AAHASH_DIR) -type f \( -name '*.cpp' -o -name '*.c' \))

# 排除 build_sketch.cpp
LIB_SRC := $(filter-out $(SRC_DIR)/build_sketch.cpp, $(ALL_SRC))

# 生成对象文件列表
LIB_OBJ := $(patsubst $(SRC_DIR)/%, obj/%, $(LIB_SRC:.c=.o))
LIB_OBJ := $(LIB_OBJ:.cpp=.o)

# 目标配置
LIB_NAME := libminhash.a
TARGET := build_sketch

# 默认目标
all: $(TARGET)

# 静态库
$(LIB_NAME): $(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $^

# 可执行文件
$(TARGET): build_sketch.o $(LIB_NAME)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# 主可执行文件的对象文件
build_sketch.o: build_sketch.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 模式规则：编译 C++ 文件
obj/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 模式规则：编译 C 文件
obj/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

# 清理
clean:
	rm -rf obj $(TARGET) $(LIB_NAME) build_sketch.o

	.PHONY: all clean
