cmake_minimum_required(VERSION 3.0)

project(RabbitYClust)

set(CMAKE_CXX_FLAGS "-std=c++17 -g -O3 ${CMAKE_CXX_FLAGS}")

if(USE_DEBUG)
	set(CMAKE_CXX_FLAGS "-D DEBUG ${CMAKE_CXX_FLAGS}")
endif()

if(USE_Timer)
	set(CMAKE_CXX_FLAGS "-D Timer ${CMAKE_CXX_FLAGS}")
endif() 

#set(CMAKE_EXE_LINKER_FLAGS "-lz ${CMAKE_EXE_LINKER_FLAGS}")

include_directories(kernels src src/aahash)

#ahash 目录中的所有源文件
file(GLOB AAHASH_SOURCES src/aahash/*)


# 收集 src 目录下所有文件
file(GLOB ALL_SRCS "src/*")
# 排除 main.cpp
list(FILTER ALL_SRCS EXCLUDE REGEX "src/main\\.cpp$")
set(DIR_SRCS ${ALL_SRCS})
#aux_source_directory(src DIR_SRCS)

# 收集 src 目录下所有 .cpp 文件
file(GLOB DIR_KERNEL_SRCS "kernels/*.cpp")
# 排除 merge-2sets
list(FILTER DIR_KERNEL_SRCS EXCLUDE REGEX "kernels/merge-2sets\\.cpp$")
#aux_source_directory(kernels DIR_KERNEL_SRCS)

# 将源文件合并到一个变量中
list(APPEND DIR_SRCS ${DIR_KERNEL_SRCS} ${AAHASH_SOURCES})

add_executable(yclust ${DIR_SRCS})

target_link_libraries(yclust z)
