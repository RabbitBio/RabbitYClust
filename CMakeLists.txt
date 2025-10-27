cmake_minimum_required(VERSION 3.0)

project(RabbitYClust)

set(CMAKE_CXX_FLAGS "-std=c++17 -g -O3 -pthread ${CMAKE_CXX_FLAGS}")

if(USE_DEBUG)
	set(CMAKE_CXX_FLAGS "-D DEBUG ${CMAKE_CXX_FLAGS}")
endif()

if(USE_VERBOSE)
	set(CMAKE_CXX_FLAGS "-D VERBOSE ${CMAKE_CXX_FLAGS}")
endif()

if(USE_TIMING)
	set(CMAKE_CXX_FLAGS "-D TIMING ${CMAKE_CXX_FLAGS}")
endif() 


find_package(OpenMP REQUIRED)

aux_source_directory(src DIR_SRCS)
aux_source_directory(src/aahash DIR_SRCS)

add_executable(yclust ${DIR_SRCS})

# 启用线程池
option(ENABLE_THREAD_POOL "Enable thread pool support" OFF)
if(ENABLE_THREAD_POOL)
	message(STATUS "Thread pool enabled")
	target_compile_definitions(yclust PRIVATE THREAD_POOL)
endif()

set(LIBCDHIT "${CMAKE_SOURCE_DIR}/libs/libcdhitlib.so")
set(LIBCDHIT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

message(STATUS "Using include path: ${LIBCDHIT_INCLUDE_DIR}")
get_filename_component(ABS_PATH ${LIBCDHIT_INCLUDE_DIR} ABSOLUTE)
message(STATUS "Absolute path: ${ABS_PATH}")

# 检查头文件路径
if(NOT EXISTS ${LIBCDHIT_INCLUDE_DIR})
    message(FATAL_ERROR "cdhit include directory not found: ${LIBCDHIT_INCLUDE_DIR}")
endif()

# 检查库文件
if(NOT EXISTS ${LIBCDHIT})
    message(FATAL_ERROR "cdhit library not found: ${LIBCDHIT}")
endif()

target_include_directories(yclust PRIVATE ${LIBCDHIT_INCLUDE_DIR})
target_link_libraries(yclust PRIVATE ${LIBCDHIT})

# 将 OpenMP 链接到目标
if(OpenMP_CXX_FOUND)
    target_link_libraries(yclust PUBLIC z OpenMP::OpenMP_CXX)
else()
	target_link_libraries(yclust z)
endif()

