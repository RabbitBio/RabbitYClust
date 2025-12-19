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


# ---- edlib ----
# 查找 edlib 动态库（通过 -DEDLIB_ROOT 指定或自动查找）
if(EDLIB_ROOT)
    set(EDLIB_INCLUDE_DIR ${EDLIB_ROOT}/include)
    find_library(EDLIB_LIBRARY NAMES edlib libedlib.so PATHS ${EDLIB_ROOT} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
    find_path(EDLIB_INCLUDE_DIR edlib.h PATHS /usr/include /usr/local/include)
    find_library(EDLIB_LIBRARY NAMES edlib libedlib.so PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
endif()

if(NOT EDLIB_INCLUDE_DIR OR NOT EDLIB_LIBRARY)
    message(FATAL_ERROR "edlib not found. Please specify -DEDLIB_ROOT=<path> or install edlib system-wide.")
else()
    target_include_directories(yclust PUBLIC ${EDLIB_INCLUDE_DIR})
    target_link_libraries(yclust PRIVATE ${EDLIB_LIBRARY})
endif()
#-------------------

if(LIBCC_ROOT)
    set(LIBCC_LIB "${LIBCC_ROOT}/lib/libcdhitlib.so")
    set(LIBCC_INCLUDE_DIR "${LIBCC_ROOT}/include")
#message(STATUS "Using include path: ${LIBCC_INCLUDE_DIR}")
#get_filename_component(ABS_PATH ${LIBCC_INCLUDE_DIR} ABSOLUTE)
#message(STATUS "Absolute path: ${ABS_PATH}")
    # 检查头文件路径
    if(NOT EXISTS ${LIBCC_INCLUDE_DIR})
        message(FATAL_ERROR "cdhit include directory not found: ${LIBCC_INCLUDE_DIR}")
    endif()
    # 检查库文件
    if(NOT EXISTS ${LIBCC_LIB})
        message(FATAL_ERROR "cdhit library not found: ${LIBCC_LIB}")
    endif()
    target_include_directories(yclust PRIVATE ${LIBCC_INCLUDE_DIR})
    target_link_libraries(yclust PRIVATE ${LIBCC_LIB})
else()
    message(FATAL_ERROR "Please specify -DLIBCC_ROOT=<path>")
endif()

# 将 OpenMP 链接到目标
if(OpenMP_CXX_FOUND)
    target_link_libraries(yclust PUBLIC z OpenMP::OpenMP_CXX)
else()
	target_link_libraries(yclust z)
endif()

find_package(TBB CONFIG REQUIRED)
target_link_libraries(yclust PRIVATE TBB::tbb)
target_link_libraries(yclust PRIVATE atomic)
target_include_directories(yclust PRIVATE ${CMAKE_SOURCE_DIR}/include)

option(OPTIMIZE_NATIVE "为本机CPU优化（-march=native）" ON)

if(OPTIMIZE_NATIVE)
  target_compile_options(yclust PRIVATE -march=native)
  #target_compile_options(bench_main PRIVATE -O3 -DNDEBUG -march=native)
endif()
