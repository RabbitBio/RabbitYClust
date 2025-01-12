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

# 将 OpenMP 链接到目标
if(OpenMP_CXX_FOUND)
    target_link_libraries(yclust PUBLIC z OpenMP::OpenMP_CXX)
else()
	target_link_libraries(yclust z)
endif()
