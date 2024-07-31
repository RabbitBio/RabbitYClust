cmake_minimum_required(VERSION 3.0)

project(RabbitYClust)

set(CMAKE_CXX_FLAGS "-g -O3 ${CMAKE_CXX_FLAGS}")

if(USE_DEBUG)
	set(CMAKE_CXX_FLAGS "-D DEBUG ${CMAKE_CXX_FLAGS}")
endif()

if(USE_Timer)
	set(CMAKE_CXX_FLAGS "-D Timer ${CMAKE_CXX_FLAGS}")
endif() 

#set(CMAKE_EXE_LINKER_FLAGS "-lz ${CMAKE_EXE_LINKER_FLAGS}")


aux_source_directory(src DIR_SRCS)

add_executable(yclust ${DIR_SRCS})

target_link_libraries(yclust z)
