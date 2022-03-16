cmake_minimum_required(VERSION 3.5)

project(compute LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
include_directories(
    /opt/xilinx/xrt/include
    /usr/lib/jvm/java-11-openjdk-amd64/include/linux
    /usr/lib/jvm/java-11-openjdk-amd64/include/
    )
link_directories(
    /opt/xilinx/xrt/lib
    )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__USE_XOPEN2K8")
set(CMAKE_EXE_LINKER_FLAGS "")

set(Sources
    compute_cl.hpp tb.cpp gatk.cpp common.h
    )
add_executable(benchmark ${Sources})
add_library(compute SHARED ${Sources} jni.cpp)
target_link_libraries(benchmark xilinxopencl)
target_link_libraries(compute xilinxopencl)
