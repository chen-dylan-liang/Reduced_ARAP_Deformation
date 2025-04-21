cmake_minimum_required(VERSION 3.29)
project(Reduced_ARAP_Deformation)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(Reduced_ARAP_Deformation
        src/ARAP.cpp
        src/ARAP.h
        src/call_back.h
        src/global_var.h
        src/helper_algebra.cpp
        src/helper_algebra.h
        src/helper_geometry.cpp
        src/helper_geometry.h
        src/helper_init.h
        src/helper_print.h
        src/image.cpp
        src/image.h
        src/main.cpp)
