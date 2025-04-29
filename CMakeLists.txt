cmake_minimum_required(VERSION 3.29)
project(Reduced_ARAP_Deformation)

set(CMAKE_CXX_STANDARD 17)

# OpenGL, GLFW, OpenCV
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(OpenCV REQUIRED)

# Eigen Setup
# add_subdirectory(external/eigen)
include_directories(
        external/eigen
        external/glad/include
        external/libigl/include
        external/opencv/modules
        ${OpenCV_INCLUDE_DIRS}
        src
)

# glad
add_library(glad STATIC external/glad/src/glad.c)

# Executable
add_executable(Reduced_ARAP_Deformation
        src/energy.h
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
        src/main.cpp
)

target_link_libraries(Reduced_ARAP_Deformation
        glad
        glfw
        ${OpenCV_LIBS}
        ${OPENGL_LIBRARIES}
)
