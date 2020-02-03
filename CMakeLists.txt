cmake_minimum_required(VERSION 3.15)
project(coda)

set(CMAKE_CXX_STANDARD 17)

# Flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories("include")

add_executable(coda coda.cpp mat.cpp utils.cpp element.cpp rsvd.cpp)
target_link_libraries(coda PUBLIC Eigen3::Eigen)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(coda PUBLIC OpenMP::OpenMP_CXX)
endif ()

