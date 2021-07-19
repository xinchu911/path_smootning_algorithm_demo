cmake_minimum_required(VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

project(Spline_generation)

find_package(OsqpEigen)
find_package(Eigen3)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

#MPCExample
add_executable(main src/main.cpp)
target_link_libraries(main OsqpEigen::OsqpEigen)