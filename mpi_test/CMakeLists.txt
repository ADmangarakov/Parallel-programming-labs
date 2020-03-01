cmake_minimum_required(VERSION 3.15)
project(mpi_test C)

set(CMAKE_C_STANDARD 11)

add_executable(mpi_test main.c)
include_directories(/usr/include/mpi/)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
target_link_libraries(mpi_test m)