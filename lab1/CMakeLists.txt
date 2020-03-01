cmake_minimum_required(VERSION 3.6)
project(lab1)

set(CMAKE_C_STANDARD 11)

add_executable(lab1 main.c Operations.h operations.c)
include_directories(/usr/include/mpi/)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
target_link_libraries(lab1 m)

