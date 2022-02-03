cmake_minimum_required(VERSION 3.21)
project(LUDecomposition C)

include_directories(/opt/homebrew/Cellar/mpich/4.0_1/include/)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -Werror")
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

add_executable(LUDecomposition main.c matrix.c matrix.h)
