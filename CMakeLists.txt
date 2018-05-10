cmake_minimum_required(VERSION 3.9)
project(vdir2DHEParSudoku C)

set(CMAKE_C_STANDARD 11)

add_executable(vdir2DHEParSudoku vdir2DParSudoku.c task_def.c task_def.h)

set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")