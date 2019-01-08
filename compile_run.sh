#!/usr/bin/env bash
mpicc -o vdir2DParSudoku.o vdir2DParSudoku.c task_def.c -Wall -Wextra -g -lm
mpirun -np 4 ./vdir2DParSudoku.o 100 60 100
