#!/bin/bash
module purge
module load gcc/14 openmpi/4.1
module load petsc-real-double/3.22
module load slepc-real-double/3.22

mpicc -O3 -I$PETSC_HOME/include -I$SLEPC_HOME/include main.c -o main  -L$PETSC_HOME/lib -lpetsc -Wl,-rpath,$PETSC_HOME/lib -L$SLEPC_HOME/lib -lslepc -Wl,-rpath,$SLEPC_HOME/lib
 
