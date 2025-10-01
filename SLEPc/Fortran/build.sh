#!/bin/bash
module purge
module load gcc/14 openmpi/4.1
module load petsc-real-double/3.22
module load slepc-real-double/3.22

mpif90 -O3 -I$PETSC_HOME/include -I$SLEPC_HOME/include ex1f.F90 -o ex1f  -L$PETSC_HOME/lib -lpetsc -Wl,-rpath,$PETSC_HOME/lib -L$SLEPC_HOME/lib -lslepc -Wl,-rpath,$SLEPC_HOME/lib
 
