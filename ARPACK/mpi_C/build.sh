module purge
module load gcc/14 openmpi/4.1
module load arpack-mpi/3.9

mpicc -O3 -I$ARPACK_HOME/include/arpack main.c -o main -L$ARPACK_HOME/lib -lparpack -Wl,-rpath,$ARPACK_HOME/lib
