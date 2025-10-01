module purge
module load gcc/14
module load arpack-serial/3.9

gfortran -O3 main.f90 -o main -L$ARPACK_HOME/lib -larpack -Wl,-rpath,$ARPACK_HOME/lib
