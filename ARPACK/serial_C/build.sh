module purge
module load gcc/13
module load arpack-serial/3.9

#-march=skylake-avx512
gcc -O3 -g -I$ARPACK_HOME/include/arpack main.c -o main -L$ARPACK_HOME/lib -larpack -Wl,-rpath,$ARPACK_HOME/lib
