module purge
module load gcc/13
module load arpack-serial/3.9

gcc -O3 -I$ARPACK_HOME/include/arpack main.c -o main -L$ARPACK_HOME/lib -larpack -Wl,-rpath,$ARPACK_HOME/lib
