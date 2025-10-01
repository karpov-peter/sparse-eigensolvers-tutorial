module purge
module load gcc/14
module load arpack-serial/3.9

g++ -O3 -I$ARPACK_HOME/include/arpack main.cpp -o main -L$ARPACK_HOME/lib -larpack -Wl,-rpath,$ARPACK_HOME/lib
