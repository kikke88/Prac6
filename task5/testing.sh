make

mpirun -n 2 ./main 2 2cub_check.in out2.data
./compare 2 2cub_check.res out2.data

mpirun -n 2 ./main 3 3cub_check.in out3.data
./compare 3 3cub_check.res out3.data

make clean
