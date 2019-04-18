#!/bin/bash
make

processes=4
qbits=5
bit=3

./gen ${qbits} n_${qbits}.data
mpirun -n ${processes} ./main ${qbits} 3 ${bit} 0 n_${qbits}.data ${qbits}_not_${bit}.data
./not ${qbits} ${bit} n_${qbits}.data not_RESULT.data
./compare ${qbits} ${qbits}_not_${bit}.data not_RESULT.data

make clean
