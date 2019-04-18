#!/bin/bash
make

processes=4
qbits=5
bit1=1
bit2=2

./gen ${qbits} n_${qbits}.data
mpirun -n ${processes} ./main ${qbits} 4 ${bit1} ${bit2} n_${qbits}.data ${qbits}_cnot_${bit1}${bit2}.data
./cnot ${qbits} ${bit1} ${bit2} n_${qbits}.data cnot_RESULT.data
./compare ${qbits} ${qbits}_cnot_${bit1}${bit2}.data cnot_RESULT.data

make clean
