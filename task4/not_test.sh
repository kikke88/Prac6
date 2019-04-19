#!/bin/bash
make
processes=4
qbits=5
./gen ${qbits} n_${qbits}.data
for (( bit=1; bit <= ${qbits}; bit++ ))
do
    mpirun -n ${processes} ${1} ./main ${qbits} 3 ${bit} 0 n_${qbits}.data ${qbits}_not_${bit}.data
    ./not ${qbits} ${bit} n_${qbits}.data not_RESULT.data
    ./compare ${qbits} ${qbits}_not_${bit}.data not_RESULT.data
done
make clean
