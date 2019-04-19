#!/bin/bash
make

processes=4
qbits=5

./gen ${qbits} n_${qbits}.data
for (( bit1=1; bit1 <= ${qbits}; bit1++ ))
do
    for (( bit2=1; bit2 <= ${qbits}; bit2++ ))
    do
        if [ ${bit1} -eq ${bit2} ]
        then
            continue
        fi
        mpirun -n ${processes} ./main ${qbits} 4 ${bit1} ${bit2} n_${qbits}.data ${qbits}_cnot_${bit1}${bit2}.data
        ./cnot ${qbits} ${bit1} ${bit2} n_${qbits}.data cnot_RESULT.data
        ./compare ${qbits} ${qbits}_cnot_${bit1}${bit2}.data cnot_RESULT.data
    done
done

make clean
