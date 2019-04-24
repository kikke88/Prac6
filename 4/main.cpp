//  #include <iostream>
//  #include <cmath>
//  #include <cstdlib>

#include "quantum_logic_gate.hpp"
/*
Z =
Hadamard_gate 0                 --1
n_Hadamard_gate 1               --0
Phase_shift_gate 2              --1 (phi = const)
NOT_gate 3                      --1
CNOT_gate 4                     --2
Controlled_Phase_shift_gate 5   --2 (phi = const)
*/

int main(int argc, char* argv[]) {  //  n, Z, k, l, infile, outfile
    int n = strtol(argv[1], NULL, 10);
    int Z = strtol(argv[2], NULL, 10);
    int k = strtol(argv[3], NULL, 10);
    int l = strtol(argv[4], NULL, 10);
    double phi = M_PI;
    const char* infile;
    const char* outfile;
    infile = argv[5];
    outfile = argv[6];
    omp_set_num_threads(4);
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {
        Quantum_vector vec1(n, size, rank);
        vec1.vector_filling(infile);
        switch (Z) {
            case 0:
                vec1.Hadamard_gate(k);
                break;
            case 1:
                vec1.n_Hadamard_gate();
                break;
            case 2:
                vec1.Phase_shift_gate(k, phi);
                break;
            case 3:
                vec1.NOT_gate(k);
                break;
            case 4:
                vec1.CNOT_gate(k, l);
                break;
            case 5:
                vec1.Controlled_Phase_shift_gate(k, l, phi);
                break;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        vec1.file_output(outfile);
    }
    MPI_Finalize();
    return 0;
}
