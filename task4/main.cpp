//  #include <iostream>
//  #include <cmath>
//  #include <cstdlib>

#include "quantum_logic_gate.hpp"

int main(int argc, char* argv[]) {  //  n, 0 infile, 1 outfile
    int n = strtol(argv[1], NULL, 10);
    const char* infile;
    const char* outfile;
    if (argc > 2) {
        if (strtol(argv[2], NULL, 10) == 0) {
            infile = argv[3];
        } else {
            outfile = argv[3];
        }
    }
    if (argc > 4) {
        if (strtol(argv[4], NULL, 10) == 0) {
            infile = argv[5];
        } else {
            outfile = argv[5];
        }
    }
    omp_set_num_threads(1);
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {
        Quantum_vector vec1(n, size, rank);
        vec1.vector_filling(infile);
        vec1.CNOT_gate(1, 2);
        MPI_Barrier(MPI_COMM_WORLD);
        vec1.file_output(outfile);
    }
    MPI_Finalize();
    return 0;
}
