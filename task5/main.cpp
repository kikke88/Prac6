#include "quantum_logic_gate.hpp"

int main(int argc, char* argv[]) {  //  n, infile, outfile, [run_conf]
    int n = strtol(argv[1], NULL, 10);
    const char* infile = argv[2];
    const char* outfile = argv[3];
    const char* run_conf = argv[4];
    omp_set_num_threads(4);
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {
        Quantum_vector vec1(n, size, rank);
        vec1.vector_filling(infile);
        vec1.Fourier_transform(outfile, run_conf);
    }
    MPI_Finalize();
    return 0;
}
