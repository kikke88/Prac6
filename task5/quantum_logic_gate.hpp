#ifndef TASK5_QUANTUM_LOGIC_GATE_HPP_
#define TASK5_QUANTUM_LOGIC_GATE_HPP_

#include <mpi.h>
#include <omp.h>
#include <complex>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

typedef std::complex<double> complexd;

class Quantum_vector {
 private:
    const int num_of_qubits;
    const int num_of_elem;
    const int size;
    const int rank;
    const int elem_in_one_proc;
    complexd* vector_1;
    complexd* vector_2;
    MPI_Datatype double_double;
    void file_manipulation(const char*, const int);
    void single_qubit_transform(const int, const complexd[4]);
    void double_qubits_transform(const int, const int, const complexd[16]);

 public:
    Quantum_vector(const int, const int, const int);
    ~Quantum_vector();
    void vector_filling();
    void vector_filling(const char*);
    void Hadamard_gate(const int);
    void n_Hadamard_gate();
    void Phase_shift_gate(const int, const double);
    void NOT_gate(const int);
    void CNOT_gate(const int, const int);
    void Controlled_Phase_shift_gate(const int, const int, const double);
    void Fourier_transform(const char*, const char*);
    void file_output(const char*);
};

Quantum_vector::Quantum_vector(const int _num_of_qubits,
                               const int _size,
                               const int _rank):
    num_of_qubits {_num_of_qubits},
    num_of_elem {1 << num_of_qubits},
    size {_size},
    rank {_rank},
    elem_in_one_proc {num_of_elem / size},
    vector_1 {new complexd[elem_in_one_proc]},
    vector_2 {new complexd[elem_in_one_proc]}
{
    MPI_Type_contiguous(2, MPI_DOUBLE, &double_double);
    MPI_Type_commit(&double_double);
}

Quantum_vector::~Quantum_vector() {
    delete[] vector_1;
    delete[] vector_2;
    MPI_Type_free(&double_double);
    //  !!!
    //  MPI_Finalize();
    //  !!!
}

void Quantum_vector::file_manipulation(const char* filename, const int flag) {
    MPI_Datatype subarr_type;
    MPI_File file;
    int start {elem_in_one_proc * rank};
    MPI_Type_create_subarray(1, &num_of_elem, &elem_in_one_proc, &start,
                             MPI_ORDER_C, double_double, &subarr_type);
    MPI_Type_commit(&subarr_type);
    if (flag) {
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,
                      MPI_INFO_NULL, &file);
        MPI_File_set_view(file, 0, double_double, subarr_type,
                          "native", MPI_INFO_NULL);
        MPI_File_read(file, vector_1, elem_in_one_proc,
                          double_double, MPI_STATUS_IGNORE);
    } else {
        MPI_File_open(MPI_COMM_WORLD, filename,
                      MPI_MODE_WRONLY | MPI_MODE_CREATE,
                      MPI_INFO_NULL, &file);
        MPI_File_set_view(file, 0, double_double, subarr_type,
                          "native", MPI_INFO_NULL);
        MPI_File_write(file, vector_1, elem_in_one_proc,
                           double_double, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&subarr_type);
}

void Quantum_vector::vector_filling() {
    unsigned int cur_time;
    if (rank == 0) {
        cur_time = time(NULL);
    }
    if (size != 1) {
        MPI_Bcast(&cur_time, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }
    double sum = 0;
    const int tmp_const = RAND_MAX / 2;
#pragma omp parallel
    {
        unsigned int seed = cur_time +
                            rank * omp_get_num_threads() +
                            omp_get_thread_num();
#pragma omp for reduction(+:sum)
        for (int i = 0; i < elem_in_one_proc; ++i) {
            vector_1[i] = complexd(rand_r(&seed) - tmp_const,
                                   rand_r(&seed) - tmp_const);
        //    vector_1[i] = complexd(i, i + 1);
            sum += norm(vector_1[i]);
        }
    }
    if (size != 1) {
        MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
    }
    sum = sqrt(sum);
#pragma omp for
    for (int i = 0; i < elem_in_one_proc; ++i) {
        vector_1[i] /= sum;
    }
}

void Quantum_vector::vector_filling(const char* input_file_name) {
    int input_flag {1};
    file_manipulation(input_file_name, input_flag);
}

void Quantum_vector::single_qubit_transform(const int k,
                                            const complexd matrix[4]) {
    int test {1 << (num_of_qubits - k)};
    int log_size {static_cast<int>(log2(size))};
    if (k > log_size) {
#pragma omp parallel for
        for (int i = 0; i < elem_in_one_proc; ++i) {
            vector_2[i] = matrix[(i & test && 1) * 2] * vector_1[i & ~(test)] +
                          matrix[(i & test && 1) * 2 + 1] * vector_1[i | test];
        }
    } else {
        complexd* recv_vec {new complexd[elem_in_one_proc]};
        int dest_rank {rank ^ 1 << (log_size - k)};
        int swap_bit {rank & 1 << (log_size - k) && 1};
        if (swap_bit == 1) {
            MPI_Send(vector_1, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD);
            MPI_Recv(recv_vec, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            swap(vector_1, recv_vec);
        } else {
            MPI_Recv(recv_vec, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(vector_1, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD);
        }
#pragma omp parallel for
        for (int i = 0; i < elem_in_one_proc; ++i) {
            vector_2[i] = matrix[swap_bit * 2] * vector_1[i] +
                          matrix[swap_bit * 2 + 1] * recv_vec[i];
        }
        if (swap_bit == 1) {
            swap(vector_1, recv_vec);
        }
        delete[] recv_vec;
    }
    swap(vector_1, vector_2);
    return;
}

void Quantum_vector::Hadamard_gate(const int k) {
    const complexd tmp(1 / sqrt(2), 0);
    const complexd H[4] {tmp, tmp, tmp, -tmp};
    single_qubit_transform(k, H);
}

void Quantum_vector::n_Hadamard_gate() {
    for (int j {1}; j < num_of_qubits + 1; ++j) {
        Hadamard_gate(j);
    }
}

void Quantum_vector::Phase_shift_gate(const int k, const double phi) {
    const complexd P_S[4] {complexd(1, 0), complexd(0, 0),
                           complexd(0, 0), complexd(cos(phi), sin(phi))};
    single_qubit_transform(k, P_S);
}

void Quantum_vector::NOT_gate(const int k) {
    const complexd NOT[4] {complexd(0, 0), complexd(1, 0),
                           complexd(1, 0), complexd(0, 0)};
    single_qubit_transform(k, NOT);
}

void Quantum_vector::double_qubits_transform(const int k,
                                             const int l,
                                             const complexd matrix[16]) {
    int log_size {static_cast<int>(log2(size))};
    int test_k {1 << (num_of_qubits - k)};
    int test_l {1 << (num_of_qubits - l)};
    if (k > log_size && l > log_size) {
        int i_k_i_l;
#pragma omp parallel for private(i_k_i_l)
        for (int i = 0; i < elem_in_one_proc; ++i) {
            i_k_i_l = (i & test_k && 1) << 1 | (i & test_l && 1);
            vector_2[i] = matrix[i_k_i_l] *
                          vector_1[i & ~(test_k | test_l)] +
                          matrix[4 | i_k_i_l] *
                          vector_1[(i & ~test_k) | test_l] +
                          matrix[8 | i_k_i_l] *
                          vector_1[(i & ~test_l) | test_k] +
                          matrix[12 | i_k_i_l] *
                          vector_1[i | test_k | test_l];
        }
    } else if (k > log_size && l <= log_size) {
        complexd* recv_vec {new complexd[elem_in_one_proc]};
        int dest_rank {rank ^ 1 << (log_size - l)};
        int swap_bit {rank & 1 << (log_size - l) && 1};
        if (swap_bit == 1) {
            MPI_Send(vector_1, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD);
            MPI_Recv(recv_vec, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            swap(vector_1, recv_vec);
        } else {
            MPI_Recv(recv_vec, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(vector_1, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD);
        }
        int i_k_i_l;
#pragma omp parallel for private(i_k_i_l)
        for (int i = 0; i < elem_in_one_proc; ++i) {
            i_k_i_l = (i & test_k && 1) << 1 | swap_bit;
            vector_2[i] = matrix[i_k_i_l] * vector_1[i & ~test_k] +
                          matrix[4 | i_k_i_l] * recv_vec[i & ~test_k] +
                          matrix[8 | i_k_i_l] * vector_1[i | test_k] +
                          matrix[12 | i_k_i_l] * recv_vec[i | test_k];
        }
        if (swap_bit == 1) {
            swap(vector_1, recv_vec);
        }
        delete[] recv_vec;
    } else if (k <= log_size && l > log_size) {
        complexd* recv_vec {new complexd[elem_in_one_proc]};
        int dest_rank {rank ^ 1 << (log_size - k)};
        int swap_bit {rank & 1 << (log_size - k) && 1};
        if (swap_bit == 1) {
            MPI_Send(vector_1, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD);
            MPI_Recv(recv_vec, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            swap(vector_1, recv_vec);
        } else {
            MPI_Recv(recv_vec, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(vector_1, elem_in_one_proc, double_double, dest_rank,
                     0, MPI_COMM_WORLD);
        }
        int i_k_i_l;
#pragma omp parallel for private(i_k_i_l)
        for (int i = 0; i < elem_in_one_proc; ++i) {
            i_k_i_l = swap_bit << 1 | (i & test_l && 1);
            vector_2[i] = matrix[i_k_i_l] * vector_1[i & ~test_l] +
                          matrix[4 | i_k_i_l] * vector_1[i | test_l] +
                          matrix[8 | i_k_i_l] * recv_vec[i & ~test_l] +
                          matrix[12 | i_k_i_l] * recv_vec[i | test_l];
        }
        if (swap_bit == 1) {
            swap(vector_1, recv_vec);
        }
        delete[] recv_vec;
    } else {
        complexd* recv_vec_0 {new complexd[elem_in_one_proc]};
        complexd* recv_vec_1 {new complexd[elem_in_one_proc]};
        complexd* recv_vec_2 {new complexd[elem_in_one_proc]};
        int place_bits_k {rank & 1 << (log_size - k) && 1};
        int place_bits_l {rank & 1 << (log_size - l) && 1};
        MPI_Group WORLD_GROUP;
        MPI_Group NEW_GROUP;
        MPI_Comm NEW_COMM;
        MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GROUP);
        int i_k_i_l;
        if (place_bits_k == 0 && place_bits_l == 0) {
            i_k_i_l = 0;
            int dest_rank_0 {rank ^ 1 << (log_size - l)};
            int dest_rank_1 {rank ^ 1 << (log_size - k)};
            int dest_rank_2 {dest_rank_0 ^ 1 << (log_size - k)};
            const int ranks[4] {rank, dest_rank_0, dest_rank_1, dest_rank_2};
            MPI_Group_incl(WORLD_GROUP, 4, ranks, &NEW_GROUP);
            MPI_Comm_create(MPI_COMM_WORLD, NEW_GROUP, &NEW_COMM);
            MPI_Bcast(vector_1, elem_in_one_proc, double_double, 0, NEW_COMM);
            MPI_Bcast(recv_vec_0, elem_in_one_proc, double_double, 1, NEW_COMM);
            MPI_Bcast(recv_vec_1, elem_in_one_proc, double_double, 2, NEW_COMM);
            MPI_Bcast(recv_vec_2, elem_in_one_proc, double_double, 3, NEW_COMM);
#pragma omp parallel for
            for (int i = 0; i < elem_in_one_proc; ++i) {
                vector_2[i] = matrix[i_k_i_l] * vector_1[i] +
                              matrix[4 | i_k_i_l] * recv_vec_0[i] +
                              matrix[8 | i_k_i_l] * recv_vec_1[i] +
                              matrix[12 | i_k_i_l] * recv_vec_2[i];
        }
        } else if (place_bits_k == 0 && place_bits_l == 1) {
            i_k_i_l = 1;
            int dest_rank_0 {rank ^ 1 << (log_size - l)};
            int dest_rank_1 {dest_rank_0 ^ 1 << (log_size - k)};
            int dest_rank_2 {rank ^ 1 << (log_size - k)};
            const int ranks[4] {dest_rank_0, rank, dest_rank_1, dest_rank_2};
            MPI_Group_incl(WORLD_GROUP, 4, ranks, &NEW_GROUP);
            MPI_Comm_create(MPI_COMM_WORLD, NEW_GROUP, &NEW_COMM);
            MPI_Bcast(recv_vec_0, elem_in_one_proc, double_double, 0, NEW_COMM);
            MPI_Bcast(vector_1, elem_in_one_proc, double_double, 1, NEW_COMM);
            MPI_Bcast(recv_vec_1, elem_in_one_proc, double_double, 2, NEW_COMM);
            MPI_Bcast(recv_vec_2, elem_in_one_proc, double_double, 3, NEW_COMM);
#pragma omp parallel for
            for (int i = 0; i < elem_in_one_proc; ++i) {
                vector_2[i] = matrix[i_k_i_l] * recv_vec_0[i] +
                              matrix[4 | i_k_i_l] * vector_1[i] +
                              matrix[8 | i_k_i_l] * recv_vec_1[i] +
                              matrix[12 | i_k_i_l] * recv_vec_2[i];
            }
        } else if (place_bits_k == 1 && place_bits_l == 0) {
            i_k_i_l = 2;
            int dest_rank_0 {rank ^ 1 << (log_size - k)};
            int dest_rank_1 {dest_rank_0 ^ 1 << (log_size - l)};
            int dest_rank_2 {rank ^ 1 << (log_size - l)};
            const int ranks[4] {dest_rank_0, dest_rank_1, rank, dest_rank_2};
            MPI_Group_incl(WORLD_GROUP, 4, ranks, &NEW_GROUP);
            MPI_Comm_create(MPI_COMM_WORLD, NEW_GROUP, &NEW_COMM);
            MPI_Bcast(recv_vec_0, elem_in_one_proc, double_double, 0, NEW_COMM);
            MPI_Bcast(recv_vec_1, elem_in_one_proc, double_double, 1, NEW_COMM);
            MPI_Bcast(vector_1, elem_in_one_proc, double_double, 2, NEW_COMM);
            MPI_Bcast(recv_vec_2, elem_in_one_proc, double_double, 3, NEW_COMM);
#pragma omp parallel for
            for (int i = 0; i < elem_in_one_proc; ++i) {
                vector_2[i] = matrix[i_k_i_l] * recv_vec_0[i] +
                              matrix[4 | i_k_i_l] * recv_vec_1[i] +
                              matrix[8 | i_k_i_l] * vector_1[i] +
                              matrix[12 | i_k_i_l] * recv_vec_2[i];
        }
        } else {
            i_k_i_l = 3;
            int dest_rank_0 {(rank ^ 1 << (log_size - l))
                             ^ 1 << (log_size - k)};
            int dest_rank_1 {rank ^ 1 << (log_size - k)};
            int dest_rank_2 {rank ^ 1 << (log_size - l)};
            const int ranks[4] {dest_rank_0, dest_rank_1, dest_rank_2, rank};
            MPI_Group_incl(WORLD_GROUP, 4, ranks, &NEW_GROUP);
            MPI_Comm_create(MPI_COMM_WORLD, NEW_GROUP, &NEW_COMM);
            MPI_Bcast(recv_vec_0, elem_in_one_proc, double_double, 0, NEW_COMM);
            MPI_Bcast(recv_vec_1, elem_in_one_proc, double_double, 1, NEW_COMM);
            MPI_Bcast(recv_vec_2, elem_in_one_proc, double_double, 2, NEW_COMM);
            MPI_Bcast(vector_1, elem_in_one_proc, double_double, 3, NEW_COMM);
#pragma omp parallel for
            for (int i = 0; i < elem_in_one_proc; ++i) {
                vector_2[i] = matrix[i_k_i_l] * recv_vec_0[i] +
                              matrix[4 | i_k_i_l] * recv_vec_1[i] +
                              matrix[8 | i_k_i_l] * recv_vec_2[i] +
                              matrix[12 | i_k_i_l] * vector_1[i];
            }
        }
        delete[] recv_vec_0;
        delete[] recv_vec_1;
        delete[] recv_vec_2;
        MPI_Comm_free(&NEW_COMM);
        MPI_Group_free(&NEW_GROUP);
        MPI_Group_free(&WORLD_GROUP);
    }
    swap(vector_1, vector_2);
    return;
}

void Quantum_vector::CNOT_gate(const int k, const int l) {
    const complexd CNOT[16] {
        complexd(1, 0), complexd(0, 0), complexd(0, 0), complexd(0, 0),
        complexd(0, 0), complexd(1, 0), complexd(0, 0), complexd(0, 0),
        complexd(0, 0), complexd(0, 0), complexd(0, 0), complexd(1, 0),
        complexd(0, 0), complexd(0, 0), complexd(1, 0), complexd(0, 0)};
    double_qubits_transform(k, l, CNOT);
}

void Quantum_vector::Controlled_Phase_shift_gate(const int k,
                                                 const int l,
                                                 const double phi) {
    const complexd C_P_S[16] {
        complexd(1, 0), complexd(0, 0), complexd(0, 0), complexd(0, 0),
        complexd(0, 0), complexd(1, 0), complexd(0, 0), complexd(0, 0),
        complexd(0, 0), complexd(0, 0), complexd(1, 0), complexd(0, 0),
        complexd(0, 0), complexd(0, 0), complexd(0, 0),
                                        complexd(cos(phi), sin(phi))};
    double_qubits_transform(k, l, C_P_S);
}

void Quantum_vector::file_output(const char* output_file_name) {
    int output_flag {0};
    file_manipulation(output_file_name, output_flag);
}

void Quantum_vector::Fourier_transform(const char* outfile_name,
                                        const char* run_conf) {
    double compute_time = 0.0;
    compute_time -= MPI_Wtime();
    for (int i = 1; i <= num_of_qubits; ++i) {
        Hadamard_gate(i);
        int k = 2;
        for (int j = i + 1; j <= num_of_qubits; ++j) {
            Controlled_Phase_shift_gate(j, i, 2 * M_PI / (1 << k++));
        }
    }
    compute_time += MPI_Wtime();
    if (rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, &compute_time, 1, MPI_DOUBLE,
                    MPI_MAX, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(&compute_time, NULL, 1, MPI_DOUBLE, MPI_MAX,
                    0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        std::ofstream ofile("time_file", std::ios::app);
        ofile << run_conf << "  " << compute_time << std::endl;
    }

    std::vector<int> bit_masks(num_of_qubits);
    bit_masks[num_of_qubits - 1] = 1;
    for (int i = num_of_qubits - 2; i >= 0; --i) {
        bit_masks[i] = bit_masks[i + 1] << 1;
    }
    int min_index = size + 1;
    std::vector<int> index_vec(elem_in_one_proc);
    for (int i = 0 ; i < elem_in_one_proc; ++i) {
        int new_index = 0;
        for (int j = 0; j < num_of_qubits; ++j) {
            new_index |= (((i + rank * elem_in_one_proc) & bit_masks[j]) && 1)
                                                                    << j;
        }
        if (new_index < min_index) {
            min_index = new_index;
        }
        index_vec[i] = new_index;
    }
    for (int i = 0; i < elem_in_one_proc; ++i) {
        index_vec[i] = (index_vec[i] - min_index) / size;
        vector_2[index_vec[i]] = vector_1[i];
    }
    swap(vector_1, vector_2);
    MPI_Datatype subarr_type;
    MPI_File file;
    int elem = 1;
    MPI_File_open(MPI_COMM_WORLD, outfile_name,
                  MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &file);
    MPI_Type_create_subarray(1, &size, &elem, &min_index, MPI_ORDER_C,
                             double_double, &subarr_type);
    MPI_Type_commit(&subarr_type);
    MPI_File_set_view(file, 0, double_double, subarr_type, "native",
                      MPI_INFO_NULL);
    MPI_File_write(file, vector_1, elem_in_one_proc, double_double,
                   MPI_STATUS_IGNORE);
    MPI_Type_free(&subarr_type);
}

#endif  // TASK5_QUANTUM_LOGIC_GATE_HPP_
