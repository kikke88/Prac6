#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include "mpi.h"

#include <omp.h>

#define GEN 0
#define FILE 1

using namespace std;
typedef complex<double> complexd;

void quantum_transformation(complexd* a, complexd* result, const int n, const double u[][2],
							const int k, const int size, const int rank, const MPI_Datatype double_double)
{
	int num_of_elem = 1 << n;
	int test = 1 << (n - k);
	int log_k = log2(size);
	if (k > log_k) {
		for (int i = 0; i < num_of_elem / size; ++i) {
			result[i] = u[i & test && 1][0] * a[i & ~(test)] + u[i & test && 1][1] * a[i | test];
		}
	} else {
		complexd* recv_vec = new complexd[num_of_elem / size];
		int dest_rank = rank ^ 1 << (log_k - k);
		int swap_bit = rank & 1 << (log_k - k) && 1;
		if (swap_bit == 1) {
			MPI_Send(a, num_of_elem / size, double_double, dest_rank, 0 , MPI_COMM_WORLD);
			MPI_Recv(recv_vec, num_of_elem / size, double_double, dest_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			complexd* tmp = a;
			a = recv_vec;
			recv_vec = tmp;
		} else {
			MPI_Recv(recv_vec, num_of_elem / size, double_double, dest_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(a, num_of_elem / size, double_double, dest_rank, 0 , MPI_COMM_WORLD);
		}
		//MPI_Sendrecv(a, num_of_elem / size, double_double, dest_rank, 0,
		//			recv_vec, num_of_elem / size, double_double, dest_rank, 0,
		//			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int i = 0; i < num_of_elem / size; ++i) {
			result[i] = u[swap_bit][0] * a[i] + u[swap_bit][1] * recv_vec[i];
		}
		if (swap_bit == 1) {
			recv_vec = a;
		}
		delete[] recv_vec;
	}
	return;
}

int main(int argc, char* argv[])//n, k, генерируешь сам 0 / файл 1, in_file, out_file
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = strtol(argv[1], NULL, 10),
		k = strtol(argv[2], NULL, 10),
		num_of_elem = 1 << n,
        mode = strtol(argv[3], NULL, 10);
	double compute_time;
	complexd* qbit_vec_local = new complexd[num_of_elem / size];
    MPI_Datatype subarr_type;
	MPI_Datatype double_double;
	MPI_Type_contiguous(2, MPI_DOUBLE, &double_double);
	MPI_Type_commit(&double_double);
	int array_of_subsizes = num_of_elem / size;
    if (mode == GEN) {
        unsigned int cur_time;
        if (rank == 0) {
            cur_time = time(NULL);
        }
		if (size != 1) {
			MPI_Bcast(&cur_time, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		}
        cur_time += rank;
    	double sum_local = 0;
    	for (int i = 0; i < num_of_elem / size; ++i) {
    		qbit_vec_local[i] = complexd(rand_r(&cur_time), rand_r(&cur_time));
    		sum_local += norm(qbit_vec_local[i]);
    	}
		if (size != 1) {
			MPI_Allreduce(MPI_IN_PLACE, &sum_local, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		}
        sum_local = sqrt(sum_local);
    	for (int i = 0; i < num_of_elem / size; ++i) {
    		qbit_vec_local[i] /= sum_local;
    	}
    } else {
        MPI_File in_file;
        int start = array_of_subsizes * rank;
        MPI_Type_create_subarray(1, &num_of_elem, &array_of_subsizes, &start, MPI_ORDER_C, double_double, &subarr_type);
        MPI_Type_commit(&subarr_type);
        MPI_File_open(MPI_COMM_WORLD, argv[4], MPI_MODE_RDONLY, MPI_INFO_NULL, &in_file);
        MPI_File_set_view(in_file, 0, double_double, subarr_type, "native", MPI_INFO_NULL);
        MPI_File_read_all(in_file, qbit_vec_local, array_of_subsizes, double_double, MPI_STATUS_IGNORE);
	}
	complexd* result_local = new complexd[num_of_elem / size];
	double tmp = 1 / sqrt(2);
	double matrix[2][2] = {{tmp, tmp}, {tmp, -tmp}};
	MPI_Barrier(MPI_COMM_WORLD);
	compute_time -= MPI_Wtime();
	quantum_transformation(qbit_vec_local, result_local, n, matrix, k, size, rank, double_double);
	compute_time += MPI_Wtime();
	if (mode == FILE) {

		MPI_File out_file;
        MPI_File_open(MPI_COMM_WORLD, argv[5], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out_file);
        MPI_File_set_view(out_file, 0, double_double, subarr_type, "native", MPI_INFO_NULL);
		MPI_File_write_all(out_file, result_local, array_of_subsizes, double_double, MPI_STATUS_IGNORE);

	} else {
/*
//вывод, случай параллельной генерации
		if (rank == 0) {
			for (int i = 0; i < num_of_elem / size; ++i) {
				cout << qbit_vec_local[i] << endl;
			}
		}
*/
	}
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	} else {
		MPI_Reduce(&compute_time, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Type_free(&double_double);
	if (mode == FILE) {
		MPI_Type_free(&subarr_type);
	}
	delete[] qbit_vec_local;
	delete[] result_local;
	if (rank == 0) {
		ofstream ofile("time_file", ios::app);
		ofile << n << "\t" << k << "\t" << size << "\t" << compute_time << endl;
	}
	MPI_Finalize();
	return 0;
}
