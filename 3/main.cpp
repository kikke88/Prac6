#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random>

#include "mpi.h"
#include "omp.h"

using namespace std;
typedef complex<double> complexd;

void noisy_matrix(double* cur_matrix, const double Xi, const double H[4], double EPS)
{
	double angle = Xi * EPS;
	double U_matrix[4] = {cos(angle), sin(angle), -sin(angle), cos(angle)};
	cur_matrix[0] = H[0] * U_matrix[0] + H[1] * U_matrix[2];
	cur_matrix[1] = H[0] * U_matrix[1] + H[1] * U_matrix[3];
	cur_matrix[2] = H[2] * U_matrix[0] + H[3] * U_matrix[2];
	cur_matrix[3] = H[2] * U_matrix[1] + H[3] * U_matrix[3];
	return;
}

void quantum_transformation(complexd* a, complexd* result, const int n, const double u[4],
							const int k, const int size, const int rank, const MPI_Datatype double_double)
{
	int num_of_elem = 1 << n;
	int test = 1 << (n - k);
	int log_k = log2(size);
	if (k > log_k) {
#pragma omp parallel for
		for (int i = 0; i < num_of_elem / size; ++i) {
			result[i] = u[(i & test && 1) * 2] * a[i & ~(test)] + u[(i & test && 1 ) * 2 + 1] * a[i | test];
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
#pragma omp parallel for
		for (int i = 0; i < num_of_elem / size; ++i) {
			result[i] = u[swap_bit * 2] * a[i] + u[swap_bit * 2 + 1] * recv_vec[i];
		}
		if (swap_bit == 1) {
			recv_vec = a;
		}
		delete[] recv_vec;
	}
	return;
}

void file_manipulation(complexd* vec, const int num_of_elem, const int size, const int rank,
						const MPI_Datatype double_double, const char* filename, bool in_out_flag)
{
	MPI_Datatype subarr_type;
	MPI_File file;
	int array_of_subsizes = num_of_elem / size;
	int start = array_of_subsizes * rank;
	MPI_Type_create_subarray(1, &num_of_elem, &array_of_subsizes, &start, MPI_ORDER_C, double_double, &subarr_type);
	MPI_Type_commit(&subarr_type);
	if (in_out_flag) {
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
		MPI_File_set_view(file, 0, double_double, subarr_type, "native", MPI_INFO_NULL);
		MPI_File_read_all(file, vec, array_of_subsizes, double_double, MPI_STATUS_IGNORE);
	} else {
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
		MPI_File_set_view(file, 0, double_double, subarr_type, "native", MPI_INFO_NULL);
		MPI_File_write_all(file, vec, array_of_subsizes, double_double, MPI_STATUS_IGNORE);
	}
	MPI_Type_free(&subarr_type);
}

void COMPLEXD_SUM(void *in, void *inout, int *len, MPI_Datatype *type)
{
	for (int i = 0; i < *len; ++i) {
		((complexd*)inout)[i] = ((complexd*)in)[i] + ((complexd*)inout)[i];
	}
}

#define GEN 0
#define FILE 1

int main(int argc, char* argv[])//n, EPS, (1, filename, out_filename) / (0, index)
{
	omp_set_num_threads(4);
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = strtol(argv[1], NULL, 10),
		num_of_elem = 1 << n;
	double EPS = strtod(argv[2], NULL),
		compute_time = 0.0;
	complexd* qbit_vec = new complexd[num_of_elem / size];
	complexd* noisy_result_1 = new complexd[num_of_elem / size];
	complexd* noisy_result_2 = new complexd[num_of_elem / size];
	complexd* normal_result_1 = new complexd[num_of_elem / size];
	complexd* normal_result_2 = new complexd[num_of_elem / size];
	MPI_Datatype double_double;
	MPI_Type_contiguous(2, MPI_DOUBLE, &double_double);
	MPI_Type_commit(&double_double);
	bool input_type;
	int num_of_starts;
	int FLAG = strtol(argv[3], NULL, 10);
	if (FLAG) {
		input_type = FILE;
		num_of_starts = 1;
	} else {
		input_type = GEN;
		num_of_starts = 120;///////
	}
	double tmp = 1 / sqrt(2);
	double Hadamard[4] = {tmp, tmp, tmp, -tmp};
	MPI_Op complexd_sum;
	MPI_Op_create(&COMPLEXD_SUM, 1, &complexd_sum);
	for (int k = 0; k < num_of_starts; ++k) {
		if (input_type == GEN) {
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
				unsigned int seed = cur_time + rank * omp_get_num_threads() + omp_get_thread_num();
#pragma omp for reduction(+:sum)
				for (int i = 0; i < num_of_elem / size; ++i) {
					qbit_vec[i] = complexd(rand_r(&seed) - tmp_const, rand_r(&seed) - tmp_const);
					sum += norm(qbit_vec[i]);
				}
			}
			if (size != 1) {
				MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			}
			sum = sqrt(sum);
#pragma omp for
			for (int i = 0; i < num_of_elem / size; ++i) {
				qbit_vec[i] /= sum;
			}
		} else {
			file_manipulation(qbit_vec, num_of_elem, size, rank, double_double, argv[4], 1);
		}
		double* noise = new double[n];
		if (rank == 0) {
			random_device rd;
			mt19937 gen(rd());
			normal_distribution<> dis(0,1);
//			uniform_real_distribution<> dis(0.0, nextafter(1.0, numeric_limits<double>::max()));
			for (int i = 0; i < n; ++i) {
				noise[i] = dis(gen);
			}
		}
		if (size != 1) {
			MPI_Bcast(noise, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		compute_time -= MPI_Wtime();
		//j = 1
		double* cur_matrix = new double[4];
		noisy_matrix(cur_matrix, noise[0], Hadamard, EPS);
		quantum_transformation(qbit_vec, noisy_result_1, n, cur_matrix, 1, size, rank, double_double);
		quantum_transformation(qbit_vec, normal_result_1, n, Hadamard, 1, size, rank, double_double);
		for (int j = 2; j < n + 1; ++j) {
			noisy_matrix(cur_matrix, noise[j - 1], Hadamard, EPS);
			quantum_transformation(noisy_result_1, noisy_result_2, n, cur_matrix, j, size, rank, double_double);
			swap(noisy_result_1, noisy_result_2);
			//quantum_transformation(normal_result_1, normal_result_2, n, Hadamard, j, size, rank, double_double);
			//swap(normal_result_1, normal_result_2);

		}
		compute_time += MPI_Wtime();

		delete[] cur_matrix;
		complexd Fidelity (0.0, 0.0);
		complexd Res_fidelity (0.0, 0.0);
		for (int j = 0; j < num_of_elem / size; ++j) {
			Fidelity += conj(noisy_result_1[j]) * normal_result_1[j];
		}
		if (size != 1) {
			MPI_Reduce(&Fidelity, &Res_fidelity, 1, double_double, complexd_sum, 0, MPI_COMM_WORLD);
		}
		if (rank == 0) {
			double res = norm(Res_fidelity);
			ofstream ofile("F_" + to_string(n) + "_" + argv[2], ios::app);
			ofile << 1 - res << endl;
			//cout << 1 - res << endl;
		}
		if (FLAG && argc == 6) {
			cout << "~~~~~~~~" << endl;
			file_manipulation(noisy_result_1, num_of_elem, size, rank, double_double, argv[5], 0);
		}
	}
/*
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	} else {
		MPI_Reduce(&compute_time, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}
	if (rank == 0) {
		ofstream ofile("time_file" + string(argv[4]), ios::app);
		ofile << argv[4] << "  " << compute_time / num_of_starts << endl;
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Type_free(&double_double);
	MPI_Op_free(&complexd_sum);
	delete[] qbit_vec;
	delete[] noisy_result_1;
	delete[] noisy_result_2;
	delete[] normal_result_1;
	delete[] normal_result_2;
	MPI_Finalize();
	return 0;
}
