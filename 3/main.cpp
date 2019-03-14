#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random>

#include "mpi.h"
#include "omp.h"

using namespace std;
//using namespace std::chrono;
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

int main(int argc, char* argv[])//n, EPS
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = strtol(argv[1], NULL, 10),
		num_of_elem = 1 << n;
	double EPS = strtod (argv[2], NULL);
	double compute_time;
	complexd* qbit_vec = new complexd[num_of_elem / size];
	MPI_Datatype double_double;
	MPI_Type_contiguous(2, MPI_DOUBLE, &double_double);
	MPI_Type_commit(&double_double);
	//int array_of_subsizes = num_of_elem / size;
    unsigned int cur_time;
    if (rank == 0) {
        cur_time = time(NULL);
    }
	if (size != 1) {
		MPI_Bcast(&cur_time, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	}
    cur_time += rank;
	double sum = 0;
	const int tmp_const = RAND_MAX / 2;
	for (int i = 0; i < num_of_elem / size; ++i) {
		qbit_vec[i] = complexd(rand_r(&cur_time) - tmp_const, rand_r(&cur_time) - tmp_const);
		sum += norm(qbit_vec[i]);
	}
	if (size != 1) {
		MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
    sum = sqrt(sum);
	for (int i = 0; i < num_of_elem / size; ++i) {
		qbit_vec[i] /= sum;
	}

	complexd* noisy_result_1 = new complexd[num_of_elem / size];
	complexd* noisy_result_2 = new complexd[num_of_elem / size];
	complexd* normal_result_1 = new complexd[num_of_elem / size];
	complexd* normal_result_2 = new complexd[num_of_elem / size];
	double tmp = 1 / sqrt(2);
	double Hadamard[4] = {tmp, tmp, tmp, -tmp};

	double* noise = new double[n];;
	if (rank == 0) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0.0, std::nextafter(1.0, std::numeric_limits<double>::max()));
		for (int i = 0; i < n; ++i) {
			noise[i] = dis(gen);
		}
	}
	if (size != 1) {
		MPI_Bcast(noise, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	//j = 1
	double* cur_matrix = new double[4];
	noisy_matrix(cur_matrix, noise[0], Hadamard, EPS);
	quantum_transformation(qbit_vec, noisy_result_1, n, cur_matrix, 1, size, rank, double_double);
	quantum_transformation(qbit_vec, normal_result_1, n, Hadamard, 1, size, rank, double_double);
	for (int j = 2; j < n + 1; ++j) {
		noisy_matrix(cur_matrix, noise[j - 1], Hadamard, EPS);//вермя надо засекать только для испорченной матрицы

		quantum_transformation(noisy_result_1, noisy_result_2, n, cur_matrix, j, size, rank, double_double);
		std::swap(noisy_result_1, noisy_result_2);
		quantum_transformation(normal_result_1, normal_result_2, n, Hadamard, j, size, rank, double_double);//
		std::swap(normal_result_1, normal_result_2);//
	}
	double Fidelity = 0.0;
	for (int j = 0; j < num_of_elem / size; ++j) {
		Fidelity += abs(conj(noisy_result_1[j]) * normal_result_1[j]);
	}
	if (size != 1) {
		if (rank == 0) {
			MPI_Reduce(MPI_IN_PLACE, &Fidelity, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce(&Fidelity, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}
	}

	if (rank == 0) {
		cout << "!!" << Fidelity << "!!" << endl;
		cout << 1 - Fidelity << endl;
	}
	//compute_time -= MPI_Wtime();
	//compute_time += MPI_Wtime();
/*
	if (mode == FILE) {
		MPI_File out_file;
        MPI_File_open(MPI_COMM_WORLD, argv[5], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out_file);
        MPI_File_set_view(out_file, 0, double_double, subarr_type, "native", MPI_INFO_NULL);
		MPI_File_write_all(out_file, result, array_of_subsizes, double_double, MPI_STATUS_IGNORE);

	} else {
//вывод, случай параллельной генерации
		if (rank == 0) {
			for (int i = 0; i < num_of_elem / size; ++i) {
				cout << qbit_vec[i] << endl;
			}
		}
	}
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	} else {
		MPI_Reduce(&compute_time, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Type_free(&double_double);

	delete[] cur_matrix;
	delete[] qbit_vec;
	delete[] noisy_result_1;
	delete[] noisy_result_2;
	delete[] normal_result_1;
	delete[] normal_result_2;
/*
	if (rank == 0) {
		ofstream ofile("time_file", ios::app);
		ofile << n << "\t" << k << "\t" << size << "\t" << compute_time << endl;
	}
*/
	MPI_Finalize();
	return 0;
}
