#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>

#include <omp.h>

#define GEN 0
#define FILE 1
using namespace std;
typedef complex<double> complexd;
const double RAND_MAX_D = RAND_MAX;


void quantum_transformation(const complexd* a, complexd* result, const int n, const double u[][2], const int k)
{
	int num_of_elem = 1 << n;
	int test = 1 << (n - k);
#pragma omp parallel for default(none) shared(a, u, result, test, num_of_elem)
	for (int i = 0; i < num_of_elem; i += 1) {
		result[i] = u[i & test && 1][0] * a[i & ~(test)] + u[i & test && 1][1] * a[i | test];
	}
}

int main(int argc, char* argv[])//n, k, генерируешь сам 0 / файл 1,имя файла
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//double prog_time1 = 0, prog_time2 = 0;
	int n = strtol(argv[1], NULL, 10),
		k = strtol(argv[2], NULL, 10),
		num_of_elem = 1 << n;
        mode = stdtol(argv[3], NULL, 10);
	complexd* qbit_vec_cur = new complexd[num_of_elem / size];
    MPI_Datatype SUBARR_TYPE;
//	double tmp = 1 / sqrt(2);
//	double matrix[2][2] = {{tmp, tmp}, {tmp, -tmp}};
    if (mode == GEN) {
        unsigned int cur_time;
        if (rank == 0) {
            cur_time = time(NULL);
        }
        MPI_Bcast(&cur_time, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        cur_time += rank;
    	double sum_cur = 0;
    	for (int i = 0; i < num_of_elem / size; i += 1) {
    		qbit_vec_cur[i] = complexd(rand_r(&cur_time), rand_r(&cur_time));
    		sum_cur += norm(qbit_vec_cur[i]);
    	}
        MPI_Allreduce(MPI_IN_PLACE, &sum_cur, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sum_cur = pow(sum_cur, 0.5);
    	for (int i = 0; i < num_of_elem / size; i += 1) {
    		qbit_vec[i] /= sum_cur;
    	}
    } else {
        MPI_File in_file;

        int array_of_subsizes = num_of_elem / size;
        int start = array_of_subsizes * rank;
        MPI_Type_create_subarray(1, &num_of_elem, &array_of_subsizes, &start, MPI_ORDER_C, MPI_DOUBLE, &SUBARR_TYPE);
        MPI_Type_commit(&SUBARR_TYPE);
        
        MPI_File_open(MPI_COMM_WORLD, argv[4], MPI_MODE_RDONLY, MPI_INFO_NULL, &in_file);
        MPI_File_set_view(in_file, 0, MPI_DOUBLE, SUBARR_TYPE, "native", MPI_INFO_NULL);
        MPI_File_read_all(in_file, qbit_vec_cur, array_of_subsizes, MPI_DOUBLE, MPI_STATUS_IGNORE);

        MPI_File out_file;
        MPI_File_open(MPI_COMM_WORLD, "OUT", MPI_MODE_WDONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out_file);
        MPI_File_set_view(out_file, 0, MPI_DOUBLE, SUBARR_TYPE, "native", MPI_INFO_NULL);
		MPI_File_write_all(out_file, /*   результирующий вектор    */, array_of_subsizes, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}

        //

    }



		quantum_transformation(qbit_vec, result, n, matrix, k);

	delete[] qbit_vec;
	delete[] result;
	ofstream ofile("out", ios::app);
	ofile << n << "\t" << k << "\t"  << thrds << "\t" <<  prog_time1 / num_of_start << "@@@ " << prog_time2 / num_of_start << "  SUM " << (prog_time1 + prog_time2) / num_of_start << endl;
	return 0;
}
