#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>

#include <omp.h>

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
/*
	for (int i = 0; i < num_of_elem; ++i) {
		cout << result[i] << endl;
	}
*/
}

int main(int argc, char* argv[])//n, k, t
{
	
	bool test = false;
	double prog_time1 = 0, prog_time2 = 0;
	int num_of_start = 1;
	int n = strtol(argv[1], NULL, 10),
		k = strtol(argv[2], NULL, 10),
		num_of_elem = 1 << n;
	int thrds = strtol(argv[3], NULL, 10);
	omp_set_num_threads(thrds);
	complexd* qbit_vec = new complexd[num_of_elem];
	complexd* result = new complexd[num_of_elem];
	double tmp = 1 / sqrt(2);
	double matrix[2][2] = {{tmp, tmp}, {tmp, - tmp}};
	for (int counter = 0; counter < num_of_start; ++counter) {
		double sum = 0;
		unsigned int thrd_seed;
		unsigned int cur_time = time(NULL);
		prog_time1 -= omp_get_wtime();
#pragma omp parallel  default(none) shared(qbit_vec, num_of_elem, sum, cur_time) private(thrd_seed)
		{
			thrd_seed = cur_time + omp_get_thread_num();
#pragma omp for reduction(+:sum)
			for (int i = 0; i < num_of_elem; i += 1) {
				qbit_vec[i] = complexd(rand_r(&thrd_seed), rand_r(&thrd_seed));
				//qbit_vec[i].real(rand_r(&thrd_seed));
				//qbit_vec[i].imag(rand_r(&thrd_seed));
				sum += norm(qbit_vec[i]);		
			}	
		}	
		sum = pow(sum, 0.5);
#pragma omp parallel for default(none) shared(qbit_vec, sum, num_of_elem)
		for (int i = 0; i < num_of_elem; i += 1) {
			qbit_vec[i] /= sum;
		}
		prog_time1 += omp_get_wtime();
		prog_time2 -= omp_get_wtime();
		quantum_transformation(qbit_vec, result, n, matrix, k);
		prog_time2 += omp_get_wtime();
	}
	delete[] qbit_vec;
	delete[] result;

	ofstream ofile("out", ios::app);
	ofile << n << "\t" << k << "\t"  << thrds << "\t" <<  prog_time1 / num_of_start << "@@@ " << prog_time2 / num_of_start << "  SUM " << (prog_time1 + prog_time2) / num_of_start << endl;	
	return 0;
}
