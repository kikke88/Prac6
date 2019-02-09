#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>

#include <omp.h>

using namespace std;
typedef complex<double> complexd;
const double RAND_MAX_D = RAND_MAX;


void quantum_transformation(const complexd* a, const int n, const int num_of_elem, const double u[][2], const int k)
{
	for (int i = 0; i < num_of_elem; ++i) {
		int tmp1 = i & ~(1 << (n - k)),
			tmp2 = i | 1 << (n - k);
		cout << i << "     " << tmp1 << "   " << tmp2 << endl;
	}
	complexd* result = new complexd[num_of_elem];
#pragma omp parallel for shared(a, u, result)// schedule(dynamic)
	for (int i = 0; i < num_of_elem; i += 1) {
		int i_k = i >> (n - k) && 1;
		result[i] = u[i_k][0] * a[i & ~(1 << (n - k))] + u[i_k][1] * a[i | 1 << (n - k)];
	}
	for (int i = 0; i < num_of_elem; ++i) {
		cout << result[i] << endl;
	}
	delete[] result;
}



int main(int argc, char* argv[])//n k
{
	int n = stoi(argv[1]),
		k = stoi(argv[2]),
		num_of_elem = pow(2, n);
	double sum = 0;
	complexd* qbit_vec = new complexd[num_of_elem];
	unsigned cur_time = time(0);
#pragma omp parallel for default(none) shared(qbit_vec, cur_time, num_of_elem) reduction(+:sum)// schedule(dynamic)	
	for (int i = 0; i < num_of_elem; i += 1) {
		srand(cur_time + i);
		qbit_vec[i].real(i);
		qbit_vec[i].imag(0);
		//qbit_vec[i].real(rand() / RAND_MAX_D);
		//qbit_vec[i].imag(rand() / RAND_MAX_D);
		sum += abs(qbit_vec[i] * qbit_vec[i]);
	}
	sum = pow(sum, 0.5);
#pragma omp parallel for default(none) shared(sum, qbit_vec, num_of_elem)// schedule(dynamic)	
	for (int i = 0; i < num_of_elem; i += 1) {
		qbit_vec[i] /= sum;	
	}
	for (int i = 0; i < num_of_elem; ++i) {
		cout << qbit_vec[i] << endl;
	}
	double tmp = 1 / sqrt(2);
	double matrix[2][2] = {{tmp, tmp},
							{tmp, - tmp}};
	quantum_transformation(qbit_vec, n, num_of_elem, matrix, k);
	delete[] qbit_vec;
	return 0;
}