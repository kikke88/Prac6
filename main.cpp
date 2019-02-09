#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>

#include <omp.h>

using namespace std;
typedef complex<double> complexd;
const double RAND_MAX_D = RAND_MAX;

int main(int argc, char* argv[])//n k
{
	int n = stoi(argv[1]),
		k = stoi(argv[2]),
		num_of_elem = pow(2, n);
	double sum = 0;
	complexd* qbit_vec = new complexd[num_of_elem];
	unsigned cur_time = time(0);
#pragma omp parallel for shared(qbit_vec, cur_time) reduction(+:sum) schedule(dynamic)	
	for (int i = 0; i < num_of_elem; i += 1) {
		srand(cur_time + i);
		qbit_vec[i].real(rand() / RAND_MAX_D);
		qbit_vec[i].imag(rand() / RAND_MAX_D);
		sum += abs(qbit_vec[i] * qbit_vec[i]);
//		cout << qbit_vec[i] << endl;
	}
	sum = pow(sum, 0.5);
//	cout << "_____" << sum << "_____" << endl;
#pragma omp parallel for shared(sum, qbit_vec) schedule(dynamic)	
	for (int i = 0; i < num_of_elem; i += 1) {
		qbit_vec[i] /= sum;	
//		std::cout << qbit_vec[i] << std::endl;
	}
	delete[] qbit_vec;
	return 0;
}