#include <iostream>
#include <complex>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

#include <omp.h>

using namespace std;
typedef complex<double> complexd;

int main(int argc, char* argv[])//n, ofile_name
{

	int n = strtol(argv[1], NULL, 10),
		num_of_elem = 1 << n;
	int thrds = 2;
	omp_set_num_threads(thrds);
	complexd* qbit_vec = new complexd[num_of_elem];
	double sum = 0;
	const int tmp_const = RAND_MAX / 2;
	unsigned int thrd_seed;
	unsigned int cur_time = time(NULL);
#pragma omp parallel  shared(qbit_vec, num_of_elem, sum, cur_time) private(thrd_seed)
	{
		thrd_seed = cur_time + omp_get_thread_num();
#pragma omp for reduction(+:sum)
		for (int i = 0; i < num_of_elem; i += 1) {
			//qbit_vec[i] = complexd(i, i * i);
			qbit_vec[i] = complexd(rand_r(&thrd_seed) - tmp_const, rand_r(&thrd_seed) - tmp_const);
			sum += norm(qbit_vec[i]);
		}
    }
	sum = pow(sum, 0.5);
#pragma omp parallel for default(none) shared(qbit_vec, sum, num_of_elem)
	for (int i = 0; i < num_of_elem; i += 1) {
		qbit_vec[i] /= sum;
	}

    ofstream ofile(argv[2], ios::binary);
    ofile.write((char *) qbit_vec, num_of_elem * sizeof(complexd));
    delete[] qbit_vec;
    return 0;
}
