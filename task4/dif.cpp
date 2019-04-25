#include <iostream>
#include <complex>
#include <fstream>
#include <float.h>


using namespace std;
typedef complex<double> complexd;

int main(int argc, char* argv[])//n file1 file2
{
	int n = strtol(argv[1], NULL, 10),
		num_of_elem = 1 << n;
	complexd* qbit_vec1 = new complexd[num_of_elem];
    complexd* qbit_vec2 = new complexd[num_of_elem];
	ifstream file1(argv[2], ios::binary);
	file1.read((char *) qbit_vec1, num_of_elem * sizeof(complexd));
    ifstream file2(argv[3], ios::binary);
	file2.read((char *) qbit_vec2, num_of_elem * sizeof(complexd));
	for (int i = 0; i < num_of_elem; ++i) {

		if (abs(qbit_vec1[i].real() - qbit_vec2[i].real()) > 10e3 * DBL_EPSILON or
            abs(qbit_vec1[i].imag() - qbit_vec2[i].imag()) > 10e3 * DBL_EPSILON) {
            std::cout << DBL_EPSILON << "  differ  " << i << std::endl;
            delete[] qbit_vec1;
            delete[] qbit_vec2;
            return 1;
        }
	}
    delete[] qbit_vec1;
    delete[] qbit_vec2;
    return 0;
}
