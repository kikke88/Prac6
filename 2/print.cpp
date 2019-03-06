#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;
typedef complex<double> complexd;

int main(int argc, char* argv[])//n
{
	int n = strtol(argv[1], NULL, 10),
		num_of_elem = 1 << n;
	complexd* qbit_vec = new complexd[num_of_elem];
	ifstream in_file(argv[2], ios::binary);
	in_file.read((char *) qbit_vec, num_of_elem * sizeof(complexd));
	cout << n << endl;
	for (int i = 0; i < num_of_elem; ++i) {
		cout << qbit_vec[i] << endl;
	}
    delete[] qbit_vec;
    return 0;
}
