#include <iostream>
#include <complex>
#include <fstream>

using namespace std;
typedef complex<double> complexd;

int main(int argc, char* argv[])//n k infile outfile
{
	int n = strtol(argv[1], NULL, 10),
		num_of_elem = 1 << n;
    int k = strtol(argv[2], NULL, 10);
	complexd* qbit_vec_orig = new complexd[num_of_elem];
	ifstream in_file(argv[3], ios::binary);
	in_file.read((char *) qbit_vec_orig, num_of_elem * sizeof(complexd));
    complexd* qbit_vec_res = new complexd[num_of_elem];
	ofstream out_file(argv[4], ios::binary);
    int test = 1 << (n - k);
	for (int i = 0; i < num_of_elem; ++i) {
        qbit_vec_res[i] = qbit_vec_orig[i ^ test];
	}
    out_file.write((char *) qbit_vec_res, num_of_elem * sizeof(complexd));
    delete[] qbit_vec_orig;
    delete[] qbit_vec_res;
    return 0;
}
