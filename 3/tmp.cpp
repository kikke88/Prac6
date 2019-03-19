#include <iostream>
#include <chrono>
#include <random>
#include <complex>

using namespace std;

typedef complex<double> complexd;

int main(int argc, char* argv[])
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, std::nextafter(1.0, std::numeric_limits<double>::max()));
	for (int i = 0; i < 10; ++i) {
		cout <<  dis(gen) << "  ";
	}
	cout << endl;
	return 0;
}
