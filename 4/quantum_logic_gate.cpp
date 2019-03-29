
#include <complex>

using namespace std;
typedef complex<double> complexd;

class Quantum_vector {
private:
    const int num_of_qubits;
    const int num_of_elem;
    const int size;//
    const int rank;
    const int elem_in_one_proc;
    MPI_Datatype double_double;//*
    complexd* vector_1;//*
    complexd* vector_2;//*
    void file_manipulation(const string, int);
public:
    Quantum_vector(int, int, int);
    ~Quantum_vector();
    void vector_filling();
    void vector_filling(string);

    void Hadamard_gate(int);


};

Quantum_vector::Quantum_vector(int _num_of_qubits, int rank = 0, int _size = 1):
    num_of_qubits { _num_of_qubits },
    num_of_elem { 2 << num_of_qubits },
    size { _size },// мб потом можно убрать
    elem_in_one_proc{ num_of_elem / size },
    vector_1{ new complexd[elem_in_one_proc] },
    vector_2{ new complexd[elem_in_one_proc] },
{
    MPI_Type_contiguous(2, MPI_DOUBLE, &double_double);
    MPI_Type_commit(&double_double);
}

Quantum_vector::~Quantum_vector()
{
    delete[] vector_1;
    delete[] vector_2;
    MPI_Type_free(&double_double);
}

void Quantum_vector::vector_filling()
{
    unsigned int cur_time;
    if (rank == 0) {
        cur_time = time(NULL);
    }
    if (size != 1) {
        MPI_Bcast(&cur_time, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }
    double sum = 0;
    const int tmp_const = RAND_MAX / 2;
#pragma omp parallel
    {
        unsigned int seed = cur_time + rank * omp_get_num_threads() + omp_get_thread_num();//
        //попробовать, можно ли определять количество используемых в программе нитей в конструкторе класса/////
#pragma omp for reduction(+:sum)
        for (int i = 0; i < elem_in_one_proc; ++i) {
            vector_1[i] = complexd(rand_r(&seed) - tmp_const, rand_r(&seed) - tmp_const);
            sum += norm(vector_1[i]);
        }
    }
    if (size != 1) {
        MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    sum = sqrt(sum);
#pragma omp for
    for (int i = 0; i < elem_in_one_proc; ++i) {
        vector_1[i] /= sum;
    }
}

void Quantum_vector::file_manipulation(const string )
void Quantum_vector::vector_filling(const string input_file_name)
{
    int input_flag = 1;
    file_manipulation(input_file_name, input_flag);
}

void Quantum_vector::Hadamard_gate(int k)
{

}
