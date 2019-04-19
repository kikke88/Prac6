COMPILE = g++
.PHONY: all clean install
all: main print gen differ not cnot
clean:
	rm -rf gen print main tmp not cnot compare *.out *.data time_file
print: print.cpp
	${COMPILE} -o print print.cpp
lint:
	cpplint --filter=-legal/copyright main.cpp
	cpplint --filter=-legal/copyright quantum_logic_gate.hpp
gen: generator.cpp
	${COMPILE} -fopenmp -O3 -o gen generator.cpp
main: main.cpp quantum_logic_gate.hpp
	mpicxx -Wall -Werror -fopenmp -o main main.cpp
test:
	./n_Had_test.sh
	./not_test.sh
	./cnot_test.sh
differ: dif.cpp
	${COMPILE} -o compare dif.cpp
not: not_check.cpp
	${COMPILE} -o not not_check.cpp
cnot:
	${COMPILE} -o cnot cnot_check.cpp
