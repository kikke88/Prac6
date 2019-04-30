COMPILE = g++
.PHONY: all clean install
all: main print gen
clean:
	rm -rf gen print main tmp not cnot compare *.out *.data time_file
print: print.cpp
	${COMPILE} -o print print.cpp
lint:
	cpplint --filter=-legal/copyright main.cpp
	cpplint --filter=-legal/copyright --extensions=hpp quantum_logic_gate.hpp
gen: generator.cpp
	${COMPILE} -fopenmp -O3 -o gen generator.cpp
main: main.cpp quantum_logic_gate.hpp
	mpicxx -Wall -Werror -fopenmp -o main main.cpp
