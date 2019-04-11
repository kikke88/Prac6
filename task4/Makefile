COMPILE = g++
.PHONY: all clean install uninstall
all: main print gen
clean:
	rm -rf gen print main tmp  *.out *.data time_file
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
	make clean && make
