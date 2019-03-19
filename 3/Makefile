COMPILE = g++
.PHONY: all clean install uninstall
all: main print gen
clean:
	rm -rf gen print main tmp  *.out *.data time_file
print: print.cpp
	${COMPILE} -o print print.cpp
gen: generator.cpp
	${COMPILE} -fopenmp -O3 -o gen generator.cpp
main: main.cpp
	mpicxx -fopenmp -o main main.cpp
