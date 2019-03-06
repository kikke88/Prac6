COMPILE = g++
.PHONY: all clean install uninstall
all: print gen main
clean:
	rm -rf gen print main tmp  *.out *.data
print: print.cpp
	$(COMPILE) -o print print.cpp
gen: generator.cpp
	$(COMPILE) -fopenmp -O3 -o gen generator.cpp
main: main.cpp
	mpicxx -O3 -o main main.cpp
