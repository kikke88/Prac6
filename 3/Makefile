COMPILE = g++
.PHONY: all clean install uninstall
all: print gen main tmp
clean:
	rm -rf gen print main tmp  *.out *.data time_file
print: print.cpp
	$(COMPILE) -o print print.cpp
gen: generator.cpp
	$(COMPILE) -fopenmp -O3 -o gen generator.cpp
main: main.cpp
	mpicxx -o main main.cpp
tmp: tmp.cpp
	${COMPILE} -o tmp tmp.cpp
