all: compile run

compile:
	g++ -O3 -o main.out main.cpp utils.cpp ising_solver.cpp -larmadillo

run:
	./main.out


allp: compilep runp

compilep:
	g++ -O3 -o main.out main.cpp utils.cpp ising_solver.cpp -larmadillo -fopenmp

runp:
	./main.out


all_testing: compile_testing run_testing

compile_testing:
	g++ -O3 -o testing.out testing.cpp utils.cpp ising_solver.cpp -larmadillo

run_testing:
	./testing.out


all_testing_basic: compile_testing_basic run_testing_basic # Runs testing.cpp without compiling the other .cpp files.

compile_testing_basic:
	g++ -o testing.out testing.cpp -larmadillo

run_testing_basic:
	./testing.out
