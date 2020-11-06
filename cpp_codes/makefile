all: compile run

compile:
	g++ -O3 -o main.out main.cpp utils.cpp ising_solver.cpp -larmadillo

run:
	./main.out



all_testing: compile_test run_test

compile_testing:
	g++ -O3 -o testing.out testing.cpp utils.cpp ising_solver.cpp -larmadillo

run_testing:
	./testing.out
