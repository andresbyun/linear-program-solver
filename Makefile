all: lp_solver

debug: lp_solver.cpp
	g++ -I . -g lp_solver.cpp -o lp_debug

lp_solver: lp_solver.cpp
	g++ -O3 -std=c++17 -I . lp_solver.cpp -o lp_solver

clean:
	rm -rf *.o lp_solver
