# Revised Simplex Method Solver
## Description
Written by Andres Byun. 

This is an implementation of a Linear Program solver using the Revised Simplex Method written in C++. 
This solution uses the library Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page) as a support for Linear Algebra operations.

This LP Solver has the following basic features:
 * Check for infeasibility.
 * Check for unboundedness.
 * Check for optimality, additionally it prints out the optimal value and the optimal assignment of the optimization variables.
 * Solve initially infeasible LPs with a Two-Phase Primal-Dual Method.
 * Cycle avoidance.

### Specific implementation
In order to solve LPs, this program initially checks for the initial primal-feasibility of an LP and then performs the Revised Simplex Method (the Linear Algebraic Simplex Method without any inverse matrix computations). This is done through a simple LU-Decomposition with partial pivoting and then solving for the vector (more on solving linear equations with LU-Decompositions: https://en.wikipedia.org/wiki/LU_decomposition#Solving_linear_equations). We check for unboundedness at this step.

If the LP is not primal-feasible, then we check for dual-feasibility, if the LP is dual-feasible we perform the Linear Algebraic Dual Simplex Method to obtain the optimal solution. This step is where we check for infeasibility as well.

If it is neither we run a two-phase primal-dual method which entails simply solving the dual with a 0-vector instead of the objective function and then using the basis returned from that function.

The solver also uses a combination of Largest Coefficient, combined with Bland's rule whenever a degenerate LP is found to avoid cycling. The check that a degenerate LP happens is simply encountering a 0 as the value of any of the variables in the basis. This is implemented in both primal and dual functions. Because Bland's Rule does not cycle, we can affirm this solver does not cycle either and must terminate in a finite number of steps. However, the performance of the solver might be affected if we encounter an LP with several degenerate variables as Bland's Rule is not performance-oriented.

Additionally, this solver uses the value of 1.0e-10 as a threshold for any numerical errors due to any faulty floating-point computation. This essentially means any absolute value that is smaller than 1.0e-10 will be treated as a 0.

And finally, to compute the optimal solution without an inverse, the solver uses the result of the dot product of the coefficients of the original objective function and the value of the optimal assigned variables.

## Compilation and Execution of the Program
To compile the program, simply run the `make` command on the terminal and the Makefile should take care of the rest.  
To run the program, simply type `./lp_solver` on the terminal and feed it an input LP. You can pipe the input as well. There are no command line arguments for this program and it only uses standard input and output through streams.

## Extra features
For marking, the specific extra features used in this program are:
1. Primal-Dual Methods (i.e. the two-phase simplex method described above)
2. Linear Algebraic Simplex Method (this implementation does not use the dictionary approach and uses Eigen as a support for linear algebraic operations)
3. Revised Simplex Method (use LU-Decompositions and solving them https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html#a49247bd2f742a46bca1f9c2bf1b19ad8 to avoid using inverse matrices)

## Citations
1. For tokenizing strings, Geeks for Geeks provided some useful code for inspiration https://www.geeksforgeeks.org/tokenizing-a-string-cpp/
2. For linear algebra operations, Eigen provided a robost library https://eigen.tuxfamily.org/index.php?title=Main_Page
