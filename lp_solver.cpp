/* CSC 445 - PROGRAMMING PROJECT
 * Written by Andres Byun
 * ============================================
 * This is a program that perform the Revised
 * Simplex Method using LU-decompositions to
 * avoid computing inverse matrices.
 */

 // Include statements
#include <iostream>
#include <string>
#include <regex>
#include <vector>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

// Using declarations
using std::abs;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::regex;
using std::vector;
using std::stof;
using std::setprecision;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Method declarations
void tokenize_line(string str, vector<string>& tok);
void primal_simplex(MatrixXd& A, VectorXd b, VectorXd c, vector<int> B, vector<int> N);
bool dual_simplex(MatrixXd& A, VectorXd b, VectorXd c, vector<int>& B, vector<int>& N, bool two_phase);
void print_optimal(MatrixXd A_B, VectorXd b, VectorXd cB, VectorXd x, unsigned int size);
MatrixXd initialize_submatrix(MatrixXd A, vector<int> indexes);
VectorXd initialize_vec_zero (unsigned int size);

// Constant to control any potential numerical errors
const double EPSILON = 1.0e-10;

// Main function
int main() {
	// Read from standard input
	string input;
	vector<vector<string>> rows;
	vector<string> cols;

	while (std::getline(cin, input)) {
		if (input.empty()) continue; // Ignore empty lines
		tokenize_line(input, cols);
		rows.push_back(cols);
	}

	// Variables for the Linear Algebraic method
	MatrixXd A(rows.size() - 1, rows.size() + cols.size() - 2);
	VectorXd b(rows.size() - 1);
	VectorXd c(rows.size() + rows.at(0).size() - 1);
	vector<int> B_indexes, N_indexes;
	bool primal_feasible = false, dual_feasible = false;	// Flags to check if the two-phase Simplex Method is needed
	int pf_count = 0, df_count = 0;				// Counters to check the flags

	// Put input into equational form
	for (int i = 0; i < rows.size(); i++) {
		for (int j = 0; j < rows.at(i).size() + rows.size() - 1; j++) {
			// Fill the c vector
			if (i == 0) {
				if (j < rows.at(0).size()) {
					c[j] = stod(rows.at(0).at(j));
					N_indexes.push_back(j);
				}
				else c[j] = 0;
				if (c[j] <= 0) {
					df_count++;
					if (df_count == c.size()) dual_feasible = true;
				}
			}
			// Fill the b vector
			else if (j == rows.at(i).size() - 1) {
				b[i - 1] = stod(rows.at(i).at(rows.at(i).size() - 1));
				if (b[i - 1] >= 0) {
					pf_count++;
					if (pf_count == b.size()) primal_feasible = true;
				}
			}
			// Fill the A matrix
			else {
				if (j < rows.at(0).size())
					A(i - 1, j) = stod(rows.at(i).at(j));
				else {
					if (i - 1 == j - rows.at(0).size() - 1) {
						B_indexes.push_back(j - 1);
						A(i - 1, j - 1) = 1;
					}
					else A(i - 1, j - 1) = 0;
				}
			}
		}
	}

	// Perform the Simplex Algorithm depending of the circumstance
	if (primal_feasible) primal_simplex(A, b, c, B_indexes, N_indexes);
	else if (dual_feasible) dual_simplex(A, b, c, B_indexes, N_indexes, false);
	else {
		VectorXd zero = initialize_vec_zero(c.size());
		if (dual_simplex(A, b, zero, B_indexes, N_indexes, true))
			primal_simplex(A, b, c, B_indexes, N_indexes);
	}

	return 0;
}

// Function to initialize sub-matrices
MatrixXd initialize_submatrix(MatrixXd A, vector<int> indexes) {
	MatrixXd mat(A.rows(), indexes.size());

	for (int i = 0; i < indexes.size(); i++)
		for (int j = 0; j < A.cols(); j++)
			if (indexes.at(i) == j)
				for (int k = 0; k < A.rows(); k++)
					mat(k, i) = A(k, j);

	return mat;
}

// Function that returs a zero vector
VectorXd initialize_vec_zero (unsigned int size) {
	VectorXd vec(size);
	for (unsigned int i = 0; i < size; i++)
		vec[i] = 0;
	
	return vec;
}

// Function for the Revised Simplex Method
void primal_simplex(MatrixXd& A, VectorXd b, VectorXd c, vector<int> B, vector<int> N) {
	MatrixXd A_B = initialize_submatrix(A, B);
	MatrixXd A_N = initialize_submatrix(A, N);

	VectorXd x = initialize_vec_zero(c.size());
	VectorXd xB = A_B.partialPivLu().solve(b);	// Solve the system A_B(xB) = b with an LU-Decomposition

	// Check initial feasibility
	for (unsigned int i = 0; i < xB.size(); i++) {
		if (abs(xB[i]) < EPSILON) xB[i] = 0;
		if (xB[i] < 0) {
			cout << "ERROR: basis is not primal-feasible" << endl;
			return;
		}
	}
	
	while (true) {
		bool is_degenerate = false;
		VectorXd cB(B.size());
		VectorXd cN(N.size());

		for (unsigned int i = 0; i < B.size(); i++) {
			cB[i] = c(B[i]);
			if (xB[i] == 0) is_degenerate = true;	// Check for degeneracies
		}
		for (unsigned int i = 0; i < N.size(); i++) {
			cN[i] = c(N[i]);
		}

		VectorXd v = A_B.transpose().partialPivLu().solve(cB);
		VectorXd zN = (A_N.transpose() * v) - cN;

		// Check optimality
		unsigned int opt_count = 0;
		int largest = -1;
		int j = -1;
		for (unsigned int i = 0; i < zN.size(); i++) {
			if (abs(zN[i]) < EPSILON) zN[i] = 0;
			if (zN[i] >= 0) {
				opt_count++;
				if (opt_count == zN.size()) {
					print_optimal(A_B, b, cB, x, N.size());
					return;
				}
			}
			// Choose entering variable
			else {
				if (j == -1) j = i;
				else if (is_degenerate && (N[j] > N[i])) j = i;		// Bland's Rule
				else if (!is_degenerate && (-zN[j] < -zN[i])) j = i;	// Largest Coefficient Rule

			}
		}

		// Choose leaving variable
		VectorXd Aj = A.col(N[j]);
		VectorXd delta_xB = A_B.partialPivLu().solve(Aj);
		int unbound_count = 0;
		double t;
		int i = -1;
		for (unsigned int k = 0; k < B.size(); k++) {
			if (abs(delta_xB[k]) < EPSILON) delta_xB[k] = 0;
			if (delta_xB[k] > 0) {
				double increment = xB[k] / delta_xB[k];
				if(i == -1 || t > increment) {
					t = increment;
					i = k;
				}
			}
			// Check unboundedness
			else {
				unbound_count++;
				if (unbound_count == B.size()) {
					cout << "unbounded" << endl;
					return;
				}
			}
		}

		// Update for the next iteration
		xB = xB - (t * delta_xB);
		xB[i] = t;

		x[B[i]] = 0;
		x[N[j]] = t;

		int temp = N[j];
		N[j] = B[i];
		B[i] = temp;

		A_B = initialize_submatrix(A, B);
		A_N = initialize_submatrix(A, N);

		for (unsigned int k = 0; k < B.size(); k++) {
			x[B[k]] = xB[k];
		}
	}
	return;
}

// Function for the Revised Dual Simplex Method
// If the value of two_phase is true, the function does not print out optimal values so the Primal Simplex method can solve it
bool dual_simplex(MatrixXd& A, VectorXd b, VectorXd c, vector<int>& B, vector<int>& N, bool two_phase) {
	MatrixXd A_B = initialize_submatrix(A, B);
	MatrixXd A_N = initialize_submatrix(A, N);

	VectorXd cB(B.size());
	VectorXd cN(N.size());
	for (unsigned int i = 0; i < B.size(); i++) {
		cB[i] = c(B[i]);
	}
	for (unsigned int i = 0; i < N.size(); i++) {
		cN[i] = c(N[i]);
	}

	VectorXd v = A_B.transpose().partialPivLu().solve(cB);
	VectorXd zN = (A_N.transpose() * v) - cN;

	for (int i = 0; i < zN.size(); i++) {
		if (zN[i] < EPSILON) zN[i] = 0;
		if (zN[i] < 0) {
			cout << "ERROR: basis is not dual-feasible" << endl;
			return false;
		}
	}

	while (true) {
		VectorXd x = initialize_vec_zero(c.size());
		VectorXd xB = A_B.partialPivLu().solve(b);
		VectorXd xN = initialize_vec_zero(N.size());

		int i = -1;
		// Check for degeneracies
		bool is_degenerate = false;
		for (int k = 0; k < N.size(); k++) {
			if (abs(zN[k]) < EPSILON) zN[k] = EPSILON;
			if (zN[k] == 0) {
				is_degenerate = true;
			}
		}

		// Check for optimality
		unsigned int opt_count = 0;
		for (unsigned int k = 0; k < xB.size(); k++) {
			if (abs(xB[k]) < EPSILON) xB[k] = 0;
			if (xB[k] >= 0) {
				x[B[k]] = xB[k]; 
				opt_count++;
				if (opt_count == xB.size()) {	
					if (!two_phase) print_optimal(A_B, b, cB, x, N.size());
					return true;
				}
			}
			else {
				// Choose leaving variable
				if (i == -1) i = k;
				else if (is_degenerate && B[k] < B[i]) i = k;		// Bland's Rule
				else if (!is_degenerate && -xB[i] < -xB[k]) i = k;	// Largest Coefficient Rule
			}
		}


		// Choose entering variable
		VectorXd u = initialize_vec_zero(B.size());
		u[i] = 1;
		VectorXd v =  A_B.transpose().partialPivLu().solve(u);
		VectorXd delta_zN = -A_N.transpose() * v;
		int j = -1;
		double s;

		for (unsigned int k = 0; k < zN.size(); k++) {
			if (abs(delta_zN[k]) < EPSILON) delta_zN[k] = 0;
			if (delta_zN[k] > 0) {
				double increase = zN[k] / delta_zN[k];
				if (j == -1 || s > increase) {
					s = increase;
					j = k;
				}
			}
		}

		// Check for unboundedness for the dual/infeasibility for the primal
		if (j == -1) {
			cout << "infeasible" << endl;
			return false;
		}

		zN = zN - (s * delta_zN);
		zN[j] = s;

		// Update for next iteration
		unsigned int temp = B[i];
		B[i] = N[j];
		N[j] = temp;

		for (unsigned int k = 0; k < B.size(); k++) {
			cB[k] = c(B[k]);
		}
		for (unsigned int k = 0; k < N.size(); k++) {
			cN[k] = c(N[k]);
		}

		A_B = initialize_submatrix(A, B);
		A_N = initialize_submatrix(A, N);
	}
	return true;
}

// Simple function that prints out the optimal output of the LP
void print_optimal(MatrixXd A_B, VectorXd b, VectorXd cB, VectorXd x, unsigned int size) {
	double optimal = (A_B.partialPivLu().solve(b)).dot(cB);
	cout << "optimal" << endl;
	cout << setprecision(7) << optimal << endl;
	for (unsigned int i = 0; i < size; i++) {
		if (abs(x[i]) < EPSILON) x[i] = 0;
		string terminator = i != size - 1 ? " " : "\n";
		cout << x[i] << terminator;
	}

	return;
}

// Function that takes in a string, tokenizes the strings with white spaces as delimiters
// Returns a vector containing the tokens
// Code skeleton found in https://www.geeksforgeeks.org/tokenizing-a-string-cpp/
void tokenize_line(string str, vector<string>& tok) {
    const regex re(R"(\s+)");

    std::sregex_token_iterator it{ str.begin(), str.end(), re, -1 };
    vector<string> tokens{ it, {} };

    // Delete empty strings
    tokens.erase(
        std::remove_if(
            tokens.begin(),
            tokens.end(),
            [](string const& s) {
                return s.size() == 0;
            }),
        tokens.end());

	tok = tokens;
}
