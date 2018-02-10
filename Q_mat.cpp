# include "./Headers/Generate_Matrix.h"
# include "./Headers/Output.h"

/* This function compares REF_LU vs QMAT */
int main (int argc, char* argv[])
{
	/* Declare variables and data structures */	
	int i, j, n, m, check, **A, **B;
	mpz_t **copy, **B_ref, **B_Q;							// Copies and exact versions of B
	REF_mat *A_ref, *Q_mat;								// Data structure for the REF Mat 
	A_ref = new REF_mat;								// Declare memory for A
	Q_mat = new REF_mat;

	/* Process the command line elements */	
	REF_LU_Options* option;								// The command line argument vector
	option = new REF_LU_Options;							// Allocate memory for option
	Options_set_defaults(option);							// Set default values
	Options_process_command_line(argc, argv, option);				// Process the command line

	/* Generation of matrix for use in the ensuing analysis */	
	n = option->dimension; m = option->dimension;					// Set the dimensions (assumed square)
	A = new int* [m];
	GM_initialize_matrix(A, n, m);							// Allocate memory for A
	GM_generate_matrix(A, n, m, option->lower, option->upper, option->seed1, option->seed2, option->density, true);	// Generate A

	/* Create any copies of input matrix */	
	REF_alloc(A_ref, n, m, option->numRHS);						// Allocate memory for REF LU A
	REF_alloc(Q_mat, n, m, option->numRHS);						// Allocate memory for Q_mat
	REF_make_copy (A, A_ref->mat, n, m);						// Create a copy of A into integer arithmetic
	REF_make_copy(A, Q_mat->mat, n, m);
		
	/* Factorization */	
	check = REF_LU_Factorization_pivots(A_ref);					// REF LU 
	if (check == 0) return 0;							// Check for errors

	mpz_t det; mpz_init(det);
	Q_Mat_get_adjunct(Q_mat, det);							// Q matrix

	/* Solve linear system */		
	if (option->numRHS > 0)								// Number of RHS vectors
	{
		B = new int* [n];
		GM_initialize_matrix(B, option->numRHS, n);				// Allocate B
		GM_generate_matrix(B, option->numRHS, n, option->lower, option->upper, option->seed1, option->seed2, option->density, false);	// Generate B
		B_ref = new mpz_t * [n];						// Allocate memory for REF RHS vectors
		for (i = 0; i < n; i++)
			B_ref[i] = new mpz_t [option->numRHS];
		REF_initialize_matrix(option->numRHS,n,B_ref);
		REF_solve(A_ref, B_ref, B);						// Solve LD^-1U x= b
		steady_clock::time_point t2_s_begin = steady_clock::now();
		B_Q = new mpz_t * [n];							// Allocate memory for Qmat RHS vectors
		for (i = 0; i < n; i++)
			B_Q[i] = new mpz_t [option->numRHS];
		REF_initialize_matrix(option->numRHS,n,B_Q);
		REF_make_copy (B, B_Q, option->numRHS, n);				// Create a copy of B into integer arithmetic
		Q_Mat_solve(Q_mat, B_Q, option->numRHS);				// Solve x = Adj(A) * b
		Q_Mat_get_x(Q_mat, option->numRHS, det);				// xq = x / det
		steady_clock::time_point t2_s_end = steady_clock::now();
		Q_mat->t_s = duration_cast<duration<float>>(t2_s_end - t2_s_begin);	// Compute solution time
	}		
	mpz_clear(det);
	
	/* Check Solutions if desired */
	if (option->check == 1 && option->numRHS > 0)
	{
		copy = new mpz_t * [m];							// Allocate memory for copy
		for (i = 0; i < n; i++)
			copy[i] = new mpz_t [n];
		REF_initialize_matrix(n, m, copy);
		REF_make_copy(A, copy, n, m);						// Copy A into int arithmetic
		A_ref->check = REF_check_x(A_ref,B_ref,copy);				// Check A * REF_x = b
		Q_mat->check = REF_check_x(Q_mat, B_ref, copy);				// Check A * Q_x = b
		REF_clear_mat(n, A_ref->m, copy);					// Clear memory of copy
	}
	
	/* Output */
	print_stats_QMat(A_ref, Q_mat, option);						// Timing statistics
	QMat_print_to_file(option, A_ref, Q_mat);					// Print results to file
	clear_QMat_memory(option, A, B, A_ref, Q_mat, B_ref, B_Q);			// Free memory	
}