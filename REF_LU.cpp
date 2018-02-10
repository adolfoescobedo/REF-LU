# include "./Headers/Generate_Matrix.h"
# include "./Headers/Output.h"

/* This programs performs an LU factorization of a randomly generated matrix. Please see README.txt for explanation */
int main (int argc, char* argv[])
{
	/* Declare variables and data structures */	
	int i, j, m, n, check, **A, **B;
	REF_mat* A_ref;								// Data structure for the REF Mat 
	QLU_mat *A_QC, *A_QD;							// Data structure for QLU Crout/Doolittle
	mpz_t **copy, **B_ref;							// Copies and exact versions of B
	mpq_t **copy_q, **B_QC, **B_QD;
	
	A_ref = new REF_mat;							// Declare memory for A
	A_QC = new QLU_mat;
	A_QD = new QLU_mat;

	/* Process the command line elements */	
	REF_LU_Options* option;							// The command line argument vector
	option = new REF_LU_Options;						// Allocate memory for option
	Options_set_defaults(option);						// Set default values
	check = Options_process_command_line(argc, argv, option);		// Process the command line
	if (check == 0) return 0;						// Check for errors
	
	/* Generation of matrix for use in the ensuing analysis */	
	n = option->dimension; m = option->dimension;				// Set the dimensions (assumed square)
	A = new int* [m];
	GM_initialize_matrix(A, n, m);						// Allocate memory for A
	GM_generate_matrix(A, n, m, option->lower, option->upper, option->seed1, option->seed2, option->density, true);	// Generate A
	if (option->numRHS > 0)							// Generate RHS vectors if necessary
	{
		B = new int* [n];
		GM_initialize_matrix(B, option->numRHS, n);			// Allocate B
		GM_generate_matrix(B, option->numRHS, n, option->lower, option->upper, option->seed1, option->seed2, option->density, false);	// Generate B
	}
	
	/* Perform factorizations */	
	if (option->REF == 1)							// Perform REF_LU Factorization
	{
		REF_alloc(A_ref, n, m, option->numRHS);				// Allocate memory for REF LU A
		REF_make_copy (A, A_ref->mat, n, m);				// Create a copy of A into integer arithmetic
		check = REF_LU_Factorization_pivots(A_ref);			// REF LU Factorization
		if (check == 0) return 0;					// Check for errors

		/* Forward and backward Sub */		
		if (option->numRHS > 0)						// Number of RHS vectors
		{
			B_ref = new mpz_t * [n];				// Allocate memory for REF RHS vectors
			for (i = 0; i < n; i++)
				B_ref[i] = new mpz_t [option->numRHS];
			REF_initialize_matrix(option->numRHS, n, B_ref);
			REF_solve(A_ref, B_ref, B);
		}		
	}
	if (option->Crout == 1)							// Crout LU
	{
		QLU_alloc(A_QC, n, m, option->numRHS);				// Allocate memory for crout A
		QLU_make_copy(A, A_QC->mat, n, m);				// Create a copy of A into rational arithmetic
		check = QLU_crout_LU_pivots(A_QC);				// Crout LU Factorization
		if (check == 0) return 0;					// Check for errors
		
		/* Forward and backward sub */
		if (option->numRHS > 0)						// Number of RHS vectors
		{
			B_QC = new mpq_t * [m];					// Allocate memory for crout RHS vectors
			for (i = 0; i < m; i++)
				B_QC[i] = new mpq_t [option->numRHS];
			QLU_initialize_matrix(A_QC->numRHS, n, B_QC);
			QLU_solve(A_QC, B_QC, B, true);				// Solve LUx = b for Crout
		}
	}
	if (option->Doolittle == 1)
	{
		QLU_alloc(A_QD, n, m, option->numRHS);				// Allocate memory for Doolittle A
		QLU_make_copy(A, A_QD->mat, n, m);				// Copy A into rational arithmetic
		check = QLU_doolittle_LU_pivots(A_QD);				// Doolittle LU Factorization
		if (check == 0) return 0;					// Check for errors

		/* Forward and backward sub */
		if (option->numRHS > 0)						// Number of RHS vectors
		{
			B_QD = new mpq_t *[m];					// Allocate memory for Doolittle RHS vectors
			for (i = 0; i < m; i++)
				B_QD[i] = new mpq_t [option->numRHS];
			QLU_initialize_matrix(option->numRHS, m, B_QD);
			QLU_solve(A_QD, B_QD, B, false);			// Solve LUx = b for Doolittle
		}
	}
	
	/* Check Solutions if desired */
	if (option->check == 1 && option->numRHS > 0)
	{
		if (option->REF == 1) 						// Check REF LU
		{
			copy = new mpz_t *[m];					// Allocate memory for copy
			for (i = 0; i < n; i++)
				copy[i] = new mpz_t [n];
			REF_initialize_matrix(n, m, copy);
			REF_make_copy(A, copy, n, m);				// Create copy
			A_ref->check = REF_check_x(A_ref, B_ref, copy);
			REF_clear_mat(n, m, copy);
		}
		if (option->Doolittle == 1 || option->Crout == 1)		// Only 1 rational copy is necesary for solution verification
		{
			copy_q = new mpq_t *[m];				// Allocate memory for copy_q
			for (i = 0; i < n; i++)
				copy_q[i] = new mpq_t [n];
			QLU_initialize_matrix(n, m, copy_q);
			QLU_make_copy(A, copy_q, n, m);				// Copy the matrix
		}
		if (option->Crout == 1) 					// Check Crout
		{
			A_QC->check = QLU_check_x(A_QC, B_QC, copy_q);
			if (option->Doolittle != 1) QLU_clear_mat(n, A_ref->m, copy_q);
		}
		if (option->Doolittle == 1)					// Check Doolittle
		{
			A_QD->check = QLU_check_x(A_QD, B_QD, copy_q);
			QLU_clear_mat(n, A_ref->m, copy_q);
		}
	}
	
	/* Output */	
	if (option->print2 == 1)						// Print x vector(s) to a file if desired
		if (option->numRHS > 0)
		{
			if (option->REF == 1)
				print_x_vector(option, A_ref->xq, A_ref->n, A_ref->numRHS);
			else if (option->Crout == 1)
				print_x_vector(option, A_QC->x, A_QC->n, A_QC->numRHS);
			else if (option->Doolittle == 1)
				print_x_vector(option, A_QD->x, A_QD->n, A_QD->numRHS);
		}
	print_matrices(option, A_ref, A_QC, A_QD, A);				// Print matrices if desired
	print_stats(A_ref, A_QC, A_QD, option);					// Timing Statistics
	print_to_file(option, A_ref, A_QC, A_QD);				// Print results to file
	clear_all_memory(option, A, B, A_ref, A_QC, A_QD, B_ref, B_QC, B_QD);	// Clear memory
}