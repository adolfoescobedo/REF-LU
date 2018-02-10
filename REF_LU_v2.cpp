# include "./Headers/Generate_Matrix.h"
# include "./Headers/Output.h" 

/* REF_LU_v2 This code performs the REF LU factorization with an input file name. Please see README.txt for information */
int main (int argc, char* argv[])
{
	/* Declare variables and data structures */	
	int i, j, m, n, i2, j2, x2, check, nnz, alg, type, prec;
	REF_mat* A_ref;								// Data structure for the REF Mat 
	QLU_mat *A_QC, *A_QD;							// Data structure for QLU Crout/Doolittle
	mpz_t **copy, **B_ref;							// Copies and exact versions of B
	mpq_t **copy_q, **B_QC, **B_QD;
	
	A_ref = new REF_mat;							// Declare memory for A
	A_QC = new QLU_mat;
	A_QD = new QLU_mat;
		
	REF_LU_Options* option;							// The command line argument vector
	option = new REF_LU_Options;						// Allocate memory for option
	Options_set_defaults(option);						// Set default values
    
	/* Read in the input matrix */
	string matrix, rhs;
	if (argc >= 3 )
	{
		matrix = argv[1];						// Read in matrix name
		rhs = argv[2];							// Read in RHS vector name
		for (i = 3; i < argc; i++)
		{
			string arg = argv[i];
			if (arg == "-a" || arg == "--alg")
			{
				alg = atoi(argv[++i]);				// Read in type of algorithm
				/* Determine which algorithm to use */
				if (alg == 1 || alg == 12 || alg == 123 || alg == 13)
					option->REF = 1;
				if (alg == 2 || alg == 12 || alg == 123 || alg == 23)
					option->Doolittle = 1;
				if (alg == 3 || alg == 13 || alg == 123 || alg == 23)
					option->Crout = 1;
			}
			else if (arg == "-t" || arg == "--type")
				option->type = atoi(argv[++i]);
			else if (arg == "-prec" || arg == "--precision")
				option->prec = atoi(argv[++i]);
			else if (arg == "-o" || arg == "--out")
				option->out = argv[++i];
			else
			{
				cout<<"\nUnknown command line input: "<<arg << "\n";
				return 0;
			}
		}
	}
	else									// Shows usage if not EXACTLY 3 command line arguments
	{
		cout<<"\nERROR: This program assumes input is a matrix and rhs vector";
		cout<<"\nIf a randomly generated matrix is desired, use REF_LU";
		cout<<"\nUSAGE: ./REF_LU_v2.exe MATRIX_FILE_NAME RHS_FILE_NAME\nFollowed by:";
		cout<<"\n\t-a,--algs \t\t[123,13,2,..]\t>default=1<\n\t-t,--type \t\t[1,2,3] \t>default=1< \n \t-o,--out\t\t[string] \t>default=default_x.out<\n";
		cout<<"\t-prec,--precision \t[int] \t\t>default=16<\n"<< endl;
		return 0;
	}
	
	/* Read in the matrix. Assumed to be in matrix market format. 
	   Refer to http://math.nist.gov/MatrixMarket/formats.html for info on matrix market */
	   
	option->numRHS = 1;							// Assumed to be using only 1 RHS vector
	string line;
	ifstream in_file;
	in_file.open(matrix);							// Open the matrix file
	getline(in_file,line);
	istringstream line_ss(line);						// Read in dimensions of matrix and number of nonzeros
	line_ss >> m;
	line_ss >> n;
	line_ss >> nnz;

	int **A, **B;								// The input matrix and RHS vector (remains int** for ease)
	A = new int* [n];
	GM_initialize_matrix(A, m, n);						// Allocate memory for A
	GM_init_mat_zeros(A, m, n);						// Set A to all zeros
	B = new int* [m];
	GM_initialize_matrix(B, option->numRHS, m);				// Allocate B
	GM_init_mat_zeros(B, m, option->numRHS);	
	for (i = 0; i < nnz; i++)						// Iterate accross all nonzeros
	{
		getline(in_file,line);
		istringstream line_ss(line);					// Input is i j x
		line_ss >> i2;
		line_ss >> j2;
		line_ss >> x2;
		A[i2-1][j2-1] = x2;						// Conversion from 1 based to 0 based
	}
	in_file.close();

	ifstream in_file2;							// RHS vector
	in_file2.open(rhs);
	for (i = 0; i < n; i++)							// Read in the values from file
	{	
		getline(in_file2,line);						// Input is strictly number. Assumes only 1 rhs vector
		istringstream line_ss2(line);
		line_ss2 >> B[i][0];
	}
	in_file2.close();
	option->check = 1;							// We will check solutions for input file

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
			B_QD = new mpq_t * [m];					// Allocate memory for Doolittle RHS vectors
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
			copy = new mpz_t * [m];					// Allocate memory for copy
			for (i = 0; i < n; i++)	
				copy[i] = new mpz_t [n];
			REF_initialize_matrix(n, m, copy);
			REF_make_copy(A, copy, n, m);				// Create copy
			A_ref->check = REF_check_x(A_ref, B_ref, copy);
			REF_clear_mat(n, m, copy);
		}
		if (option->Doolittle == 1 || option->Crout == 1)		// Only 1 rational copy is necesary for solution verification
		{
			copy_q = new mpq_t * [m];				// Allocate memory for copy_q
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
	if (option->REF == 1)
		print_x_vector(option, A_ref->xq, A_ref->n, A_ref->numRHS);
	else if (option->Crout == 1)
		print_x_vector(option, A_QC->x, A_QC->n, A_QC->numRHS);
	else if (option->Doolittle == 1)
		print_x_vector(option, A_QD->x, A_QD->n, A_QD->numRHS);
	print_stats(A_ref, A_QC, A_QD, option);					// Timing Statistics
	clear_all_memory(option, A, B, A_ref, A_QC, A_QD, B_ref, B_QC, B_QD);	// Free memory
}