#include "REF_LU_config.h"

/* This function initializes the matrix A  */
void GM_initialize_matrix (int** A, int n, int m)
{
	for (int i = 0; i < m; i++)
		A[i] = new int [n];
}

/* Initializes a matrix of ints to zero */
void GM_init_mat_zeros (int** A, int m, int n)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = 0;
}

/* This function clears the memory of the int matrix A */
void GM_clear_matrix(int** A, int m)
{
	for (int i = 0; i < m; i++)
		delete[] A[i];
	delete[] A;
}

/* This function randomly generates the entries in A */
void GM_generate_matrix(int** A, int n, int m, int lb, int ub, int seed1, int seed2, double density, bool diag)
{
	int i, j, index1, index2, temp, *row_index_rand;
	double frac_pos;							// Fraction of entries expected to be positive
	int nnz = ceil(density*n*m);						// Number of nonzeros
	GM_init_mat_zeros(A, m, n);						// Initialize A to all zeros
	
	/* Determine the ratio of positive to negative entries */
	if (lb >= 0) frac_pos = 1;
	else if (ub <= 0) frac_pos = 0;
	else frac_pos = (double) ub/(ub-lb);
	
	/* First, we ensure the matrix is nonsingular by placing a nonzero element along each column without repeating rows */
	srand(seed1);
	row_index_rand = new int [m];		
	for (i = 0; i < m; i++)							// Populate row_index_rand with the natural order first
		row_index_rand[i] = i;
	
	if (!diag)								// Generate the nonzeros to ensure nonsingularity RANDOMLY
		for (i = 0; i < 2*m; i++)
		{
			index1 = rand() % m;					// Generate the random indices
			index2 = rand() % m;
			temp = row_index_rand[index1];				// Swap
			row_index_rand[index1] = row_index_rand[index2];
			row_index_rand[index2] = temp;
		}
	srand(2*seed1);								// Generate new set of random numbers
	for (i = 0; i < min(n,m); i++)						// Populate the entries
	{
		index1 = row_index_rand[i];
		if ( (double) rand() / (RAND_MAX) < frac_pos)			// Check if this number should be positive or negative
			A[index1][i] = 1;
		else
			A[index1][i] = -1;
		nnz-=1;								// Decrement nnz
	}
	delete[] row_index_rand;						// Clear memory associated with row_index_rand

	/* Symbolically place the new values */
	while (nnz > 0)
	{
		index1 = rand() % m;						// Randomly generate the indices
		index2 = rand() % n;
		if ( A[index1][index2] == 0)					// Place a symbolic entry
		{
			if ( (double) rand() / (RAND_MAX) < frac_pos)		// Check if number should be positive or negative
				A[index1][index2] = 1;
			else
				A[index1][index2] = -1;
			nnz -= 1;						// Decrement nnz
		}
	}
	
	/* Generate the numerical entries */
	for (j = 0; j < n; j++)
		for (i = 0; i < m; i++)
		{ 
			if (A[i][j] < 0)					// Populate with a negative value
				A[i][j] = -1 * ( rand() % (-1*lb)+1);
			else if (A[i][j] > 0)					// Populate with positive value
				A[i][j] = rand() % ub + 1;
		}
}