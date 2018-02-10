# include "REF_LU_config.h"

typedef struct QLU_mat			// Fully dense matrix stored as full precision integers for REF LU
{
	int n;				// Number of columns of the matrix 
	int m;				// Number of rows of the matrix
	int numRHS;			// Number of RHS vectors
	mpq_t** mat;			// Numerical matrix. Size of matrix is m*n
	mpq_t** x;			// Solution of linear system
	int* p;				// Row pointers
	duration<float> t_f;		// Factorization time
	duration<float> t_s;		// Solution time
	duration<float> t_tot;		// Total time
	int check;			// 0 if correct
} QLU_mat;

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Primary Routines: Memory management, initialization, etc-----------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function initializes an mpq matrix of size m*n. This function must be called for all mpq matrices created */
void QLU_initialize_matrix (int n, int m, mpq_t** A)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mpq_init(A[i][j]);
}

/* This function clears the memory associated with a mpq matrix. */
void QLU_clear_mat (int n, int m, mpq_t** A)
{
	for (int i = 0; i < m; i++)			// Clear the memory associated with each element
		for (int j = 0; j < n; j++)
			mpq_clear(A[i][j]);
	for (int i = 0; i < m; i++)			// Delete the vectors of A from memory
		delete[] A[i];	
	delete[] A;					// Delete A
}

/* This function allocates a QLU_mat matrix. */
void QLU_alloc (QLU_mat* A, int n, int m, int numRHS)
{
	A->n = n; A->m = m; A->numRHS = numRHS; 	// Assign Constants
	A->p = new int [m];
	A->mat = new mpq_t * [m];			// Declare memory for matrix
	for (int i = 0; i < m; i++)
	{
		A->mat[i] = new mpq_t [n];		// Declare matrix
		A->p[i] = i;				// Initialize row permutation
	}
	QLU_initialize_matrix(n, m, A->mat);
	if (numRHS > 0)					// Allocate memory for x vectors
	{
		A->x = new mpq_t* [m];
		for (int i = 0; i < m; i++)
			A->x[i] = new mpq_t [numRHS];
		QLU_initialize_matrix(numRHS, m, A->x);
	}
}

/* This function clears the memory associated with the QLU_mat A */
void QLU_delete(QLU_mat* A)
{
	delete[] A->p;
	QLU_clear_mat(A->n, A->m, A->mat);
	if (A->numRHS > 0) 
		QLU_clear_mat(A->numRHS,A->n,A->x); 
	delete A;
}

/* This function creates a copy of the input int matrix A */
void QLU_make_copy(int** A, mpq_t** A_c, int n, int m)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mpq_set_si(A_c[i][j], A[i][j], 1);		// Set A_c[i][j] = A[i][j] / 1
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------LU Factorization---------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function selects the kth pivot element */
int QLU_select_pivot(QLU_mat* A, int k)
{
	int pivot = -1, i, temp;
	
	/* Select pivot */
	i = k;								// Begin the search at A[k][k]
	while (i < A->m)
	{
		if (mpq_sgn(A->mat[A->p[i]][k]) != 0)			// Check if A[i][k] is nonzero
		{
			pivot = i;					// Select i as the pivot
			i = A->m+1;					// Break the loop
		}
		i+=1;
	}
	if (pivot == -1) return pivot;					// Matrix is singular
	
	/* Swap row locations */
	temp = A->p[pivot];
	A->p[pivot] = A->p[k];
	A->p[k] = temp;
	return 0;
}

/* This function performs Crout LU. If A[k][k] doesn't exist, the next nonzero is chosen in column k */
int QLU_crout_LU_pivots(QLU_mat* A)
{
	steady_clock::time_point t_f_begin = steady_clock::now();
	int i, j, k, check;
	mpq_t temp, sum;
	mpq_inits(temp, sum, NULL);
	for (k = 0; k < min(A->m,A->n); k++)
	{
		for (i = k; i < A->m; i++) 
		{
			mpq_set_ui(sum, 0,1);						// Default sum to 0
			for (j = 0; j < k; j++) 
			{
				mpq_mul(temp, A->mat[A->p[i]][j], A->mat[A->p[j]][k]);	// temp = A[i][j] * A[j][k]
				mpq_add(sum, sum, temp);				// sum = sum + temp
			}
			mpq_sub(A->mat[A->p[i]][k], A->mat[A->p[i]][k], sum);		// A[i][k] = A[i][k] - sum
		}
		check = QLU_select_pivot(A,k);						// Select the pivot element
		if (check == -1)							// Check for errors
		{
			cerr<<"\n*****ERROR Matrix is singular!****";
			return 0;							// Failure
		}
		for (i = k; i < A->m; i++) 
		{
			mpq_set_ui(sum, 0,1);						// Default sum to 0
			for (j = 0; j < k; j++) 
			{
				mpq_mul(temp, A->mat[A->p[k]][j], A->mat[A->p[j]][i]);	// temp = A[k][j]*A[j][i]
				mpq_add(sum, sum, temp);				// sum = sum + temp
			}
			if (i != k)
			{
				mpq_sub(A->mat[A->p[k]][i], A->mat[A->p[k]][i], sum);	// A[k][i] = A[k][i] - sum
				mpq_div(A->mat[A->p[k]][i], A->mat[A->p[k]][i], A->mat[A->p[k]][k]);	// A[k][i] = A[k][i] / A[k][k]
			}
		}
	}
	mpq_clear(temp); mpq_clear(sum);						// Free memory with temp and sum
	steady_clock::time_point t_f_end = steady_clock::now();
	A->t_f = duration_cast<duration<float>>(t_f_end - t_f_begin);			// Crout Factorization time
	return 1;									// Success
}

/* This function performs Doolittle LU. If A[k][k] doesn't exist, the next nonzero is chosen in column k */
int QLU_doolittle_LU_pivots(QLU_mat* A)
{
	steady_clock::time_point t_f_begin = steady_clock::now();
	int i, j, k, check;
	mpq_t temp; mpq_init(temp);								// Temporary product term
	for (k = 0; k < min(A->n,A->m); k++)
	{
		check = QLU_select_pivot(A,k);							// Select the kth pivot
		if (check == -1)
		{
			cerr<<"\n*****ERROR Matrix is singular!****";
			return 0;								// Failure
		}
		for (i = k+1; i < A->m; i++)
		{
			mpq_div(A->mat[A->p[i]][k], A->mat[A->p[i]][k], A->mat[A->p[k]][k]);	// A[i][k] = A[i][k] / A[k][k]
			for (j = k+1; j < A->n; j++)
			{
				mpq_mul(temp, A->mat[A->p[i]][k], A->mat[A->p[k]][j]);		// temp = A[i][k] * A[k][j]
				mpq_sub(A->mat[A->p[i]][j], A->mat[A->p[i]][j], temp);		// A[i][j] = A[i][j] - temp
			}
		}
	}
	mpq_clear(temp);									// Free memory associated with temp
	steady_clock::time_point t_f_end = steady_clock::now();
	A->t_f = duration_cast<duration<float>>(t_f_end - t_f_begin);				// Doolittle Factorization time
	return 1;										// Success
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Forward/Backward Sub-----------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* QLU Forward Substitution */
void QLU_forward_sub(mpq_t** Y, mpq_t** B, QLU_mat* A, bool triag)
{
	int i, j, k;
	mpq_t prod; mpq_init(prod);								// Temporary product term
	for (k = 0; k < A->numRHS; k++)								// Iterate accross all RHS vectors
		for (i = 0; i < A->n; i++)
		{
			mpq_set(Y[A->p[i]][k], B[A->p[i]][k]);					// Y[i][k] = B[i][k]
			for (j = 0; j < i; j++)
			{
				mpq_mul(prod, A->mat[A->p[i]][j], Y[A->p[j]][k]);		// prod = A[i][j] * Y[j][k]
				mpq_sub(Y[A->p[i]][k], Y[A->p[i]][k], prod);			// Y[i][k] = Y[i][k] - prod
			}
			if (!triag)
				mpq_div(Y[A->p[i]][k], Y[A->p[i]][k], A->mat[A->p[i]][i]);	// Y[i][k] = Y[i][k] / A[i][i]
		}
	mpq_clear(prod);									// Free memory associated with prod
}

/* QLU Forward Substitution */
void QLU_backward_sub(mpq_t** Y, QLU_mat* A, bool triag)
{
	int i, j, k;
	mpq_t prod; mpq_init(prod);									// Intermediate prod term
	for (k = 0; k < A->numRHS; k++)									// Iterate across all RHS vectors
		for (i = A->n-1; i >= 0; i--)
		{
			mpq_set(A->x[A->p[i]][k], Y[A->p[i]][k]);					// x[i][k] = Y[i][k]
			for(j = i+1; j < A->n; j++)
			{
				mpq_mul(prod,  A->mat[A->p[i]][j], A->x[A->p[j]][k]);			// prod = A[i][j] * A[j][k]
				mpq_sub(A->x[A->p[i]][k], A->x[A->p[i]][k], prod);			// A[i][k] = A[i][k] - prod
			}
			if(!triag)
				mpq_div(A->x[A->p[i]][k], A->x[A->p[i]][k], A->mat[A->p[i]][i]);	// A[i][k] = A[i][k] / A[i][i]
		}
	mpq_clear(prod);										// Free memory associated with prod
}

/* Solve the linear system */
void QLU_solve(QLU_mat* A, mpq_t** B, int** B_orig, bool crout)
{
	steady_clock::time_point t_s_begin = steady_clock::now();
	QLU_make_copy (B_orig, B, A->numRHS, A->n);			// Create a copy of B into rational arithmetic
	mpq_t** Y;							// Declare intermediate array Y
	Y = new mpq_t* [A->n];
	for (int i = 0; i < A->n; i++)
		Y[i] = new mpq_t [A->numRHS];
	QLU_initialize_matrix(A->numRHS, A->n, Y);
	if (crout)							// Crout solve
	{
		QLU_forward_sub(Y, B, A, false);
		QLU_backward_sub(Y, A, true);
	}
	else								// Doolittle Solve
	{
		QLU_forward_sub(Y, B, A, true);
		QLU_backward_sub(Y, A, false);
	}
	QLU_clear_mat(A->numRHS, A->n, Y);
	steady_clock::time_point t_s_end = steady_clock::now();
	A->t_s = duration_cast<duration<float>>(t_s_end - t_s_begin);	// Solution time
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Solution Verification----------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function checks to make sure that A*x=b */
int QLU_check_x(QLU_mat* A, mpq_t** B, mpq_t** orig)
{
	int check, i, j, k, n = A->n, numRHS = A->numRHS, m = A->m, y = 0;
	mpq_t temp; mpq_init(temp);						// Temp product
	mpq_t** AX;								// Declare Ax
	AX = new mpq_t* [m];
	for (i = 0; i < m; i++)
	{
		AX[i] = new mpq_t[numRHS];
		for (j = 0; j < numRHS; j++)
			mpq_init(AX[i][j]);
	}
	
	for (i = 0; i < m; i++)							// Ax = orig * x
		for (j = 0; j < numRHS; j++)
			for (k = 0; k < n; k++)
			{
				mpq_set(temp,orig[A->p[i]][k]);			// temp = orig[i][k]
				mpq_mul(temp, temp, A->x[A->p[k]][j]);		// temp = temp * x[k][j]
				mpq_add(AX[A->p[i]][j], AX[A->p[i]][j], temp);	// Ax[i][j] = Ax[i][j] + temp
			}

	for (i = 0; i < m; i++)							// Compare
		for (j = 0; j < numRHS; j++)
		{
			check = mpq_equal(AX[i][j], B[i][j]);			// Does Ax[i][j] = B[i][j]?
			if (check == 0) y = 1;					// Not equal
		}
	mpq_clear(temp);							// Clear memory associated with temp
	QLU_clear_mat(numRHS,n,AX);						// Clear memory of Ax
	return y;
}