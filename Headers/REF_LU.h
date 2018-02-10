# include "REF_LU_config.h"

typedef struct REF_mat			// Fully dense matrix stored as full precision integers for REF LU
{
	int n;				// Number of columns of the matrix 
	int m;				// Number of rows of the matrix
	int numRHS;			// Number of RHS vectors
	mpz_t** mat;			// Numerical matrix. Size of matrix is m*n
	mpz_t** x;			// Numerator to solution of linear system. Size is m*numRHS
	mpq_t** xq;			// Full solution of linear system
	int* p;				// Row pointers
	duration<float> t_f;		// Factorization time
	duration<float> t_s;		// Solution time
	duration<float> t_tot;		// Total time
	int check;			// 0 if correct
} REF_mat;

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Primary Routines: Memory management, initialization, etc-----------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/


/* This function initializes an mpz matrix of size m*n. This function must be called for all mpz matrices created */
void REF_initialize_matrix (int n, int m, mpz_t** A)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mpz_init(A[i][j]);
}

/* This function clears the memory associated with a mpz matrix */
void REF_clear_mat (int n, int m, mpz_t** A)
{
	for (int i = 0; i < m; i++)			// Clear the memory associated with each element
		for (int j = 0; j < n; j++)
			mpz_clear(A[i][j]);
	for (int i = 0; i < m; i++)			// Delete the vectors of A from memory
		delete[] A[i];
	delete[] A;					// Delete A
}

/* This function clears the memory associated with a mpq matrix */
void REF_clear_mpq_mat (int n, int m, mpq_t** A)
{
	for (int i = 0; i < m; i++)			// Clear memory associated with each element
		for (int j = 0; j < n; j++)
			mpq_clear(A[i][j]);
	for (int i = 0; i < m; i++)			// Delete the vectors of A
		delete[] A[i];
	delete[] A;					// Delete A
}

/* This function allocates a REF_mat matrix */
void REF_alloc (REF_mat* A, int n, int m, int numRHS)
{
	A->n = n; A->m = m; A->numRHS = numRHS;		// Set constants
	A->p = new int [m];
	A->mat = new mpz_t *[m];
	for (int i = 0; i < m; i++)			// Allocate memory for A matrix and A->p
	{
		A->mat[i] = new mpz_t [n];
		A->p[i] = i;
	}
	REF_initialize_matrix(n, m, A->mat);		// Initialize the mpz object A->mat
	if (numRHS > 0)					// Allocate RHS if they exist
	{
		A->x = new mpz_t* [m];
		A->xq = new mpq_t* [m];
		for (int i = 0; i < n; i++)
		{
			A->x[i] = new mpz_t [numRHS];
			A->xq[i] = new mpq_t [numRHS];
			for (int j = 0; j < numRHS; j++) 
				mpq_init(A->xq[i][j]);
		}
		REF_initialize_matrix(numRHS,m,A->x);
	}
}

/* This function clears the memory associated with the REF_mat A */
void REF_delete(REF_mat* A)
{
	delete[] A->p;				// Delete the row permutation vector
	REF_clear_mat(A->n, A->m, A->mat);	// Delete the A->mat structure
	if (A->numRHS > 0) 			// Delete RHS dependent matrices
	{
		REF_clear_mat(A->numRHS, A->n, A->x); 
		REF_clear_mpq_mat(A->numRHS, A->n, A->xq);
	}
	delete A;				// Delete A
}

/* This function creates a copy of the input int matrix A*/
void REF_make_copy(int** A, mpz_t** A_c, int n, int m)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mpz_set_si(A_c[i][j],A[i][j]);			// Set A_c[i][j] = A[i][j]
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------LU Factorization---------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function selects the kth pivot element and places it in position A->mat[k][k] */

int REF_select_pivot(REF_mat* A, int k)
{
	int pivot = -1, i, temp;
	
	/* Select pivot */
	i = k;								// Begin the search at A[k][k]
	while (i < A->m)
	{
		if (mpz_sgn(A->mat[A->p[i]][k]) != 0)			// Check if A[i][k] is nonzero
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

/* This function performs REF LU with pivoting. If A[k][k] does not exist, the next nonzero is chosen in column k */

int REF_LU_Factorization_pivots(REF_mat* A)
{
	steady_clock::time_point t_f_begin = steady_clock::now();
	int i, j, k, n = A->n, m = A->m, check;
	mpz_t temp; mpz_init(temp);									// Temporary product term
	for (k = 0; k < min(m,n); k++)									// Iterate accross all IPGE iterations
	{
		check = REF_select_pivot(A, k);								// Select the kth pivot. Default is A[k][k]
		if (check == -1)
		{
			cerr<<"\n*****ERROR Matrix is singular!****";
			return 0;									// Failure
		}
		for (i = k+1; i < m; i++)								// Operations on the sub matrix
			for (j = k+1; j < n; j++)
			{
				mpz_mul(A->mat[A->p[i]][j], A->mat[A->p[i]][j], A->mat[A->p[k]][k]); 	// A[i][j] = A[i][j] * A[k][k]
				mpz_mul(temp, A->mat[A->p[i]][k], A->mat[A->p[k]][j]);			// temp = A[i][k] * A[k][j]
				mpz_sub(A->mat[A->p[i]][j], A->mat[A->p[i]][j], temp);			// A[i][j] = A[i][j] - temp
				if (k > 0)								// A[i][j] = A[i][j] / A[k-1][k-1]
					mpz_divexact(A->mat[A->p[i]][j], A->mat[A->p[i]][j], A->mat[A->p[k-1]][k-1]); 
			}
	}
	mpz_clear(temp);										// Free memory associated with temp
	steady_clock::time_point t_f_end = steady_clock::now();
	A->t_f = duration_cast<duration<float>>(t_f_end - t_f_begin);					// Compute Factorization time
	return 1;											// Success!
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Forward and Back Substitution--------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function performs REF forward substitution and stores the solution in the intermediate array Y */
void REF_forward_sub(mpz_t** Y, mpz_t** B, REF_mat* A, int numRHS, int n)
{
	mpz_t temp; mpz_init(temp);								// Intermediate product term
	for (int k = 0; k < numRHS; k++)							// Iterate accross all RHS vectors
	{
		mpz_set(Y[A->p[0]][k],B[A->p[0]][k]);						// Y[0][k] = B[0][k]
		for (int i = 1; i < n; i++)
		{
			mpz_set(Y[A->p[i]][k], B[A->p[i]][k]);					// Initialization: Y[i][k] = B[i][k]
			for (int r = 0; r < i; r++)						// r to follow convention of paper REF Factorization
			{
				mpz_mul(Y[A->p[i]][k], A->mat[A->p[r]][r], Y[A->p[i]][k]);	// Y[i][k] = A[r][r]*Y[i][k]
				mpz_mul(temp, A->mat[A->p[i]][r], Y[A->p[r]][k]);		// temp = A[i][r] * Y[r][k]
				mpz_sub(Y[A->p[i]][k], Y[A->p[i]][k], temp);			// Y[i][k] = Y[i][k] - temp
				if (r > 0)							// Y[i][k] = Y[i][k] / A[r-1][r-1]
					mpz_divexact(Y[A->p[i]][k], Y[A->p[i]][k], A->mat[A->p[r-1]][r-1]);
			}
		}
	}
	mpz_clear(temp);									// Free memory associated with temp
}

/* This function performs REF Backward substitution and stores the solution in A->x */
void REF_back_sub(REF_mat* A, mpz_t** Y, int numRHS, int n)
{
	mpz_t prod, det;									// Temporary mpz objects
	mpz_inits(prod, det, NULL);
	mpz_set(det, A->mat[A->p[n-1]][n-1]);							// Set the determinant
	for (int k = 0; k < numRHS; k++)							// Iterate accross all right hand side vectors
		for (int i = n-1; i>=0; i--)
		{
			mpz_mul(A->x[A->p[i]][k], det, Y[A->p[i]][k]);				// Scaling: x[i][k] = det*Y[i][k]
			for (int j = i+1; j < n; j++)
			{
				mpz_mul(prod, A->mat[A->p[i]][j], A->x[A->p[j]][k]);		// prod = A[i][j]*A[j][k]
				mpz_sub(A->x[A->p[i]][k], A->x[A->p[i]][k], prod);		// x[i][k] = x[i][k] - prod
			}
			mpz_divexact(A->x[A->p[i]][k], A->x[A->p[i]][k], A->mat[A->p[i]][i]); 	// x[i][k] = x[i][k] / A[i][i]
		}
	mpz_clear(prod); mpz_clear(det);							// Free memory associated with prod and det
}

/* This function gets the full solution to the linear system in rational arithmetic */
void REF_get_soln(REF_mat* A, int n, int numRHS)
{
	mpq_t det; mpq_init(det);				// Set the determinant
	mpq_set_z(det, A->mat[A->p[n-1]][n-1]);
	
	for (int k = 0; k < numRHS; k++)			// Set xq = x / det
		for (int i = 0; i < n; i++)
		{
			mpq_set_z(A->xq[i][k],A->x[i][k]);	// convert to rational arithmetic xq[i][k] = x[i][k]
			mpq_div(A->xq[i][k],A->xq[i][k],det);	// xq[i][k] = xq[i][k] / det
		}
	mpq_clear(det);						// Free memory associated with det
}

/* This function solves the linear system using REF LU */
void REF_solve(REF_mat* A, mpz_t** B_ref, int** B)
{
		steady_clock::time_point t_s_begin = steady_clock::now();
		REF_make_copy (B, B_ref, A->numRHS, A->n);			// Create a copy of B into integer arithmetic
		mpz_t** Y;							// Allocate intermediate array Y
		Y = new mpz_t* [A->n];
		for (int i = 0; i < A->n; i++)
			Y[i] = new mpz_t [A->numRHS];
		REF_initialize_matrix(A->numRHS, A->n, Y);			// Initialize Y
		REF_forward_sub(Y, B_ref, A, A->numRHS, A->m);			// Perform Forward Sub
		REF_back_sub(A, Y, A->numRHS, A->n);				// Perform backwards sub
		REF_clear_mat(A->numRHS, A->n, Y);				// Free memory associated with Y
		REF_get_soln(A, A->n, A->numRHS);				// Obtian rational solution
		steady_clock::time_point t_s_end = steady_clock::now();
		A->t_s = duration_cast<duration<float>>(t_s_end - t_s_begin);	// REF_LU Solution time
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Edmond's Q Matrices------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function implements Edmond's Q matrices. Note: No Pivoting is performed */
void Q_Mat_get_adjunct(REF_mat* A, mpz_t det)
{
	steady_clock::time_point t2_f_begin = steady_clock::now();
	mpz_t temp, piv, prevpiv;					// Intermediate terms
	mpz_inits(temp, piv, prevpiv, NULL);
	mpz_set_ui(prevpiv,1);
	for (int k = 0 ; k < A->m; k++)
	{
		for (int i = 0; i < A->m; i++)
			if (i != k)
				for (int j = 0; j < A->n; j++)
					if (j != k)
					{
						mpz_mul(A->mat[A->p[i]][j], A->mat[A->p[i]][j], A->mat[A->p[k]][k]);	// A[i][j] = A[i][j] * A[k][k]
						mpz_mul(temp, A->mat[A->p[k]][j], A->mat[A->p[i]][k]);			// temp = A[k][j] * A[i][k]
						mpz_sub(A->mat[A->p[i]][j], A->mat[A->p[i]][j], temp);			// A[i][j] = A[i][j] - temp
						if (k > 0)
							mpz_divexact(A->mat[A->p[i]][j], A->mat[A->p[i]][j], prevpiv);	// A[i][j] = A[i][j] / A[k-1][k-1]
					}
		
		mpz_set(piv, A->mat[A->p[k]][k]);			// Set the current pivot
		for (int i = 0; i < A->m; i++)				// Negate all elements in the pivot column except the pivot
			if (i != k)
				mpz_neg(A->mat[A->p[i]][k], A->mat[A->p[i]][k]);
		
		mpz_set(A->mat[A->p[k]][k],prevpiv);			// Set A[k][k] = previous pivot
		mpz_set(prevpiv, piv);					// Update previous pivot
	}
	mpz_set(det, prevpiv);						// Set the determinant
	mpz_clear(temp); mpz_clear(piv); mpz_clear(prevpiv);		// Free memory
	steady_clock::time_point t2_f_end = steady_clock::now();
	A->t_f = duration_cast<duration<float>>(t2_f_end - t2_f_begin);	// Compute Factorization time
}

/* Multiply A->mat * b for Ax */
void Q_Mat_solve(REF_mat* A, mpz_t** B, int numRHS)
{
	mpz_t temp;		// Temporary calculations constant
	mpz_init(temp);
	for (int i = 0; i < A->m; i++)
		for (int j = 0; j < numRHS; j++)
			for (int k = 0; k < A->n; k++)
			{
				mpz_mul(temp, A->mat[A->p[i]][k], B[A->p[k]][j]);	// temp = A[i][k] * B[k][j]
				mpz_add(A->x[A->p[i]][j], A->x[A->p[i]][j], temp);	// x[i][j] = x[i][j] + temp
			}
	mpz_clear(temp);	// Free memory from temp
}

/* Finalize x by dividing by determinant */
void Q_Mat_get_x(REF_mat* A, int numRHS, mpz_t det)
{
	mpq_t det2;							// Rational version of determinant
	mpq_init(det2);
	mpq_set_z(det2, det);						// Convert det to rational
	for (int k = 0; k < numRHS; k++)
		for (int i = 0; i < A->m; i++)
		{
			mpq_set_z(A->xq[i][k], A->x[i][k]);		// Convert x into rational arithmetic
			mpq_div(A->xq[i][k], A->xq[i][k], det2);	// x[i][k] = x[i][k] / det2
		}
	mpq_clear(det2);						// Free memory
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Solution Check-----------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* This function checks to make sure that A*x=b */
int REF_check_x(REF_mat* A, mpz_t** B, mpz_t** orig)
{
	int check, i, j, k, n = A->n, numRHS = A->numRHS, m = A->m, y = 0;
	mpq_t temp; mpq_init(temp);
	mpq_t** AX;								// Declare Ax matrix
	AX = new mpq_t* [m];
	for (i = 0; i < m; i++)
	{
		AX[i] = new mpq_t [numRHS];
		for (j = 0; j < numRHS; j++)
			mpq_init(AX[i][j]);
	}
	
	for (i = 0; i < m; i++)							// Ax = orig * A->xq
		for (j = 0; j < numRHS; j++)
			for (k = 0; k < n; k++)
			{
				mpq_set_z(temp,orig[A->p[i]][k]);		// temp = orig[i][k]
				mpq_mul(temp, temp, A->xq[A->p[k]][j]);		// temp = temp * xq[k][j]
				mpq_add(AX[A->p[i]][j], AX[A->p[i]][j], temp);	// AX[i][j] = AX[i][j] + temp
			}

	for (i = 0; i < m; i++)							// Compare Ax with B
		for (j = 0; j < numRHS; j++)
		{
			mpq_set_z(temp, B[i][j]);				// temp = B[i][j]
			check = mpq_equal(AX[i][j], temp);			// nonzero if equal
			if (check == 0) y = 1;					// Not equal
		}
	mpq_clear(temp);							// Free memory associated with temp
	REF_clear_mpq_mat(numRHS,n,AX);						// Free memory associated with AX
	return y;
}