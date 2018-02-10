# include "REF_LU_config.h"
# include "REF_LU_Options.h"
# include "QLU.h"
# include "REF_LU.h"

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Printing to File/Screen--------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* Converts x to a double precision matrix for output */
void REF_convert_x_to_doub(mpq_t** A, double **x, int n, int numRHS)
{
	for (int i = 0; i < n; i++)
	{
		x[i] = new double [numRHS];
		for (int j = 0; j < numRHS; j++)
			x[i][j] = mpq_get_d(A[i][j]);
	}
}

/* Converts x to an MPFR matrix for output */
void REF_convert_x_to_mpfr(mpq_t** A, mpfr_t **x_mpfr, int n, int numRHS, int prec)
{
	int h;
	for (int i = 0; i < n; i++)
	{
		x_mpfr[i] = new mpfr_t [numRHS];
		for (int j = 0; j < numRHS; j++)
		{
			mpfr_init2(x_mpfr[i][j], prec);
			h = mpfr_set_q(x_mpfr[i][j], A[i][j], MPFR_RNDN);
		}
	}
}

/* This function prints summary statistics for REF LU code */
void print_stats(REF_mat* A, QLU_mat* C, QLU_mat* D, REF_LU_Options* option)
{
	if (option->numRHS > 0)				// Different set of summary stats depending on presence of RHS vectors
	{
		if (option->REF == 1)			// REF LU Stats
		{
			A->t_tot = A->t_f+A->t_s;	// Total factorization time
			cout<<"\nREF LU:\nFactorization Time (seconds): \t"<< A->t_f.count() ;
			cout<<"\nSolution Time (seconds): \t" << A->t_s.count() ;
			cout<<"\nTotal Time (seconds): \t\t" << A->t_tot.count() << "\n";
			if (option->check == 1)
			{
				if (A->check == 0)
					cout<<"REF LU is correct \n";
				else cout<<"REF LU is incorrect \n";
			}
		} 
		if (option->Crout == 1)			// Crout Stats
		{
			C->t_tot = C->t_f+C->t_s;	// Total factorization time
			cout<<"\nCrout LU:\nFactorization Time (seconds): \t"<< C->t_f.count();
			cout<<"\nSolution Time (seconds): \t" << C->t_s.count();
			cout<<"\nTotal Time (seconds): \t\t" << C->t_tot.count() << "\n";
			if (option->check == 1)
			{
				if (C->check == 0)
					cout<<"Crout LU is correct \n";
				else cout<<"Crout LU is incorrect \n";
			}
		}
		if (option->Doolittle == 1)		// Doolittle Stats
		{
			D->t_tot = D->t_f+D->t_s;	// Total factorization time
			cout<<"\nDoolittle LU:\nFactorization Time (seconds): \t"<< D->t_f.count();
			cout<<"\nSolution Time (seconds): \t" << D->t_s.count();
			cout<<"\nTotal Time (seconds): \t\t" << D->t_tot.count() << "\n";
			if (option->check == 1)
			{
				if (D->check == 0)
					cout<<"Doolittle LU is correct \n";
				else cout<<"Doolittle LU is incorrect \n";
			}
		}
	}
	else
	{
		if (option->REF == 1)			// REF Stats
		{
			cout<<"\nREF LU:\nFactorization Time (seconds): \t"<< A->t_f.count();
			cout<<"\nNo RHS Vectors";
			cout<<"\nTotal Time (seconds): \t\t" << A->t_f.count() << "\n";
		} 
		if (option->Crout == 1)			// Crout stats
		{
			cout<<"\nCrout LU:\nFactorization Time (seconds): \t"<< C->t_f.count();
			cout<<"\nNo RHS Vectors";
			cout<<"\nTotal Time (seconds): \t\t" << C->t_f.count() << "\n";
		}
		if (option->Doolittle == 1)		// Doolittle stats
		{
			cout<<"\nDoolittle LU:\nFactorization Time (seconds): \t"<< D->t_f.count();
			cout<<"\nNo RHS Vectors";
			cout<<"\nTotal Time (seconds): \t\t" << D->t_f.count() << "\n";
		}
	}	
}

/* Summary Stats for QMat */
void print_stats_QMat(REF_mat* A, REF_mat* Q, REF_LU_Options* option)
{
	if (option->numRHS > 0)			// Different set of summary stats depending on presence of RHS vectors
	{
		if (option->print2 == 1)		// Print x vector to file
			cout<<"\nSolution output to file named \t" << option->out;
		A->t_tot = A->t_f+A->t_s;	// Total factorization time
		cout<<"\nREF LU:\nFactorization Time (seconds): \t"<< A->t_f.count() ;
		cout<<"\nSolution Time (seconds): \t" << A->t_s.count() ;
		cout<<"\nTotal Time (seconds): \t\t" << A->t_tot.count() << "\n";
		if (option->check == 1)
		{
			if (A->check == 0)
				cout<<"REF LU is correct \n";
			else cout<<"REF LU is incorrect \n";
		}
		
		Q->t_tot = Q->t_f+Q->t_s;	// Total factorization time
		cout<<"\nQMat:\nConstruction Time (seconds): \t"<< Q->t_f.count() ;
		cout<<"\nSolution Time (seconds): \t" << Q->t_s.count() ;
		cout<<"\nTotal Time (seconds): \t\t" << Q->t_tot.count() << "\n";
		if (option->check == 1)
		{
			if (Q->check == 0)
				cout<<"Q_mat is correct \n";
			else cout<<"Q_mat is incorrect \n";
		}
	}
	else
	{
		cout<<"\nREF LU:\nFactorization Time (seconds): \t"<< A->t_f.count();
		cout<<"\nNo RHS Vectors";
		cout<<"\nTotal Time (seconds): \t\t" << A->t_f.count() << "\n";
		
		cout<<"\nQMat:\nConstruction Time (seconds): \t"<< Q->t_f.count();
		cout<<"\nNo RHS Vectors";
		cout<<"\nTotal Time (seconds): \t\t" << Q->t_f.count() << "\n";
	}	
}

/* This function prints an mpz matrix onto the screen */	
void print_mpz_mat(mpz_t** A, int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		cout<<"\n";
		for (int j = 0; j < n; j++)
			cout<<A[i][j]<<" ";
	}
}

/* This function prints an mpq matrix onto the screen */	
void print_mpq_mat(mpq_t** A, int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		cout<<"\n";
		for (int j = 0; j < n; j++)
			cout<<A[i][j]<<" ";
	}
}

/* This function prints an int matrix onto the screen */	
void print_int_mat(int** A, int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		cout<<"\n";
		for (int j = 0; j < n; j++)
			cout<<A[i][j]<<" ";
	}
}

/* This function prints out output from the main program */
void print_matrices(REF_LU_Options* option, REF_mat* A, QLU_mat* C, QLU_mat* D, int** orig)
{
	int i, j, n = A->n, m = A->m;
	print_options(option);			// Summary stats (i.e., what algorithm chosen)
	if (option->print == 1)			// If printing is desired
	{
		cout<<"\nInput matrix is: ";
		print_int_mat(orig, n, m);
		if (option->REF == 1)		// Print REF mat
		{
			cout<<"\nREF LU is :";
			print_mpz_mat(A->mat, n, m);
		}
		if (option->Doolittle == 1)	// Print Doolittle matrix
		{
			cout<<"\nDoolittle LU is: ";
			print_mpq_mat(D->mat, n, m);
		}
		if (option->Crout == 1)		// Print crout matrix
		{
			cout<<"\nCrout LU is: ";
			print_mpq_mat(C->mat, n, m);
		}
		if (option->numRHS > 0)		// Print x
		{
			cout<<"\nSolution is: ";
			print_mpq_mat(A->xq, A->numRHS, n);
			cout<<"\n";
		}
	}
}

/* Print the x vector to a file */
void print_x_vector(REF_LU_Options* option, mpq_t** A, int n, int numRHS)
{
	string filename = option->out;
	std::ofstream output(filename);
	if (option->type == 1)			// Print as rational numbers
	{
		output<<"Rational Solution:\n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < numRHS; j++)
			{
				output << A[i][j] << " ";
			}
			output<<"\n";
		}
	}
	else if (option->type == 2)		// Print in double precision
	{
		output<<"Double Solution:\n";
		double **x;
		x = new double *[n];
		REF_convert_x_to_doub(A, x, n, numRHS);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < numRHS; j++)
			{
				output << x[i][j] << " ";
			}
			output << "\n";
		}
		for (int i = 0; i < n; i++)
			delete[] x[i];
		delete[] x;
	}
	else 					// Print in fixed precision
	{
		output << "Fixed precision of size " << option->prec << "\n";
		mpfr_t** x;
		x = new mpfr_t *[n];
		REF_convert_x_to_mpfr(A, x, n, numRHS, option->prec);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < numRHS; j++)
			{
				char* out = NULL;
				mpfr_asprintf(&out, "%.*Rf", option->prec, x[i][j]);
				output << out << " ";
				mpfr_free_str(out);
			}
			output << "\n";
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < numRHS; j++)
			{
				mpfr_clear(x[i][j]);
			}
			delete[] x[i];
		}
		delete[] x;
	}
}

/* This function prints the data to file */
void print_to_file(REF_LU_Options* option, REF_mat* A, QLU_mat* C, QLU_mat* D)
{
	ofstream runStats;							// Output file
	char head[200];
	if (!std::ifstream(option->file))					// Print out the header of the file (if file is original)
		sprintf(head, "n\tDensity\tSeed1\tSeed2\t[lb,ub]\tQC(s)\tQC_s(s)\tQC_c\tQD(s)\tQD_s(s)\tQD_c\tIP(s)\tIP_s(s)\tIP_c\tQC/QD\tQC/QD_S\tQC/IP\tQC/IP_S\tQD/IP\tQD/IP_S\n ");
	else sprintf(head," ");
	runStats.open(option->file, ios::out | ios::app);
	runStats << head;
	runStats<<option->dimension <<  "|" << option->numRHS << "\t" << option->density << "\t" << option->seed1 	// Print out problem parameters
		<< "\t" << option->seed2 << "\t[" << option->lower << "," << option->upper << "]\t";
	if (option->Crout == 1)							// Print crout related stats
	{
		runStats <<  C->t_f.count() << "\t";		// IMPORTANT: Set precision as high as desired for number of decimal places/sig digits
		if (option->numRHS > 0) runStats << C->t_s.count() << "\t";
		else runStats << "\t";
		if (option->check == 1 && option->numRHS > 0) 
		{
			if (C->check == 0) runStats << 1;
			else runStats << 0;
		}
		runStats << "\t";
	}
	else runStats << "\t\t\t";
	if (option->Doolittle == 1)						// Print Doolittle related stats
	{
		runStats << D->t_f.count() << "\t" ;
		if (option->numRHS > 0) runStats << D->t_s.count() << "\t";
		else runStats << "\t";
		if (option->check == 1 && option->numRHS > 0) 
		{
			if (D->check == 0) runStats << 1;
			else runStats << 0;
		}
		runStats << "\t";
	}
	else runStats << "\t\t\t";
	if (option->REF == 1)							// Print REF LU related stats
	{
		runStats << A->t_f.count() << "\t";
		if (option->numRHS > 0) runStats << A->t_s.count() << "\t";
		else runStats << "\t";
		if (option->check == 1 && option->numRHS > 0)
		{
			if (A->check == 0)runStats << 1;
			else runStats << 0;
		}
		runStats << "\t";
	}
	else runStats << "\t\t\t";
	if (option->Doolittle == 1 && option->Crout == 1)			// Doolittle vs Crout
	{
		runStats << (double) C->t_f.count() / D->t_f.count() << "\t";
		if (option->numRHS > 0) runStats << (double) C->t_s.count() / D->t_s.count() << "\t";
		else runStats << "\t";
	}
	else runStats << "\t\t";
	if (option->REF == 1 && option->Crout == 1)				// REF vs Crout
	{
		runStats << (double) C->t_f.count() / A->t_f.count() << "\t";
		if (option->numRHS > 0) runStats << (double) C->t_s.count() / A->t_s.count() << "\t";
		else runStats << "\t";
	}
	else runStats << "\t\t";
	if (option->Doolittle == 1 && option->REF == 1)				// REF vs Doolittle
	{
		runStats << (double) D->t_f.count() / A->t_f.count() << "\t";
		if (option->numRHS > 0) runStats << (double) D->t_s.count() / A->t_s.count() << "\t";
		else runStats << "\t";
	}
	else runStats << "\t\t";
	runStats << "\n"; 
	runStats.close();
}

/* This function prints the data to file for QMat */
void QMat_print_to_file(REF_LU_Options* option, REF_mat* A, REF_mat* Q)
{
	ofstream runStats;		// Output file
	char head[150];
	if (!std::ifstream(option->file))
		sprintf(head, "n\tDensity\tSeed1\tSeed2\t[lb,ub]\tQ(s)\tQ_s(s)\tQ_c\tIP(s)\tIP_s(s)\tIP_c\tQ/IP\tQ/IP_S\n ");
	else sprintf(head," ");
	
	runStats.open(option->file, ios::out | ios::app);
	runStats << head;
	
	runStats<<option->dimension <<  "|" << option->numRHS << "\t" << option->density << "\t" << option->seed1 
		<< "\t" << option->seed2 << "\t[" << option->lower << "," << option->upper << "]\t";
	
	runStats <<  Q->t_f.count() << "\t";		// IMPORTANT: Set precision as high as desired for number of decimal places/sig digits
	if (option->numRHS > 0) runStats << Q->t_s.count() << "\t";	// Qmat related stats
	else runStats << "\t";
	if (option->check == 1 && option->numRHS > 0)
	{
		if (Q->check == 0)runStats << 1;
		else runStats << 0;
	}
	runStats << "\t";

	runStats << A->t_f.count() << "\t";				// REF LU Stats
	if (option->numRHS > 0) runStats << A->t_s.count() << "\t";
	else runStats << "\t";
	if (option->check == 1 && option->numRHS > 0)
	{
		if (A->check == 0)runStats << 1;
		else runStats << 0;
	}
	runStats << "\t";
	
	runStats << (double) Q->t_f.count() / A->t_f.count() << "\t";	// QMAT vs REF
	if (option->numRHS > 0) runStats << (double) Q->t_s.count() / A->t_s.count() << "\t";
	else runStats << "\t";
	runStats << "\n";
	runStats.close();
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Final Memory Management--------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* Free all memory to conclude a run of the code */
void clear_all_memory(REF_LU_Options* option, int** A, int** B, REF_mat* R, QLU_mat* C, QLU_mat* D, mpz_t** B_ref, mpq_t** B_QC, mpq_t** B_QD)
{
	GM_clear_matrix(A, R->m);		// Clear A
	if (option->numRHS > 0)			// Clear RHS vectors
	{
		GM_clear_matrix(B, R->m);
		if (option->REF == 1) REF_clear_mat(R->numRHS, R->n, B_ref);
		if (option->Crout == 1) QLU_clear_mat(C->numRHS, R->n, B_QC);
		if (option->Doolittle == 1) QLU_clear_mat(D->numRHS, R->n, B_QD);
	}
	if (option->REF == 1) REF_delete(R);	// Clear exact data structures
	else delete R;
	if (option->Crout == 1) QLU_delete(C);
	else delete C;
	if (option->Doolittle == 1) QLU_delete(D);
	else delete D;
	delete option;				// Clear command line option memory
}

/* This function clears memory for QMAT */
void clear_QMat_memory(REF_LU_Options* option, int** A, int** B, REF_mat* R, REF_mat* Q, mpz_t** B_ref, mpz_t** B_Q)
{
	int n = R->n;
	GM_clear_matrix(A,n);			// Clear A
	if (option->numRHS > 0)			// Clear RHS vectors
	{
		GM_clear_matrix(B,n);
		REF_clear_mat(R->numRHS,n,B_ref);
		REF_clear_mat(R->numRHS,n,B_Q);
	}
	REF_delete(R);				// Delete Exact data structures
	REF_delete(Q);
	delete option;				// Clear command line option memory
}