# include "REF_LU_config.h"

typedef struct REF_LU_Options		// Command line/problem options
{
	int REF;			// 1 if REF LU is chosen
	int Crout;			// 1 if Crout LU is chosen
	int Doolittle;			// 1 if Doolittle LU is chosen
	int numRHS;			// Number of right hand side vectors
	int dimension;			// Size of matrix
	int seed1;			// The first seed used to generate nonzero pattern
	int seed2;			// The second seed used to generate numbers
	int lower;			// Lower bound of randomly generated numbers
	int upper;			// Upper bound of randomly generated numbers
	double density;			// Density of the matrix to be generated
	int print;			// Print out matrix on screen
	int print2;			// 1 if solutions are output to a file
	int type;			// 1: Rational, 2: Double, 3: MPFR variable precision 
	int prec;			// MPFR Precision if necessary
	int check;			// Check if solution is correct
	string file;			// Filename
	string out;			// Name of output file containing x
} REF_LU_Options;

/* This function sets the defaults for the randomly generated matrices. These default parameters are defined in the REF_LU_config.h file */
void Options_set_defaults(REF_LU_Options* option)
{
	option->REF 		= DEFAULT_REF;	
	option->Crout		= DEFAULT_Crout;
	option->Doolittle 	= DEFAULT_Doolittle;
	option->numRHS		= DEFAULT_numRHS;
	option->dimension 	= DEFAULT_dimension;
	option->seed1 		= DEFAULT_seed1;
	option->seed2 		= DEFAULT_seed2;
	option->lower 		= DEFAULT_lower;
	option->upper 		= DEFAULT_upper;
	option->density		= DEFAULT_density;
	option->print 		= DEFAULT_print;
	option->check 		= DEFAULT_check;
	option->file 		= DEFAULT_name;
	option->print2		= DEFAULT_print2;
	option->type 		= DEFAULT_types;
	option->out 		= DEFAULT_out;
	option->prec 		= DEFAULT_prec;
}

/* This function will show the usage of the program. No input arguments */
void Options_show_usage()
{
	cerr << "\nUsage: \n\nOptions:\n\t-h,--help\n\t-a,--algs \t[123,13,2,..]\t>default=1<\n";
	cerr << "\t\t1:REF_LU, 2:Doolittle, 3:Crout\n\t(Can specify multiple w/o spaces in ascending order)\n\n";
	cerr << "\t-b,--num_RHS \t\t[int]\t\t>default=50<\t\t\n\t-c,--check_sol\t\t[1,0]\t\t>default=0<\t\n\t-d,--density\t\t[double]\t>default=1.0<\t\n";
	cerr << "\t-f,--file_out \t\t[string]\t>default_output.out<\t\n\t-n,--num_vars \t\t[int]\t\t>default=100<\n\t-p,--print \t\t[1,0]\t\t>default=0<\t\n";
	cerr << "\t-r,--range \t\t[int] [int]\t>default=-99,99<\n\t-s,--seeds \t\t[int] [int]\t>default=11,1400<\n" ;
	cerr << "\t-p2,--print2 \t\t[1,0] \t\t>default=1< \n\t-t,--type \t\t[1,2,3] \t>default=1< \n \t-o,--out\t\t[string] \t>default=default_x.out<\n";
	cerr << "\t-prec,--precision \t[int] \t\t>default=16<\n"<< endl; 
}

/* This function prints the options used. */
void print_options(REF_LU_Options* option)
{
	cout<< "\n==============Command Options===================\nAlgorithm Chosen: ";
	if (option->REF == 1)
		cout<< "REF LU ";
	if (option->Crout == 1)
		cout<< "Crout LU ";
	if (option->Doolittle == 1)
		cout<< "Doolittle LU";
	cout<< "\n\nDimension: \t\t"<< option->dimension;
	cout<< "\nNum RHS: \t\t"<< option->numRHS;
	cout<< "\nSeeds: \t\t\t" << option->seed1 << "," << option->seed2;
	cout<< "\nRange: \t\t\t[" << option->lower << "," << option->upper << "]";
	cout<< "\nDensity: \t\t" << option->density;
	cout<< "\nFilename: \t\t" << option->file;
	cout<< "\nOutput x name: \t\t" << option->out;
	cout<< "\nMPFR prec (if used): \t" << option->prec;
	cout<<"\n================================================\n";
}

/* This function processes the command line arguments */
int Options_process_command_line( int numArgs, char* args[], REF_LU_Options* options)
{
	int alg_choice;
	for (int i = 1; i < numArgs; ++i)
	{
		string arg = args[i]; 
		if (arg == "-h" || arg == "--help") 		// Check if requested help
		{
			Options_show_usage();			// Show the usage of the program
			return 0;				// Terminate the program
		}
		else if (arg == "-a" || arg == "--algs")	// Determine the algorithm to be used
		{
			alg_choice = atoi(args[++i]);
			if (alg_choice == 1 || alg_choice == 12 || alg_choice == 13 || alg_choice == 123) 	// REF LU
				options->REF = 1;
			if (alg_choice == 2 || alg_choice == 12 || alg_choice == 23 || alg_choice == 123)	// Doolittle
				options->Doolittle = 1;
			if (alg_choice == 3 || alg_choice == 13 || alg_choice == 23 || alg_choice == 123)	// Crout
				options->Crout = 1;
		}
		else if (arg == "-b" || arg == "--num_RHS")	// Determine the number of right hand side vectors to be generated
		{
			options->numRHS = atoi(args[++i]);
			if (options->numRHS < 0)
			{
				cerr << "\n*** Error: Number of RHS vectors must be non-negative ***\n";
				return 0;
			}
		}
		else if (arg == "-c" || arg == "--check_sol")	// Flag to check the solution
		{
			++i;
			if (atoi(args[i]) != 0 && atoi(args[i]) != 1)
			{
				cerr << "\n*** Error: Check input restricted to 0 or 1 ***\n";
				return 0;
			}
			options->check = 1;
		}
		else if (arg == "-d" || arg == "--density")	// Determine the density of the matrix
		{
			options->density = atof(args[++i]);
			if (options->density < 0 || options->density > 1)
			{
				cerr << "\n*** Error: Density must be in interval [0,1] ***\n";
				return 0;
			}
		}
		else if (arg == "-o" || arg == "--out")		// Determine the x output file
		{
			options->out = args[++i];
			if (options->file.size() < 1)
			{
				cerr << "\n*** Error: No file string input ***\n";
				return 0;
			}
		}
		else if (arg == "-f" || arg == "--file_out")	// Determine the output file
		{
			options->file = args[++i];
			if (options->file.size() < 1)
			{
				cerr << "\n*** Error: No file string input ***\n";
				return 0;
			}
		}
		else if (arg == "-n" || arg == "--num_vars")	// Determine the dimensions of the matrix
		{
			options->dimension = atoi(args[++i]);
			if (options->dimension < 1)
			{
				cerr << "\n*** Error: Number of variables must be positive ***\n";
				return 0;
			}
		}
		else if (arg == "-p" || arg == "--print")	// Determine if matrix is printed onto the screen
		{
			++i;
			if (atoi(args[i]) != 0 && atoi(args[i]) != 1)
			{
				cerr << "\n*** Error: Print input restricted to 0 or 1 ***\n";
				return 0;
			}
				options->print = 1;
		}
		else if (arg == "-p2" || arg == "--print2")	// Determine if x is printed to a file
		{
			++i;
			if (atoi(args[i]) != 0 && atoi(args[i]) != 1)
			{
				cerr << "\n*** Error: Print2 input restricted to 0 or 1 ***\n";
				return 0;
			}
				options->print2 = 1;
		}
		else if (arg == "-prec" || arg == "--precision")	// Determine MPFR precision
		{
			options->prec = atoi(args[++i]);
			if (options->prec < 2 )
			{
				cerr << "\n*** Error: Precision must be >= 2 ***\n";
				return 0;
			}
		}
		else if (arg == "-t" || arg == "--type")	// Determine how x is printed to a file
		{
			options->type = atoi(args[++i]);
			if (atoi(args[i]) < 1 || atoi(args[i]) >3 )
			{
				cerr << "\n*** Error: type input restricted to 1, 2, or 3 ***\n";
				return 0;
			}
				options->print2 = 1;
		}
		else if (arg == "-r" || arg == "--range")	// Determine the range of entries in the matrix
		{
			options->lower = atoi(args[++i]);
			options->upper = atoi(args[++i]);
			if(options->lower > options->upper)
			{
				cerr << "\n*** Error: Range lower bound must be less than or equal to upper bound ***\n";
				return 0;					
			}
		}
		else if (arg == "-s" || arg == "--seeds")	// Determine the seeds used to randomly generate matrices
		{
			options->seed1 = atoi(args[++i]);
			options->seed2 = atoi(args[++i]);
		}
		else
		{
			cerr << "\n*** Error: Unknown command line parameter: " << arg <<"\n";
			return 0;				// Terminate program
		}
	}
	return 1;
	cout << endl;
}