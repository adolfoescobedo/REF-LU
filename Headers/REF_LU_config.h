/*	This software package exactly solves a dense system of linear equations using either a REF_LU 
	factorization, Crout LU factorization, or Doolittle LU factorization. Accompanies the paper:
	
	"Efficient Validation of Basic Solutions via Integer-Preserving Factorization"

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Authors--------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

		Adolfo R. Escobedo 
	Erick Moreno-Centeno Christopher Lourenco
	
****** Contact info: adRes@asu.edu *******************************************************************
	
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Copyright------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

	This software is under the GNU General Public License
	See license.txt for license info.
	This software is Copyright by Adolfo Escobedo, Email: adRes@asu.edu

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------DISCLAIMER-----------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

** This software comes with no implied warranty, use it at your own risk *****************************
** ONLY for dense matrices. Sparse matrices may lead to undefined behavior/errors with code **********

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------C++ Libraries--------------------------------------------------------------
------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------*/

# include <stdlib.h>	// Allows the use of many standard C++ functions (atof, rand, etc)
# include <iostream>	// Allows for cout, cerr etc
# include <fstream>	// Allows for file operations (open, close, etc)
# include <chrono>	// Allows for use of steady clock class, which is designed for time intervals
# include <math.h>	// Allows many mathematical operations (ceil, etc)
# include <cmath>	// Allows many matheamtical operations
# include <sstream>	// Allows for input/output from files
# include <gmp.h>	// Allows use of GMP
# include <gmpxx.h>	// Allows use of GMP
# include <mpfr.h>	// Allows use of MPFR library
# include <random>	// For randomly generated matrices
# include <iomanip> 	// Setting precision on cout

using namespace std;
using namespace std::chrono;

/*----------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Default Parameters---------------------------------------------------------
------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------*/

#define DEFAULT_REF 1
#define DEFAULT_Crout 0
#define DEFAULT_Doolittle 0
#define DEFAULT_numRHS 50
#define DEFAULT_dimension 100
#define DEFAULT_seed1 11
#define DEFAULT_seed2 1400
#define DEFAULT_lower -99
#define DEFAULT_upper 99
#define DEFAULT_density 1.0
#define DEFAULT_print 0
#define DEFAULT_check 0
#define DEFAULT_name "default_output.out"
#define DEFAULT_print2 1
#define DEFAULT_types 1
#define DEFAULT_prec 16
#define DEFAULT_out "default_x.out"