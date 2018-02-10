This code is used to solve a dense system of linear equations exactly using the REF LU Factorization:

*****************************************************************Typing make in this folder creates three executables*****************************************************************

1) REF_LU.exe: This executable is used to randomly generate a matrix (A) and RHS vector(s) (b) of desired size and solve Ax=b exactly using REF LU, Doolittle LU, and/or Crout LU

2) REF_LU_v2.exe: This executable reads in a matrix and rhs vector that it then solves exactly Ax=b using REF LU, Doolittle LU, and/or Crout LU

3) Q_mat.exe: This executable randomly generates a matrix and RHS vector(s) of desired size and solves Ax=b exactly using REF_LU and Edmond's Q Matrices

A description of each executable is given below:

1) REF_LU.exe:
usage: ./REF_LU.exe

followed by the following options:

        -h,--help
        -a,--algs       [123,13,2,..]   >default=1<
                1:REF_LU, 2:QD_LU, 3:QC_LU
        (Can specify multiple w/o spaces in ascending order)

        -b,--num_RHS      [int]           >default=50<            (Number of RHS vectors to solve w/factorizations)
        -c,--check_sol    [1,0]           >default=0<             (Verify if subsitution outputs are correct)
        -d,--density      [double]        >default=1.0<	
        -f,--file_out     [string]        >default_output.out<    (File to save particular runs)
        -n,--num_vars     [int]           >default=100<		
        -p,--print        [1,0]           >default=0<             (Prints working matrix after each iteration)
        -r,--range        [int] [int]     >default=-99,99<
        -p2 --print2      [1, 0]          >default=1<             (Prints x vectors in a file)
        -t --type         [1,2,3]         >default=1<             (how x is printed. 1: Rational, 2: double, 3: MPFR)
        -o --out          [string]        >default_x.out<         (Filename for x vectors)
        -prec --precision [int]         >default=16<              (MPFR precision)
        -s,--seeds        [int] [int]     >default=11,1400<       (1st seed for nonzero-entries placement, 
                                                                  second for entry values)
							

2) REF_LU_v2.exe:
usage: ./REF_LU_v2.exe MATRIX_FILE_NAME RHS_VECTOR_NAME 
        followed by the flags: -a/--algs, -t/--type, -o/--out, -prec/--precision as defined above
	where:
		MATRIX_FILE_NAME: The direction to a file containing a fully dense matrix stored in matrix market format. Information about matrix market format can be found at: 
		http://math.nist.gov/MatrixMarket/formats.html
		
		RHS_VECTOR_NAME: Direction a file containing a SINGLE rhs vector
				
	****NOTE: DO NOT use this code for sparse matrices, it may lead to undefined behavior and definitely will lead to insufficient memory management and very slow run time****
		
3) ./Q_mat.exe 
usage: ./Q_mat.exe

followed by the following options:

-h,--help
        
        -b,--num_RHS    [int]           >default=50<		(Number of RHS vectors to solve w/factorizations)
        -c,--check_sol  [1,0]           >default=0<		(Verify if subsitution outputs are correct)
        -d,--density    [double]        >default=1.0<			
        -f,--file_out   [string]        >default_output.out<	(File to save particular runs)
        -n,--num_vars   [int]           >default=100<		
        -p,--print      [1,0]           >default=0<		(Prints working matrix after each iteration)
        -r,--range      [int] [int]     >default=-99,99<
        -s,--seeds      [int] [int]     >default=11,1400<	(1st seed for nonzero-entries placement, 
								second for entry values)
