#-------------------------------------------------------------------
# Compiler selection and flags
#--------------------------------------------------------------------
CC  = gcc
CCC = g++
debug = -g -o
optimize = -o

#--------------------------------------------------------------------
#Library flags
#--------------------------------------------------------------------
CCLSTD = -std=c++11
CCLGMP =  -lgmp -lgmpxx
CCLMPFR = -lmpfr
#
#

#--------------------------------------------------------------------
# Applications
#--------------------------------------------------------------------- 
APPS = REF_LU REF_LU_v2 Q_mat 
Debug = debug1 debug2 debug3

all: $(APPS)
debug: $(Debug)
	
REF_LU: REF_LU.cpp 
	$(CCC) REF_LU.cpp $(optimize) REF_LU.exe $(CCLMPFR) $(CCLGMP) $(CCLSTD)

REF_LU_v2: REF_LU_v2.cpp 
	$(CCC) REF_LU_v2.cpp $(optimize) REF_LU_v2.exe $(CCLMPFR) $(CCLGMP) $(CCLSTD)
	
Q_mat: Q_mat.cpp 
	$(CCC) Q_mat.cpp $(optimize) Q_mat.exe $(CCLMPFR) $(CCLGMP) $(CCLSTD)

debug1: REF_LU.cpp 
	$(CCC) REF_LU.cpp $(debug) debug_REF.exe $(CCLMPFR) $(CCLGMP) $(CCLSTD)

debug2: REF_LU_v2.cpp 
	$(CCC) REF_LU_v2.cpp $(debug) debug_REF_v2.exe $(CCLMPFR) $(CCLGMP) $(CCLSTD)

debug3: Q_mat.cpp 
	$(CCC) Q_mat.cpp $(debug) debug_Q.exe $(CCLMPFR) $(CCLGMP) $(CCLSTD)

#--------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------
clean:
	rm *.exe
