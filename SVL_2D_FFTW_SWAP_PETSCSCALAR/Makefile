# The CFLAGS variable sets compile flags for gcc:
#  -o        Filename to specify the name of the executable to be something other than the default 'a.out'
#  -g        Compile with debug information
#  -lm       To link in the math library (it contains functions like sqrt, cos, sin, etc)
#  -lfftw3   To link to the fftw3 library. 
#  -Wall     Give verbose compiler warnings, turn on all compiler warnings (this is usually a good idea):
#  -O0       Do not optimize generated code.This is the default.
#  -O1       Optimize.Optimizing compilation takes somewhat more time, and a lot more memory for a large function. 
#  -O2       Optimize even more. GCC performs nearly all supported optimizations that do not involve a space-speed tradeoff.
#            The compiler does not perform loop unrolling or function inlining when you specify '-O2'. 
#            As compared to '-O', this option increases both compilation time and the performance of the generated code. 
#  -O3       Optimize yet more. '-O3' turns on all optimizations specified by '-O2' and also turns on the '-finline-functions' and '-frename-registers' options.
#  -Os       Optimize for size. '-Os' enables all '-O2' optimizations that do not typically increase code size. It also performs further optimizations designed to reduce code size.
#  -std=c99  Use the C99 standard language definition, ISO C99. Other option are C89, gnu89, gnu99
# They are add on the follow way:
# CFLAGS = -g -Wall -O0 -std=c99
# LIBS = -lm, -lfftw3 

# Compile:
#  >> export PETSC_DIR=~/Desktop/PETSC/petsc-3.5.4
#  >> export PETSC_ARCH=linux-gnu-complex
# Compile:
#  >> make SVL_PETSC
# Execute:
#  >> mpirun -np 4 out
# All previos step can be done with SVL.sh 
# >> bash SVL.sh
CFLAGS	         = 
FFLAGS	         = 
CPPFLAGS         = 
FPPFLAGS         =
LOCDIR           =  					# Working folder
EXAMPLESC        =  					# *.c or *.cpp file names here
EXAMPLESF        =
MANSEC           = 
LIBS             = -lfftw3 -lm 		        	# serial fftw libraries
MPI_LIBS         = -lfftw3_mpi -lfftw3 -lm 		# parallel fftw libraries	
TRASH	         :=  *.*~  *~  *.o 

# SVL_PETSC_DFT_01.c SVL_PETSC_DFT_02.c SVL_PETSC_DFT_03.c
# Note: PETSC main -> SVL_PETSC_DFT_MAIN.c

SOURCE_PETSC     := SVL_PETSC_2D_DFT_MAIN.c\
		    SVL_PETSC_2D_TRIANGLE_CELL.c\
                    SVL_PETSC_2D_FFTW_SWAP_METHOD_1.c\
                    SVL_PETSC_2D_FFTW_SWAP_METHOD_2.c\
                    SVL_PETSC_2D_FFTW_SWAP_METHOD_3.c\
                    SVL_PETSC_2D_FFTW_SWAP_METHOD_4.c\
		    SVL_PETSC_1D_FFTW.c\
		    SVL_PETSC_2D_FFTW.c\
		    SVL_PETSC_2D_TRANSPOSE_REAL.c\
		    SVL_PETSC_2D_TRANSPOSE_COMPLEX.c\
		    SVL_PETSC_2D_IFFTW.c\
		    SVL_PETSC_2D_SWAP_QUADRANTS.c\
		    SVL_PETSC_2D_TRUNCATE_FFTW_ARRAY.c\
		    SAVE_1D_To_2D_ARRAY_REAL.c\
		    SAVE_1D_To_2D_ARRAY_COMPLEX.c\
		    #SVL_PETSC_2D_TRANSPOSE.c\
		    #SVL_PETSC_2D_ZERO_CELL.c\
		    #SVL_PETSC_2D_UNIT_CELL.c\
		    #SVL_PETSC_2D_GRATING_AND_IMPROVEMENTS.c\
		    #SVL_PETSC_2D_GRADING_VECTOR.c\
		    #SVL_PETSC_2D_ELIMINATED_GRATING_ACCORD_THEIR_AMPLITUD.c\
		    #SVL_PETSC_2D_IDENTIFIED_COLLINEAR_PLANAR_GRATING.c\
                    #SVL_PETSC_2D_IMPLEMENT_IMPROVEMENTS.c\
		    #SVL_PETSC_1D_FFTW_X_AXIS.c\
                    #SVL_PETSC_1D_FFTW_Y_AXIS.c\
		    #SVL_PETSC_2D_Y_ROW_ARRANGE.c\
		    #SVL_PETSC_1D_SWAP_X_AXIS.c\
		    #SVL_PETSC_1D_SWAP_Y_AXIS.c\
		    #SVL_PETSC_2D_TRUNCATE_FFTW_SPATIAL_HARMONICS.c\
		    #SVL_PETSC_Orientation_Function.c\
		    #SVL_PETSC_fdder.c\
		    #SVL_PETSC_LOOP.c\
		    #SVL_PETSC_Orientation_Vector.c\
		    #SVL_PETSC_Cartesian_To_Polar.c\
		    #SVL_PETSC_Rotation.c\
		    #SVL_PETSC_Polar_To_Cartesian.c\
                    #SVL_PETSC_RHS.c\
		    #SVL_PETSC_Print_Real.c\
		    #SVL_PETSC_Print_Complex.c\
                    #SVL_PETSC_Lattice_Spacing_Function.c\
	            #SVL_PETSC_Spacing.c\

OBJECTS_PETSC    := $(SOURCE_PETSC:.c=.o)
EXECUTABLE_PETSC := OUTPUT_PETSC

# Version 3.5.4: These Makefiles lines must be updated every time you Update your petsc version.
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# Version 3.7.6: These Makefiles lines must be updated every time you Update your petsc version.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

SVL_PETSC: $: $(OBJECTS_PETSC) chkopts
	-${CLINKER} -o $(EXECUTABLE_PETSC) $(OBJECTS_PETSC) ${PETSC_LIB} ${LIBS}
	 ${RM} $(OBJECTS_PETSC) $(TRASH)
	 #mpirun -np 4 ./OUTPUT_PETSC -ksp_type cg -ksp_converged_reason > log.txt
	 mpirun -np 4 ./OUTPUT_PETSC > log.txt


