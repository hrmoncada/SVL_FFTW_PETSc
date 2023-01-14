# Compile:
#  >> export PETSC_DIR=~/Desktop/PETSC/petsc-3.5.4
#  >> export PETSC_ARCH=linux-gnu-complex
# Compile:
#  >> make SVL_PETSC
# Execute:
#  >> mpirun -np 4 out
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

SOURCE_PETSC     := SVL_PETSC_3D_DFT_MAIN_ORNL.c\
		SVL_PETSC_3D_UNIT_CELL.c\
		SVL_PETSC_3D_FFTW_SWAP_METHOD_1.c\
                    SVL_PETSC_3D_FFTW_SWAP_METHOD_2.c\
		SVL_PETSC_1D_FFTW.c\
		SVL_PETSC_3D_FFTW.c\
		SVL_PETSC_3D_SWAP_QUADRANTS.c\
		SVL_PETSC_3D_TRANSPOSE_1_COMPLEX.c\
		SVL_PETSC_3D_TRANSPOSE_2_COMPLEX.c\
		SVL_PETSC_3D_TRANSPOSE_3_COMPLEX.c\
		SVL_PETSC_3D_TRANSPOSE_4_COMPLEX.c\
		SVL_PETSC_3D_IFFTW.c\
		SVL_PETSC_3D_TRUNCATED_FFTW_ARRAY.c\
                    SAVE_1D_TO_3D_ARRAY_REAL.c\
		SAVE_1D_TO_3D_ARRAY_COMPLEX.c\
 		#SVL_PETSC_3D_GRADING_VECTOR.c\
                    #SVL_PETSC_3D_IMPLEMENT_IMPROVEMENTS.c\
		#SVL_PETSC_3D_ELIMINATE_GRATING_ACCORD_THEIR_AMPLITUD.c\
		#SVL_PETSC_3D_IDENTIFIED_COLLINEAR_PLANAR_GRATING.c\
                    #SVL_PETSC_3D_SPHERICAL_SPATIAL_VARIANT.c\
                    #SVL_PETSC_3D_SPHERICAL_ORIENTATION_FUNCTION.c\
                    #SVL_PETSC_3D_SPHERICAL_LATTICE_SPACING_FUNCTION.c\
                    #SVL_PETSC_3D_CYLINDRICAL_SPATIAL_VARIANT.c\
                    #SVL_PETSC_3D_CYLINDRICAL_ORIENTATION_FUNCTION.c\
                    #SVL_PETSC_3D_CYLINDRICAL_LATTICE_SPACING_FUNCTION.c\
                    #SVL_PETSC_3D_IMPROVEMENTS.c\
		#SVL_PETSC_3D_FDDER.c\
		#SVL_PETSC_3D_LOOP.c\
		#SVL_PETSC_3D_ORIENTATION_VECTOR.c\
                    #SVL_PETSC_3D_SPHERICAL_TRANSLATION.c\
		#SVL_PETSC_3D_CARTESIAN_TO_SPHERICAL.c\
		#SVL_PETSC_3D_SPHERICAL_ROTATION.c\
	          #SVL_PETSC_3D_SPHERICAL_SPACING.c\
		#SVL_PETSC_3D_SPHERICAL_TO_CARTESIAN.c\
                    #SVL_PETSC_3D_CYLINDRICAL_TRANSLATION.c\
		#SVL_PETSC_3D_CARTESIAN_TO_CYLINDRICAL.c\
		#SVL_PETSC_3D_CYLINDRICAL_ROTATION.c\
	          #SVL_PETSC_3D_CYLINDRICAL_SPACING.c\
		#SVL_PETSC_3D_CYLINDRICAL_TO_CARTESIAN.c\
                    #SVL_PETSC_3D_RHS.c\

                                        
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
	 mpirun -np 4 ./OUTPUT_PETSC -ksp_type cg -ksp_converged_reason > log.txt
	 #mpirun -np 4 ./OUTPUT_PETSC > log.txt

