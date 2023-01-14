/***************************************************/
/*     TRUCATE THE SPATIAL HARMONICS               */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_2D_TRUNCATE_FFTW_ARRAY"
PetscErrorCode SVL_2D_TRUNCATE_FFTW_ARRAY(PetscInt Nx, PetscInt Ny, PetscInt NM, PetscInt NN, PetscScalar* U, PetscScalar* AMN) {
  PetscInt i, j;
  PetscInt m0, n0, n1, n2, m1, m2;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 5: SET TRUNCATION ARRAY OF FFTW SPATIAL HARMONICS\n");CHKERRQ(ierr);

  //double complex *TA  =  malloc(NM * NN * sizeof(double complex)); 
  PetscScalar  TA[NM*NN];
/*
Loop in C: For i = 0 ; i < 16  (Start in 0 until 15),  for example
       A[0] A[1] A[2] A[3] A[4] A[5] A[6] A[7] A[8] A[9] A[10] A[11] A[12] A[13] A[14] A[15]
                       m1                       m0                           m2
                       n1                       n0                           n2
*/

                     // floor: return the last integral value less than o equal to x
  m0 = floor(Nx/2);  // x-axis middle point, ex: floor of 12.4 is 12.0
  n0 = floor(Ny/2);  // y-axis middle point, ex: floor of 12.4 is 12.0

  m1 = m0 - floor(NM/2); // ex: floor of 12.6 is 12.0
  m2 = m0 + floor(NM/2); // ex: floor of 13.1 is 13.0
  n1 = n0 - floor(NN/2); // ex: floor of -2.3 is -3.0
  n2 = n0 + floor(NN/2); // ex: floor of -3.8 is -4.0
   
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n         Truncated Harmonic Values         \n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"(m0, m1, m2, n0, n1, n2) = (%d, %d, %d, %d, %d, %d)\n", m0, m1, m2, n0, n1, n2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n         Vector Column AMN(:) \n");CHKERRQ(ierr);
 
/* Cutting and translating the Highest Harmonics Array from the Harmonics array*/ 
  for (i = m1 ; i <= m2; ++i) {
      for (j = n1 ; j <= n2; ++j) {
           TA[(i-m1)*(n2 - n1 +1) + (j-n1)] =  U[i * Ny + j];
       }	
   }

/* Spatial Harmonics - Rewritting TA as a Column Array, AMN = TA(:) */
   for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           AMN[i * NN + j]  =   TA[j * NN + i]; 
      }	
   } 

/* Save Truncate Array */   
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_TRUNC_FFTW_REAL", "OUTPUT_TRUNC_FFTW_IMAG", NM, NN, AMN);CHKERRQ(ierr);

   PetscFunctionReturn(0); 
} //END FUNCTION TRUNCATE FFTW SPATIAL HARMONIC 
