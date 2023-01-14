/***************************************************/
/* SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY  */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_2D_FFTW_SWAP_METHOD_3"
PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_3(PetscInt Nx, PetscInt Ny, PetscReal* U, PetscScalar* AC) {

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 3: SET FFTW & SWAP \n");CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 3 : FORWARD FFTW & SWAP               */
/*             1. 1D FFT on each Vertical Column             */
/*             2. Transpose 2D Matrix                        */   
/*             3. 1D FFT on each Vertical Column             */
/*************************************************************/
  PetscReal   U_array[Ny];
  PetscScalar AC_array[Ny];
  PetscInt    i, j;

/*************************************************************/
/*   Ny x 1D Vertical columns arrangment - FFTW along x-axis */
/*************************************************************/
  for (i = 0; i < Nx; i++) { // Loop through the rows.
// Build a vector of size Nz
      for (j = 0; j < Ny; j++) { // Loop through the columns .  
          U_array[j] = U[i + j * Nx] ;   // Vertical columns  arrangment 
      }
// Get FFTW for the vector         
      ierr = SVL_1D_FFTW(Nx, Ny, U_array, AC_array);CHKERRQ(ierr);      

// Store the results in a new vector array.
      for (j = 0; j < Nx; j++) { // Loop through the rows.
          AC[i + j * Nx] = AC_array[j];
      }
  }

/*************************************************************/
/*                          Transpose                        */
/*************************************************************/
  ierr = SVL_2D_TRANSPOSE_COMPLEX(Nx, Ny, AC);CHKERRQ(ierr); 

/*************************************************************/
/*   Ny x 1D Vertical columns arrangment - FFTW along y-axis */
/*************************************************************/
  for (i = 0; i < Nx; i++) { // Loop through the rows.
// Build a vector of size Nz
      for (j = 0; j < Ny; j++) { // Loop through the column 
          U_array[j] = AC[i + j * Nx] ;   // Vertical columns  arrangment 
      }
// Get FFTW for the vector         
      ierr = SVL_1D_FFTW(Nx, Ny, U_array, AC_array);CHKERRQ(ierr);      

// Store the results in a new vector array.
      for (j = 0; j < Nx; j++) { // Loop through the rows.
          AC[i + j * Nx] = AC_array[j];
      }
  }

/*************************************************************/
/*             For a total of (Nx x Ny) x 1D FFTW            */
/*************************************************************/
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL", "OUTPUT_FFTW_IMAG", Nx, Ny, AC);CHKERRQ(ierr);

/*************************************************************/
/*                      BACKWARD FFTW                        */
/*             FAST FOURIER TRANSFORM ON THE WEST            */
/*************************************************************/
   PetscScalar Inv_U[Nx * Ny];
   ierr = SVL_2D_IFFTW(Nx, Ny, AC, Inv_U);CHKERRQ(ierr);
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_INV_FFTW_REAL", "OUTPUT_INV_FFTW_IMAG" ,Nx, Ny, Inv_U);CHKERRQ(ierr);

/*************************************************************/
/*                           SWAP                            */
/*************************************************************/
/*                            2D                             */
/*                          1    2                           */ 
/*                          4    3                           */
/*      SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY       */
/*                          3    4                           */ 
/*                          2    1                           */
/*************************************************************/
   ierr = SVL_2D_SWAP_QUADRANTS(Nx, Ny, AC);CHKERRQ(ierr);
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_SWAP_REAL","OUTPUT_SWAP_IMAG", Nx, Ny, AC);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
