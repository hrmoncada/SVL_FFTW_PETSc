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
# define __FUNCT__ "SVL_2D_SWAP_QUADRANTS"
PetscErrorCode SVL_2D_SWAP_QUADRANTS(PetscInt Nx, PetscInt Ny, PetscScalar* U) {
  PetscInt i, j;
  PetscScalar tmp13, tmp24;
  PetscInt N2x = Nx/2;
  PetscInt N2y = Ny/2;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n	STEP 2.3: SET FFTW_SWAP\n");CHKERRQ(ierr);

  PetscScalar U_SWAP[Nx * Ny]; // SWAP ARRAY

  for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
           U_SWAP[i * Ny + j] =  U[i * Ny + j]; 
       }
   }


for(i = 0; i < N2x; ++i) {
     for(j = 0; j < N2y; ++j) {
	tmp13 = U_SWAP[i * Ny + j];
	U_SWAP[i * Ny + j] = U_SWAP[(i + N2x)*Ny + (j + N2y)] ;
	U_SWAP[(i + N2x)*Ny + (j + N2y)] = tmp13;

	tmp24 = U_SWAP[(i * Ny + N2x) + j];
	U_SWAP[(i * Ny + N2x) + j] = U_SWAP[i * Ny + j + N2y*Ny];
	U_SWAP[i * Ny + j + N2y*Ny] = tmp24;
      }
  }


  for (i = 0; i < Nx; i++) {
       for (j = 0; j < Ny; j++) {
           U[i * Ny + j] = U_SWAP[i * Ny + j];
       }	
  }

  PetscFunctionReturn(0);
}/* END FUNCTION SVL_SWAPQUADRANTS */
