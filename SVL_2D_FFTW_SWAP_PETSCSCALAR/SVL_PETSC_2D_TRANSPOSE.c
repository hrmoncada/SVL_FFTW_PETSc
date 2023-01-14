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
# define __FUNCT__ "SVL_2D_TRANSPOSE"
PetscErrorCode SVL_2D_TRANSPOSE(PetscInt Nx, PetscInt Ny, PetscScalar* U) {
  PetscInt i, j;  
  PetscScalar Transp[Nx * Ny];
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n	STEP 3.1 : TRANSPOSE\n");CHKERRQ(ierr);

// Example Transpose: (0 1 2 | 3 4 5 | 6 7 8) ==> (0 3 4 | 1 4 7 | 2 5 8)
  for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) { 
          Transp[i + j * Nx] = U[i * Ny + j]; 
       }
   }

  for (i = 0; i < Nx; i++) {
       for (j = 0; j < Ny; j++) {
           U[i * Ny + j] = Transp[i * Ny + j];
       }	
  }

  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
