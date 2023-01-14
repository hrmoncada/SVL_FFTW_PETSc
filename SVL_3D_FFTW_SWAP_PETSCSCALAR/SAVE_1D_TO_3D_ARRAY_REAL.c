/*************************************************************/
/*      INPUT 2D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY       */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

#undef __FUNCT__
#define __FUNCT__ "SAVE_1D_TO_3D_ARRAY_REAL"
PetscErrorCode SAVE_1D_TO_3D_ARRAY_REAL(const char* desc, PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscReal* U)  {
  PetscInt       i, j, k;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

 /* open file and write header */
  FILE *fp; 
  fp = fopen(desc, "w"); // save mat into a file as Square Unit Cell of zeros

  for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) { 
           for (k = 0; k < Nz; ++k) {     
	       ierr = PetscFPrintf(PETSC_COMM_WORLD,fp,"%5.2f   ", U[i * Ny * Nz +  j * Nz + k]);CHKERRQ(ierr);
           }
           ierr = PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");CHKERRQ(ierr);
       }	
       ierr = PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");CHKERRQ(ierr);
  }
  fclose(fp);

  PetscFunctionReturn(0);
}
