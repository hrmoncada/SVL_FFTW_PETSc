/***************************************************/
/*          SVL FILL FRACTION FUNCTION             */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_SPHERICAL_LATTICE_SPACING_FUNCTION"
PetscErrorCode SVL_3D_SPHERICAL_LATTICE_SPACING_FUNCTION(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz,  PetscReal* RSQ, PetscReal* PER) {
  PetscInt  i, j, k;
  PetscReal a = 1.0;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;


/* Find the threshold : Highest harmonics value */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 5.2: LATTICE SPACING FUNCTION \n");CHKERRQ(ierr);

  for (i = 0; i < New_Nx; ++i) {
     for (j = 0; j < New_Ny; ++j)  { 
         for (k = 0; k < New_Nz; ++k)  { 
             //PER[i * New_Ny * New_Nz +  j * New_Nz + k] = RSQ[i * New_Ny * New_Nz + j * New_Nz + k];// Used to check data transfer, it work
             PER[i * New_Ny * New_Nz +  j * New_Nz + k] = 1.0 + exp(-RSQ[i * New_Ny * New_Nz + j * New_Nz + k]/(2*a*a));
         }
     }    
  }

  ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_PER", New_Nx, New_Ny, New_Nz, PER);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} /*    END SVL FILL FRACTION FUNCTION    */
