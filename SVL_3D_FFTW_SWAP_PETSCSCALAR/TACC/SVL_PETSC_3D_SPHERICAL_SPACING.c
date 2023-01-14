/***************************************************/
/*          SVL FILL FRACTION FUNCTION             */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_SPHER_SPACING"
PetscErrorCode SVL_3D_SPHER_SPACING (PetscInt New_Nx, PetscInt New_Ny,  PetscInt New_Nz, PetscReal* RHO, PetscReal* PER) {
  PetscInt i, j, k;

  PetscFunctionBeginUser;
  for (i = 0; i < New_Nx; ++i){
      for (j = 0; j < New_Ny; ++j){
          for (k = 0; k < New_Nz; ++k){
              RHO[i * New_Nx * New_Ny + j * New_Nz + k] = RHO[i * New_Nx * New_Ny + j * New_Nz + k] / PER[i * New_Nx * New_Ny + j * New_Nz + k];  
          }
      }    
  }

   PetscFunctionReturn(0); 
} /*    END SVL FILL FRACTION FUNCTION    */
