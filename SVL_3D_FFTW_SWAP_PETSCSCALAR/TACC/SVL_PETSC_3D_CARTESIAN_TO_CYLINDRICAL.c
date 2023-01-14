/*************************************************************/
/*                   CARTESIAN TO POLAR                      */
/*************************************************************/ 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_CART2CYLIN"
PetscErrorCode SVL_3D_CART2CYLIN(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscReal* Kx, PetscReal* Ky, PetscReal* Kz, PetscReal* RHO, PetscReal* TH){
  PetscInt i, j, k;
  //PetscErrorCode     ierr;

  PetscFunctionBeginUser;

  for (i = 0; i < New_Nx; i++) {
      for (j = 0; j < New_Ny; j++) {
          for (k = 0; k < New_Nz; k++) { 
                RHO[i * New_Nx * New_Ny + j * New_Nz + k] = sqrt(Kx[i * New_Nx * New_Ny + j * New_Nz + k] * Kx[i * New_Nx * New_Ny + j * New_Nz + k]
                                                          + Ky[i * New_Nx * New_Ny + j * New_Nz + k] * Ky[i * New_Nx * New_Ny + j * New_Nz + k]);

                TH[i * New_Nx * New_Ny + j * New_Nz + k]  = atan2(Ky[i * New_Nx * New_Ny + j * New_Nz + k], Kx[i * New_Nx * New_Ny + j * New_Nz + k]);
            }
       }    
  }
  PetscFunctionReturn(0); 
 }  /*   END CART2POL     */
