/*************************************************************/
/*                     IMPLEMENT_ATTRIBUTES                  */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_CYLINDRIC_SPATIAL_VARIANT"
PetscErrorCode SVL_3D_CYLINDRIC_SPATIAL_VARIANT(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscReal New_dx, PetscReal New_dy, PetscReal New_dz, PetscReal* X, PetscReal* Y, PetscReal* Z, PetscReal* RSQ, PetscReal* THETA, PetscReal* ZZ, PetscReal* PER) {

  PetscErrorCode ierr;
 
  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 5: IMPLEMENT ATTRIBUTES\n");CHKERRQ(ierr);

/*************************************************************/
/* Attribute 1 - ORIENTATION FUNCTION THETA(r) & ZH(r)       */
/*************************************************************/ 
  ierr = SVL_3D_CYLINDRICAL_ORIENTATION_FUNCTION(New_Nx, New_Ny, New_Nz, New_dx, New_dy, New_dz, X, Y, Z, RSQ, THETA); CHKERRQ(ierr);

/*************************************************************/
/* Attribute 2 - LATTICE SPACING FUNCTION  PER(r)            */
/*************************************************************/ 
  ierr = SVL_3D_CYLINDRICAL_LATTICE_SPACING_FUNCTION(New_Nx, New_Ny, New_Nz, RSQ, PER);CHKERRQ(ierr);

  PetscFunctionReturn(0);

} //END FUNCTION ELIMINATION_GRATING_ACCORD_AMPLITUD
