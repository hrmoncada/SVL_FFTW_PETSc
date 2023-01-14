/*************************************************************/
/*                   CYLINDRICAL TRANSLATION                 */
/*************************************************************/
/*************************************************************/
/*                  CARTESIAN TO CYLINDRICAL                 */
/*                            AND                            */
/*                  CYLINDRICAL TO CARTESIAN                 */
/*************************************************************/ 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_CYLINDRICAL_TRANSLATION"
PetscErrorCode SVL_3D_CYLINDRICAL_TRANSLATION(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscReal* Kx, PetscReal* Ky, PetscReal* Kz, PetscReal* THETA, PetscReal* PER) {

  PetscReal       TH[New_Nx * New_Ny * New_Nz],  RHO[New_Nx * New_Ny * New_Nz];
  PetscErrorCode ierr;
 
  PetscFunctionBeginUser;

/*************************************************************/
/*   Attribute 1 - Lattice Orientation Function              */
/*   K field is oriented with angles THETA(r) & VARPHI(r)    */ 
/*************************************************************/  
/* TRANSFORMS CARTESIAN COORDINATES TO POLAR COORDINATES, output: TH and RHO */
      ierr = SVL_3D_CART2CYLIN(New_Nx, New_Ny, New_Nz, Kx, Ky, Kz, RHO, TH); CHKERRQ(ierr); // (x, y, z) -> (r, thera, phi)

/* PERFORME ANGLE ROTATION */
      ierr = SVL_3D_CYLIN_ROTATION(New_Nx, New_Ny, New_Nz, TH, THETA); CHKERRQ(ierr);

/***************************************************************/
/* Attribute 2 - Lattice Spacing Function  RHO = RHO/PER       */ 
/***************************************************************/
      ierr = SVL_3D_CYLIN_SPACING(New_Nx, New_Ny, New_Nz, RHO, PER);CHKERRQ(ierr);

/* TRANSFORMS THE POLAR COORDINATE TO 2D CARTESIAN COORDINATES, output: Kx and Ky */    
      ierr = SVL_3D_CYLIN2CART(New_Nx, New_Ny, New_Nz, RHO, TH, Kx, Ky, Kz);CHKERRQ(ierr); //  (r, thera, phi) -> (x, y, z)

  PetscFunctionReturn(0);
} //END FUNCTION ELIMINATION_GRATING_ACCORD_AMPLITUD
