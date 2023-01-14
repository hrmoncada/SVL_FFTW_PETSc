/***************************************************/
/* SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY  */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_FFTW_SWAP"
PetscErrorCode SVL_3D_FFTW_SWAP(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U) {

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 4: SET FFTW & SWAP \n");CHKERRQ(ierr);

/*************************************************************/
/*                   FORWARD FFTW & SWAP                     */
/*             FAST FOURIER TRANSFORM ON THE WEST            */
/*************************************************************/

/* FFTW & SPAW along x-axis */
   ierr = SVL_1D_FFTW_X(Nx, Ny, Nz, U);CHKERRQ(ierr);
   ierr = SVL_1D_SAWP_X(Nx, Ny, Nz, U);CHKERRQ(ierr);

/* FFTW & SPAW along y-axis */
   ierr = SVL_1D_FFTW_Y(Nx, Ny, Nz, U);CHKERRQ(ierr);
   ierr = SVL_1D_SAWP_Y(Nx, Ny, Nz, U);CHKERRQ(ierr);

/* FFTW& SPAW  along z-axis */
   ierr = SVL_1D_FFTW_Z(Nx, Ny, Nz, U);CHKERRQ(ierr);
   ierr = SVL_1D_SAWP_Z(Nx, Ny, Nz, U);CHKERRQ(ierr);

/************** END FORWARD FFTW & SWAP  *********************/ 

  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
