/***************************************************/
/*                  FFTW                           */
/*     FAST FOURIER TRANSFORM ON THE WEST          */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <fftw3.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_1D_FFTW"
PetscErrorCode SVL_1D_FFTW(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U, PetscScalar* AC) {
/* Input and Output Variables */ 
  // PetscInt         i, j, k;
   PetscInt         k, N; 
   Vec              x, y;
   //PetscScalar      a;
    PetscBool        view = PETSC_FALSE; // PETSC_TRUE allo to see result
   PetscErrorCode   ierr;

   PetscFunctionBeginUser;

/* U is 3D array, The 3D can be split in 2D slide planes along and axis */
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 2.1: FORWARD FFTW\n");CHKERRQ(ierr);

/* arrays size */
   N = Nz;

/* FFTW forward plan */ 
   fftw_plan  f_plan;
   fftw_complex    *U_in, *U_out;

/* Set fftw Input and Output array pointers */ 
   U_in   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // The input is Real
   U_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // The output is Complex 

/* Create vector x with size N.  Set x -> U_in */
   ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_in, &x);CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) x, "Real Space vector");CHKERRQ(ierr);

/* Create vector x with size N. Set y -> U_out */
   ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_out, &y);CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) y, "Frequency space vector");CHKERRQ(ierr);

/* Set up fftw plan  */
   f_plan = fftw_plan_dft_1d(N, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);
  
   /*PetscInt   dim[1];
   dim[0] = N; // N = Nz
   f_plan = fftw_plan_dft(1, dim, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);
 */
/* Initialize Real space vector(Input data) : The data in the in/out arrays is overwritten during FFTW_MEASURE planning, so planning should be done before the input is initialized by the user. */
   PetscScalar  *x_array;
   ierr = VecGetArray(x, &x_array);CHKERRQ(ierr);

/* x_array -> U */
   for (k = 0; k < Nz; k++) { x_array[k] = U[k]; }

/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray().  x -> x_array -> U */
   ierr = VecRestoreArray(x, &x_array);CHKERRQ(ierr);

/* view = PETSC_TRUE allow to see x values */
   if (view){ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}

/* Execute forward fftw plan, output file U_out */
   fftw_execute(f_plan);

/* Destroy plan*/
    fftw_destroy_plan(f_plan);

/* Free spaces */
   fftw_free(U_in);  
   fftw_free(U_out); 

/*************************************************************/
/*                END FUNCTION FORWARD FFTW                  */
/*************************************************************/

/* View real input and fftw complex output*/
   PetscScalar  *xa, *ya;
   ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
   ierr = VecGetArray(y,&ya);CHKERRQ(ierr);

// Output data forward fftw 
   for (k = 0; k < Nz; k++) { AC[k] = ya[k]; }

/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray(). */
    ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
    ierr = VecRestoreArray(y,&ya);CHKERRQ(ierr);

/* Free spaces */
   ierr = VecDestroy(&x);CHKERRQ(ierr);
   ierr = VecDestroy(&y);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} /* END FUNCTION SL_FFTW */

