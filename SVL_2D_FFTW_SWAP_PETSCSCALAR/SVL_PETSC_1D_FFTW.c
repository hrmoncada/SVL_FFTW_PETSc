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
# include "SVL_PETSC_Header.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_1D_FFTW"
PetscErrorCode SVL_1D_FFTW(PetscInt Nx, PetscInt Ny, PetscReal* U, PetscScalar* AC) {
/* Input and Output Variables */ 
   //PetscInt         i, j; 
   PetscInt         j, N;
   Vec              x, y;
   //PetscScalar      a;
   PetscBool        view = PETSC_FALSE; // PETSC_TRUE allo to see result
   PetscErrorCode   ierr;

   PetscFunctionBeginUser;

   ierr = PetscPrintf(PETSC_COMM_WORLD,"\n	STEP 2.1: FORWARD FFTW\n");CHKERRQ(ierr);

/* arrays size */
   N = Ny;

/* FFTW forward plan */ 
    fftw_plan     f_plan;
    fftw_complex  *U_in, *U_out;

/* Set fftw Input and Output array pointers */ 
    U_in   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    U_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

/* Create vector x with size N.  Set x -> U_in */
   ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_in, &x);CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) x, "Real Space vector");CHKERRQ(ierr);

/* Create vector x with size N. Set y -> U_out */
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_out, &y);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) y, "Frequency space vector");CHKERRQ(ierr);

/* Set up fftw plan  */
   f_plan = fftw_plan_dft_1d(N, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);
  
   /*PetscInt   dim[1];
   dim[0] = N; // N = Ny
   f_plan = fftw_plan_dft(1, dim, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);
 */
/* Initialize Real space vector(Input data) : The data in the in/out arrays is overwritten during FFTW_MEASURE planning, so planning should be done before the input is initialized by the user. */
   PetscScalar  *x_array;
   ierr = VecGetArray(x, &x_array);CHKERRQ(ierr);

/* x_array -> U */
   for (j = 0; j < Ny; j++) {x_array[j] = U[j]; }
/*   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
          x_array[i * Ny + j] = U[i * Ny + j]; 
       } 
   }
*/

/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray().  x -> x_array -> U */
   ierr = VecRestoreArray(x, &x_array);CHKERRQ(ierr);

/*view = PETSC_TRUE allow to see x values */
   if (view) {ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}

/* Execute forward fftw plan, output file U_out */
    fftw_execute(f_plan);
 
/*************************************************************/
/*                END FUNCTION FORWARD FFTW                  */
/*************************************************************/

/*************************************************************/
/*                           VIEW                            */
/*************************************************************/
/* View real input and fftw complex output*/
    PetscScalar  *xa, *ya;
    ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
    ierr = VecGetArray(y,&ya);CHKERRQ(ierr);

// Output data forward fftw 
   for (j = 0; j < Ny; j++) { AC[j] = ya[j];}
 
/*    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            U[i * Ny + j] =  ya[i * Ny + j];
        }
    }
*/
/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray(). */
    ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
    ierr = VecRestoreArray(y,&ya);CHKERRQ(ierr);

/* Destroy plan*/
    fftw_destroy_plan(f_plan);

/* Free spaces */
    fftw_free(U_in);  ierr = VecDestroy(&x);CHKERRQ(ierr);
    fftw_free(U_out); ierr = VecDestroy(&y);CHKERRQ(ierr);

    PetscFunctionReturn(0); 
} /* END FUNCTION SL_FFTW */

