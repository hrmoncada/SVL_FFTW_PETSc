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
# define __FUNCT__ "SVL_2D_IFFTW"
PetscErrorCode SVL_2D_IFFTW(PetscInt Nx, PetscInt Ny, PetscScalar* U, PetscScalar* Inv_U) {
/* Input and Output*/ 
    PetscInt i, j; 
    fftw_complex    *U_in2, *U_out2; 
    PetscInt         N = Nx * Ny;
    Vec              y2, z; 
    PetscScalar      a;
    PetscBool        view = PETSC_FALSE;
    PetscErrorCode   ierr;
  
    PetscFunctionBeginUser;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n	STEP 2.2: BACKWARD FFTW\n");CHKERRQ(ierr);

/* FFTW backward plan */ 
    fftw_plan  b_plan;

/* Set Input and Output */ 
    U_in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *  N);
    U_out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_in2, &y2);CHKERRQ(ierr); 
    ierr = PetscObjectSetName((PetscObject) y2, "Frequency space vector");CHKERRQ(ierr);

    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_out2, &z);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) z, "Reconstructed vector");CHKERRQ(ierr);

/*Set up plan  */
    b_plan = fftw_plan_dft_2d(Nx, Ny, U_in2, U_out2, FFTW_BACKWARD, FFTW_ESTIMATE);

/* Initialize Real Space Vector (Input data): x_array -> x */
   PetscScalar  *y2_array;
   ierr = VecGetArray(y2, &y2_array);CHKERRQ(ierr);

/* x_array -> U */
   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
          y2_array[i * Ny + j] = U[i * Ny + j]; 
       } 
   }

/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray().  x -> x_array -> U */
   ierr = VecRestoreArray(y2, &y2_array);CHKERRQ(ierr);

   if (view){ierr = VecView(y2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}

/* Execute backward fftw plan, output file U_out2 */
    fftw_execute(b_plan);

/* FFTW computes an unnormalized DFT, thus z = N*x */
    a = 1.0/(PetscReal) N;
    ierr = VecScale(z,a);CHKERRQ(ierr);

/*************************************************************/
/*                END FUNCTION BACKWARD FFTW                 */
/*************************************************************/
/*************************************************************/
/*                           VIEW                            */
/*************************************************************/
/* view complex output*/
   PetscScalar  *ya2, *za; 
   ierr = VecGetArray(y2,&ya2);CHKERRQ(ierr);
   ierr = VecGetArray(z,&za);CHKERRQ(ierr);

// output data forward fftw  
  for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
          Inv_U[i * Ny + j] =  za[i * Ny + j]; 
      }
  }
/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray(). */
  ierr = VecRestoreArray(y2,&ya2);CHKERRQ(ierr);
  ierr = VecRestoreArray(z,&za);CHKERRQ(ierr);

/* Destroy plan*/
  fftw_destroy_plan(b_plan);

/* Free spaces */
  fftw_free(U_in2);     ierr = VecDestroy(&y2);CHKERRQ(ierr);
  fftw_free(U_out2);    ierr = VecDestroy(&z);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} /* END FUNCTION SL_FFTW */

