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
# define __FUNCT__ "SVL_3D_FFTW"
PetscErrorCode SVL_3D_FFTW(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U, PetscScalar* AC) {
/* Input and Output Variables */ 
   PetscInt         N, i, j, k;
   Vec              x, y;
   //PetscScalar      a;
    PetscBool        view = PETSC_FALSE; // PETSC_TRUE allo to see result
   PetscErrorCode   ierr;

   PetscFunctionBeginUser;

/* U is 3D array, The 3D can be split in 2D slide planes along and axis */
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 2.1: FORWARD FFTW\n");CHKERRQ(ierr);

/* Arrays size */
   N = Nx * Ny * Nz;

/* FFTW forward plan */ 
   fftw_plan     f_plan;
   fftw_complex  *U_in, *U_out;

/* Set fftw Input and Output array pointers */ 
    U_in   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // The input is Real
    U_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // The output is Complex 

/* VecCreateSeqWithArray: Creates a standard,sequential array-style vector, where the user provides the array space to store the vector values.
	#include "petscvec.h" 
	PetscErrorCode  VecCreateSeqWithArray(MPI_Comm comm, PetscInt bs, PetscInt n, const PetscScalar array[], Vec *V)

	Collective on MPI_Comm
	Input Parameter
		comm	- the communicator, should be PETSC_COMM_SELF
		bs	- the block size
		n	- the vector length
		array	- memory where the vector elements are to be stored.
	Output Parameter
		*V        -the vector 
  PetscObjectSetName: Sets a string name associated with a PETSc object.
*/

/* Create vector x with size N.  Set x -> U_in */
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) N,(const PetscScalar*) U_in, &x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Real Space vector");CHKERRQ(ierr);

/* Create vector x with size N. Set y -> U_out */
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 1,(PetscInt) 2 * N,(const PetscScalar*) U_out, &y);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) y, "Frequency space vector");CHKERRQ(ierr);

/* Set up fftw plan  */
    //f_plan = fftw_plan_dft_3d(Nx, Ny, Nz, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);
    PetscInt   dim[3];
    dim[0] = Nx;
    dim[1] = Ny;
    dim[2] = Nx;
    f_plan = fftw_plan_dft(3, dim, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);
    
/* Initialize Real space vector(Input data) : The data in the in/out arrays is overwritten during FFTW_MEASURE planning, so planning should be done before the input is initialized by the user. */
   PetscScalar  *x_array;
   ierr = VecGetArray(x, &x_array);CHKERRQ(ierr);

/* x_array -> U */
   for (i = 0; i < Nx; i++) {
       for (j = 0; j < Ny; j++) {
           for (k = 0; k < Nz; k++) { 
	     x_array[i * Ny * Nz + j * Nz + k] = U[i * Ny * Nz + j * Nz + k]; 
           }
       } 
   }

/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray().  x -> x_array -> U */
   ierr = VecRestoreArray(x, &x_array);CHKERRQ(ierr);

/*view = PETSC_TRUE allow to see x values */
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
/* VecGetArray: Returns a pointer to a contiguous array that contains this processor's portion of the vector data. For the standard PETSc vectors, VecGetArray() returns a pointer to the local data array and does not use any copies. If the underlying vector data is not stored in a contiguous array this routine will copy the data to a contiguous array and return a pointer to that. You MUST call VecRestoreArray() when you no longer need access to the array.

	#include "petscvec.h"   
	PetscErrorCode VecGetArray(Vec x,PetscScalar **a)

	Logically Collective on Vec
	Input Parameter
		x  -the vector 
	Output Parameter
		a  -location to put pointer to the array 

*/

/* View real input and fftw complex output, *xa, *ya are pointer to  vectors x, y */
   PetscScalar  *xa, *ya;
   ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
   ierr = VecGetArray(y,&ya);CHKERRQ(ierr);

// Output data forward fftw  
   for (i = 0; i < Nx; i++) {
       for (j = 0; j < Ny; j++) {
           for (k = 0; k < Nz; k++) { 
               AC[i * Ny * Nz +  j * Nz + k] = ya[i * Ny * Nz +  j * Nz + k]; 
           }
       }
   }

/* VecRestoreArray() : Copy the data back into the underlying vector data structure from the array obtained with VecGetArray(). */
    ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
    ierr = VecRestoreArray(y,&ya);CHKERRQ(ierr);

/* VecDestroy: Destroys a vector.*/
   ierr = VecDestroy(&x);CHKERRQ(ierr);
   ierr = VecDestroy(&y);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} /* END FUNCTION SL_FFTW */

