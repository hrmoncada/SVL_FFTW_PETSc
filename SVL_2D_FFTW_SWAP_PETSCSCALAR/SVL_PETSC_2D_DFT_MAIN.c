/***********************************************************/
/*         SPATIAL VARIANT LATTICE                         */
/***********************************************************/
static char help[] = "Spatially Variant Lattice Project\n";
/*PETSC libraries*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <complex.h>
# include <fftw3.h>
# include <petsc.h>
# include <petscksp.h>
# include <petscmat.h>
# include <petscvec.h>
# include "SVL_PETSC_Header.h"

/*************************************************************/
/*                     START MAIN PROGRAM                    */
/*************************************************************/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {
  PetscErrorCode     ierr;
  PetscMPIInt        size, rank;
  PetscInt           Nx = 256 , Ny  = Nx; //256
  //PetscReal          Lx = 1.0 , Ly  = 1.0;
  PetscInt           NM = 11 , NN  = NM;  // Number of spatial Harmonics
  PetscInt           NK = NM * NN; // Total number of FFTW spatial harmonics NK = = NM * NN for the 2D

/***************************************************/
/*                  START PETSC                    */
/***************************************************/
  ierr  = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  ierr  = MPI_Comm_rank(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
  ierr  = MPI_Comm_size(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"         SPACIAL VARIANT LATTICE       \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");

/* Star-t Elapse time*/
  clock_t tic = clock();

/*************************************************************/
/*                CREAT UNIT MATRIX U                        */
/*************************************************************/  
   //ierr = SVL_ZERO_CELL(Nx, Ny, U);CHKERRQ(ierr);
   //ierr = SAVE_1D_To_2D_ARRAY_Complex("OUTPUT_ZERO_CELL_REAL", "OUTPUT_ZERO_CELL_IMAG" ,Nx, Ny, U);CHKERRQ(ierr);

/*************************************************************/
/*                  ADD TRIANGLE                             */
/* INSERTING AN EQUILATERAL TRIANGLE INTO A SQUARE UNIT CELL */
/*************************************************************/   
   PetscReal   U[Nx * Ny];
   ierr = SVL_2D_TRIANGLE(Nx, Ny, U);CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 1 : FORWARD FFTW & SWAP               */
/*             1. 1D FFT on each Vertical Column             */
/*             2. 1D FFT on each Horizontal Row              */
/*************************************************************/
  PetscScalar AC[Nx * Ny];
  //ierr = SVL_2D_FFTW_SWAP_METHOD_1(Nx, Ny, U, AC);CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 2 : FORWARD FFTW & SWAP               */
/*             1. Transpose 2D Matrix                        */ 
/*             2. 1D FFT on each Horizontal Row              */
/*             3. Transpose 2D Matrix                        */   
/*             4. 1D FFT on each Horizontal Row              */
/*************************************************************/
 ierr = SVL_2D_FFTW_SWAP_METHOD_2(Nx, Ny, U, AC);CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 3 : FORWARD FFTW & SWAP               */
/*             1. 1D FFT on each Vertical Column             */
/*             2. Transpose 2D Matrix                        */   
/*             3. 1D FFT on each Vertical Column             */
/*************************************************************/
 //ierr = SVL_2D_FFTW_SWAP_METHOD_3(Nx, Ny, U, AC);CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 4 : FORWARD FFTW & SWAP               */
/*             1. 2D FFT on the 2D array                     */
/*************************************************************/
 //ierr = SVL_2D_FFTW_SWAP_METHOD_4(Nx, Ny, U, AC);CHKERRQ(ierr);

/*************************************************************/
/*               TRUCATE THE SPATIAL HARMONICS               */
/*************************************************************/
   PetscScalar  AMN[NK]; // TRUNCATED ARRAY OF FFTW SPATIAL HARMONICS
   ierr = SVL_2D_TRUNCATE_FFTW_ARRAY(Nx, Ny, NM, NN, AC, AMN);CHKERRQ(ierr);

   //ierr = SVL_2D_TRUNCATE_FFTW_SPATIAL_HARMONIC(Nx, Ny, NM, NN, U, AMN);CHKERRQ(ierr);


/*************************************************************/
/*                      END MAIN LOOP                       */
/*************************************************************/

/*End Elapsed time */
    clock_t toc = clock(); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n=============THE END===========\n");
    ierr = PetscFinalize();
    return 0;
}





