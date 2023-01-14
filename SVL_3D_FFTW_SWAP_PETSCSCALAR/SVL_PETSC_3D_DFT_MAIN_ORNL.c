/***********************************************************/
/*         SPATIAL VARIANT LATTICE                         */
/***********************************************************/
static char help[] = "Spatially Variant Lattice 3D-Project\n";
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
# include "SVL_PETSC_3D_HEADER.h"

/*************************************************************/
/*              START MAIN PROGRAM Program                   */
/*************************************************************/
struct data  {    
    char vname[50];
    char equal[5];
};
 
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {
  char ch;
  const char *filename = "SVL.data";
  int  i = 0, Array[9];
  struct data input;

/* Input Data */   
  FILE *infile;
  infile = fopen (filename, "r");

  if (infile == NULL){
     fprintf(stderr, "\nError opening file\n");
     exit (1);
  }
/* Read Input Data */
  do {
        if(i < 9){
           fscanf(infile, "%s %s %d\n",input.vname, input.equal, &Array[i]);
           //printf (" %d %s %s %d \n", i, input.vname, input.equal, Array[i]);
        } else if (i == 9) {
           fscanf(infile, "%s %s %c\n",input.vname, input.equal, &ch);
           //printf ("%d %s %s %c \n", i, input.vname, input.equal, ch);
        } 
        i++;
  } while (i < 11);

  fclose(infile);
/*************************************************************/
/*                   LOAD INITIAL DATA VALUES                */
/*************************************************************/
  PetscErrorCode     ierr;
  PetscMPIInt        size, rank;
  PetscReal          a = Array[0];//1.0;                       // shrink : make smaller in size or amount; to cause to shrink or contract; reduce.    
  PetscInt           Nx = Array[1], Ny = Nx, Nz = Nx ;    // Unit Cell dimensions or Unit cell resolution
  PetscReal          Sx = Array[2], Sy = Sx, Sz = Sx;
  PetscReal          dx = Sx/Nx, dy = Sy/Ny, dz = Sz/Nz;
  PetscScalar        U[Nx * Ny * Nz];                // UNIT CELL ARRAY
  PetscScalar        AC[Nx * Ny * Nz];               // FFTW SPATIAL HARMONICS
/* Truncate array, NM x NN x NP*/
  PetscInt           NM = Array[3], NN = NM, NP = NM;      // Truncate array size, Number of Spatial Harmonics
  PetscInt           NK = NM * NN * NP;              // Total Number of Spatial Harmonics on the 3D Truncation Array
  PetscScalar        AMNP[NK];                       // TRUNCATED ARRAY OF FFTW SPATIAL HARMONICS
  //PetscReal          KX[NK], KY[NK], KZ[NK];         // GRATING VECTORS ARRAYS
/* Number of Unit cells used to build the Lattice Grid, NPx x NPy x NPz  */
  //PetscInt           NPx = Array[4], NPy = NPx, NPz = NPx;  //9, 11;  // Carefull, NP? = 9, 11, etc didn't run on the Desktop Computer
/* Number of points per Unit Cell on the Lattice Grid NGP*/
  //PetscInt           NGP = Array[5];                        
/* Switch Coordinate System CYLINDRIC (C) & SPHERICAL(S) */
  //char Switch_Cood_Sys = ch;//'C'; //For print : %s is for strings (char*) and %c is for single characters (char).
/* New Matrix size - (Number of points for each unit cell on the lattice grid) * (Number of Unit cells used to build the Lattice Grid)*/ 
  //PetscInt           New_Nx = NGP * NPx; //update rows size matrix Nx
  //PetscInt           New_Ny = NGP * NPy; //update columns size matrix Ny 
  //PetscInt           New_Nz = NGP * NPz; //update heights size matrix Nz
/* New Matrix size - Step size */
  //PetscReal          New_dx = NPx*(Sx/New_Nx); //update step size dx
  //PetscReal          New_dy = NPy*(Sy/New_Ny); //update step size dy
  //PetscReal          New_dz = NPz*(Sz/New_Nz); //update step size dz
/* FDDER Operator - Size, step size and Boundary conditions */
  //PetscInt           NS[3]  = {New_Nx, New_Ny, New_Nz}; /* Size */
  //PetscReal          RES[3] = {New_dx, New_dy, New_dz}; /* Step size */ 
  //PetscInt           BC[3]  = {Array[6], Array[7], Array[8]};//{1, 1, 1};                /* Boundary Conditions : 0 = Dirichlet, -1 = Periodic, +1 = Neumann */ 
  //PetscInt           M      = NS[0] * NS[1] * NS[2];         /* Operator size M =  New_Nx * New_Ny * New_Nz;  */      

/***************************************************/
/*              START PETSC                        */
/***************************************************/
  ierr  = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  ierr  = MPI_Comm_rank(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
  ierr  = MPI_Comm_size(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"         SPACIAL VARIANT LATTICE       \n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"============= START PROGRAM ===========\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");CHKERRQ(ierr);

// Basic data:q
 /*  ierr = PetscPrintf(PETSC_COMM_WORLD,"  CYLINDRIC (C) & SPHERICAL(S)\n");CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"      SWITCH COORDINATE SYSTEM = %c\n", Switch_Cood_Sys);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"                UNIT CELL SIZE = (%d, %d, %d)\n", Nx, Ny, Nz);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"   NUMBER OF SPATIAL HARMONICS = (%d, %d, %d)\n", NM, NN, NP);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"NUMBER OF UNIT CELL ON LATTICE = (%d, %d, %d)\n", NPx, NPy, NPz);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"NUMBER OF POINTS PER UNIT CELL = %d\n", NGP);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"          OPERATOR MATRIX SIZE = (%d, %d, %d)\n", New_Nx, New_Ny, New_Nz);CHKERRQ(ierr);  
   ierr = PetscPrintf(PETSC_COMM_WORLD,"         FDDER OPERATOR SIZE M = %d\n", M);CHKERRQ(ierr); 
   ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");CHKERRQ(ierr);
*/
/* Start Elapse time - Full program */
  clock_t tic_start_program = clock();

/*************************************************************/
/*                BUILD 3D UNIT CELL DIVICE                  */
/*************************************************************/
  ierr = SVL_3D_UNIT_CELL(Nx, Ny, Nz, U, a, dx, dy, dz);CHKERRQ(ierr); 

/*************************************************************/
/*              METHOD 1 : FORWARD FFTW & SWAP               */
/*              1. 3D FFT over a 3D array                    */
/*************************************************************/
 ierr = SVL_3D_FFTW_SWAP_METHOD_1(Nx, Ny, Nz, U, AC);CHKERRQ(ierr);

/******************************************************************************/
/*                         METHOD 2 : FORWARD FFTW & SWAP                     */ 
/*                         1. 1D FFT on each vector Height (i)                */
/*                         2. 1D FFT on each vector Column (j)                */
/*                         3. 1D FFT on each vector Row    (k)                */                                                                        
/*                                                                            */
/* 1. Transpose 1: ALONG Y                                                    */                           
/*    1D FFT ALONG Y                                                          */
/*                                                                            */
/* 2. Transpose 2: ALONG Z                                                    */
/*    1D FFT ALONG Z                                                          */
/*                                                                            */ 
/* 3. Transpose 3: ALONG X                                                    */
/*    1D FFT ALONG X                                                          */
/*    Repeat Transpose 3: ALONG Z                                             */  
/******************************************************************************/
 //ierr = SVL_3D_FFTW_SWAP_METHOD_2(Nx, Ny, Nz, U, AC);CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 3 : FORWARD FFTW & SWAP               */
/*             1. 2D FFT a long slap (Rows (k), Columns (j)) */
/*             2. 1D FFT a long row Height (i)               */
/*************************************************************/
 //ierr = SVL_3D_FFTW_SWAP_METHOD_3(Nx, Ny, Nz, U, AC);CHKERRQ(ierr);

/*************************************************************/
/*               TRUCATE THE SPATIAL HARMONICS               */
/*************************************************************/
  ierr = SVL_3D_TRUNCATE_FFTW_ARRAY(Nx, Ny, Nz, NM, NN, NP, AC, AMNP);CHKERRQ(ierr);

/*************************************************************/
/*                        GRADING VECTOR                     */
/*************************************************************/
  //ierr = SVL_3D_GRADING_VECTOR(NM, NN, NP, Sx, Sy, Sz, KX, KY, KZ); CHKERRQ(ierr);

/*************************************************************/
/*                  IMPLEMENT IMPROVEMENTS                   */
/*************************************************************/
/* AMNP[NK] array is the input of size NK = NM*NN*NP and output of size NK << NM*NN*NP,
   INPUT  : Contain the truncated fftw of highest harmonics 
   OUTPUT : Contain a reduce number of highest harmonic mode from FFTW amplitudes and collinear

   NK = NM*NN*NP will get update, &NK parameters is passing by reference, the new value of NK << NM*NN*NP, must smaller
   NK is used to count the number of spatial harmonics amplitudes than pass the threshold and do not fall below that threshold
   In short, We eliminate all the gratings vectors and their associates spatial harmonic amplitudes that fall below the threshold spatial harmonic. */

  //ierr = SVL_3D_IMPLEMENT_IMPROVEMENTS(NM, NN, NP, &NK, KX, KY, KZ, AMNP);CHKERRQ(ierr); 
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Improvements - number of iteration = %d\n", NK);CHKERRQ(ierr);

/*************************************************************/
/*                    DEFINE SPATIAL VARIANT                */
/*   AAttribute 1 - ORIENTATION FUNCTION                     */
/*                 C = CYLINDRICAL (RHO, THETA. ZH)          */
/*                 S = SPHERICAL   (RHO, THETA, VARPHI)      */
/*  Attribute 2 - LATTICE SPACING FUNCTION  PER(r)           */
/*************************************************************/
/* Rotation Coordinates arrays  
  PetscReal       X[New_Nx * New_Ny * New_Nz], Y[New_Nx * New_Ny * New_Nz], Z[New_Nx * New_Ny * New_Nz];
  PetscReal     RSQ[New_Nx * New_Ny * New_Nz];
  PetscReal   THETA[New_Nx * New_Ny * New_Nz];
  PetscReal  VARPHI[New_Nx * New_Ny * New_Nz];
  PetscReal     PER[New_Nx * New_Ny * New_Nz];
  PetscReal      ZZ[New_Nx * New_Ny * New_Nz];
*/
/* Switch Coordinate System CYLINDRIC (C) & SPHERICAL(S) 
  switch(Switch_Cood_Sys) {
         case 'C':
           ierr = SVL_3D_CYLINDRIC_SPATIAL_VARIANT(New_Nx, New_Ny, New_Nz, New_dx, New_dy, New_dz, X, Y, Z, RSQ, THETA, ZZ, PER);CHKERRQ(ierr);
           break;
         case 'S':
           ierr = SVL_3D_SPHERICAL_SPATIAL_VARIANT(New_Nx, New_Ny, New_Nz, New_dx, New_dy, New_dz, X, Y, Z, RSQ, THETA, VARPHI, PER);CHKERRQ(ierr);
           break;
        default:
           ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTRY AGAIN!!!, YOU DID NOT IMPLEMENTED ANY ATTRIBUTES \n");CHKERRQ(ierr);
           break;
  }*/

/*************************************************************/
/*     FDDER - Finite-Difference DERivative Operators        */
/*************************************************************/ 
/* Matrix, Mat Array D operator */
  //Mat   D;                                       

/* Create operator array D = [DX / DY/ DZ] with size (rows, columns) = (3M, M) */ 
  //ierr = MatCreate(PETSC_COMM_WORLD,&D);CHKERRQ(ierr);
  //ierr = MatSetSizes(D,PETSC_DECIDE,PETSC_DECIDE, 3*M, M);CHKERRQ(ierr); /* change the number: 2 for 2D and 3 for 3D */
  
/* Set up format opeartor  D = [DX / DY/ DZ] */ 
  //ierr = MatSetFromOptions(D);CHKERRQ(ierr); 
  //ierr = MatSetUp(D);CHKERRQ(ierr);

/* Fill operator D */
  //ierr = SVL_3D_FDDER(NS, BC, RES, D);CHKERRQ(ierr);

/************************************************************************/
/*                   Show the operator D on screen                      */
/************************************************************************/    
   /*PetscInt   m, n;
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D on MAIN\n----------------\n");
   ierr = MatView(D,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   ierr = MatGetSize(D,&m,&n);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nD[cols, rows] = [%d %d]\n", m, n);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nM = %d\n", M);*/

/*************************************************************/
/*                     START MAIN LOOP                       */
/*************************************************************/
/* Start Elapse time - LOOP program */
  clock_t tic_start_loop = clock();

  //ierr = SVL_3D_LOOP(Switch_Cood_Sys, M, New_Nx, New_Ny, New_Nz, NK, AMNP, KX, KY, KZ, THETA, VARPHI, PER, D);CHKERRQ(ierr);

/*************************************************************/
/*                      END MAIN LOOP                       */
/*************************************************************/

/************************************************************************/
/*                           Exiting the program                        */
/************************************************************************/
/* "break" is a keyword, which causes an immediate exit from the switch */
/*  or loop (for, while or do),                                         */ 
/*  "exit() is a standard library function, which terminates program   */
/*  execution when it is called.                                        */
/************************************************************************/
/* PetscInt aa = 1; // 1 Exit program , 2 Continue program 
   if(aa == 1) {  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Exiting the program....\n");CHKERRQ(ierr);
      exit(1);
   } 
*/
/*************************************************************/
/*              END PROGRAM and Find The Elapsed Time        */
/*************************************************************/
    clock_t toc = clock(); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n============= END PROGRAM ===========\n");
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n    Elapsed SERIAL: %f seconds\n", (double)(tic_start_loop - tic_start_program) / CLOCKS_PER_SEC);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n    Elapsed LOOP  : %f seconds\n", (double)(toc - tic_start_loop) / CLOCKS_PER_SEC);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n    Elapsed FULL  : %f seconds\n", (double)(toc - tic_start_program) / CLOCKS_PER_SEC);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n=====================================\n");
    ierr = PetscFinalize();
    return 0;
}





