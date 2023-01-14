/*************************************************************/
/*      FDDER - Finite-Difference DERivative Operators       */
/*************************************************************/
/*                      OPERATOR FDDER                       */
/*************************************************************/
/* [DX, D2X, DY, D2Y, DZ, D2Z] = fdder(NS,RES,BC);

 This function generates matrix derivative operators
 for scalar functions on a collocated grid.

 INPUT ARGUMENTS
 ================
 NS        [Nx Ny Nz] size of grid
 RES       [dx dy dz] grid resolution
 BC        [BCx BCy BCz] Boundary Conditions
           0 = Dirichlet, -1 = Periodic, +1 = Neumann

 OUTPUT ARGUMENTS
 ================
 DX        First-order derivative with respect to x
 D2X       Second-order derivative with respect to x
 DY        First-order derivative with respect to y
 D2Y       Second-order derivative with respect to y
 DZ        First-order derivative with respect to a
 D2Z       Second-order derivative with respect to z
*/ 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_FDDER"
PetscErrorCode SVL_3D_FDDER(PetscInt* NS, PetscInt* BC, PetscReal* RES, Mat D){
  //Mat            DX, DY, DZ;
  Mat            D2X, D2Y, D2Z;  
  //Mat            D;
  PetscErrorCode ierr;
  PetscInt       i, k;
  //PetscMPIInt    size,rank;
  PetscInt       col1[2], col2[3];
  PetscScalar    val1[2], val2[3];
  //PetscBool      nonzeroguess = PETSC_FALSE;
  PetscMPIInt    n1, n2; //n3;
  PetscMPIInt    ny, nz; //nx; 
  PetscMPIInt    mxlo, mxhi, mylo, myhi, mzlo, mzhi;
  PetscScalar    CX, C2X, CY, C2Y, CZ, C2Z;


/* Finite Difference DERivative */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 12: FDDER (Finite Difference DERivative)\n");CHKERRQ(ierr);

/***************************************************/
/*       HANDLE INPUT AND OUTPUT ARGUMENTS         */
/***************************************************/
//int NS[3]  = {3, 3, 3};
//int RES[3] = {1, 1, 1};
//int BC[3]  = {1, 1, 1};   /* 0 = Dirichlet, -1 = Periodic, +1 = Neumann */ 

/* DETERMINE SIZE OF GRID */
 PetscInt  Nx = NS[0];
 PetscInt  Ny = NS[1];
 PetscInt  Nz = NS[2];

/* COMPUTE SIZE OF MATRICES */
 PetscInt  M = Nx * Ny * Nz;

/* DIFFERENCE QUOTIENT APPROXIMATION */
 CX  = (double) 1/(2*RES[0]);
 C2X = (double) 1/(RES[0]*RES[0]);
 CY  = (double) 1/(2*RES[1]);
 C2Y = (double) 1/(RES[1]*RES[1]);
 CZ  = (double) 1/(2*RES[2]);
 C2Z = (double) 1/(RES[2]*RES[2]);
/***************************************************/
/*                  BEGIN PETSC                    */
/***************************************************/
  PetscFunctionBeginUser;

/***************************************************/
/*              START PETSC                        */
/***************************************************/
  //PetscInitialize(&argc,&args,(char*)0,help);
  //ierr  = MPI_Comm_rank(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
  //ierr  = MPI_Comm_size(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

  //if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  //ierr = PetscOptionsGetInt(NULL,"-M",&M,NULL);CHKERRQ(ierr);
  //ierr = PetscOptionsGetBool(NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);

/* MatCreate : Creates a matrix where the type is determined from either a call to MatSetType()
   or from the options database with a call to MatSetFromOptions(). The default matrix type is
   AIJ, using the routines MatCreateSeqAIJ() or MatCreateAIJ() if you do not set a type in the
   options database. If you never call MatSetType() or MatSetFromOptions() it will generate an
   error when you try to use the matrix. 
   Performance tuning note:  For problems of substantial size, preallocation of matrix memory 
   is crucial for attaining good performance. See the matrix chapter of the users manual for details. 
   MatSetSizes : Sets the local and global sizes, and checks to determine compatibility */

/* Create the matrix DX & D2X*/ 
  /*ierr = MatCreate(PETSC_COMM_WORLD,&DX);CHKERRQ(ierr);
  ierr = MatSetSizes(DX,PETSC_DECIDE,PETSC_DECIDE,M,M);CHKERRQ(ierr);*/
  ierr = MatCreate(PETSC_COMM_WORLD,&D2X);CHKERRQ(ierr);
  ierr = MatSetSizes(D2X,PETSC_DECIDE,PETSC_DECIDE,M,M);CHKERRQ(ierr);

/* Create the matrix DY & D2Y */ 
  /*ierr = MatCreate(PETSC_COMM_WORLD,&DY);CHKERRQ(ierr);
  ierr = MatSetSizes(DY,PETSC_DECIDE,PETSC_DECIDE,M,M);CHKERRQ(ierr);*/
  ierr = MatCreate(PETSC_COMM_WORLD,&D2Y);CHKERRQ(ierr);
  ierr = MatSetSizes(D2Y,PETSC_DECIDE,PETSC_DECIDE,M,M);CHKERRQ(ierr);

/* Create the matrix DZ & D2Z */ 
 /* ierr = MatCreate(PETSC_COMM_WORLD,&DZ);CHKERRQ(ierr);
  ierr = MatSetSizes(DZ,PETSC_DECIDE,PETSC_DECIDE,M,M);CHKERRQ(ierr);*/
  ierr = MatCreate(PETSC_COMM_WORLD,&D2Z);CHKERRQ(ierr);
  ierr = MatSetSizes(D2Z,PETSC_DECIDE,PETSC_DECIDE,M,M);CHKERRQ(ierr);

/* Create the matrix D = [DX / DY / DZ] */ 
//  ierr = MatCreate(PETSC_COMM_WORLD,&D);CHKERRQ(ierr);
//  ierr = MatSetSizes(D,PETSC_DECIDE,PETSC_DECIDE,3*M,M);CHKERRQ(ierr);

/* MatSetFromOptions : Creates a matrix where the type is determined from the options database. Generates a
   parallel MPI matrix if the communicator has more than one processor. The default matrix type is AIJ,
   using the routines MatCreateSeqAIJ() and MatCreateAIJ() if you do not select a type in the options database   
   MatSetUp: Sets up the internal matrix data structures for the later use.  */

/* Set up DX & D2X */ 
  /*ierr = MatSetFromOptions(DX);CHKERRQ(ierr); 
  ierr = MatSetUp(DX);CHKERRQ(ierr);*/
  ierr = MatSetFromOptions(D2X);CHKERRQ(ierr); 
  ierr = MatSetUp(D2X);CHKERRQ(ierr);

/* Set up DY & D2Y*/ 
  /*ierr = MatSetFromOptions(DY);CHKERRQ(ierr); 
  ierr = MatSetUp(DY);CHKERRQ(ierr);*/
  ierr = MatSetFromOptions(D2Y);CHKERRQ(ierr); 
  ierr = MatSetUp(D2Y);CHKERRQ(ierr);

/* Set up DZ & D2Z*/ 
  /*ierr = MatSetFromOptions(DZ);CHKERRQ(ierr); 
  ierr = MatSetUp(DZ);CHKERRQ(ierr);*/
  ierr = MatSetFromOptions(D2Z);CHKERRQ(ierr); 
  ierr = MatSetUp(D2Z);CHKERRQ(ierr);

/* Set up D = [DX / DY / DZ]*/ 
 // ierr = MatSetFromOptions(D);CHKERRQ(ierr); 
 // ierr = MatSetUp(D);CHKERRQ(ierr);

/***************************************************/
/*                Assemble matrix                  */
/*                   DX and D2X                    */
/***************************************************/
// CREATE MATRIX IF DIMENSION EXISTS

// Build Default DX D2X
if (Nx > 1 ) { 
  val1[0] = -1.0*CX;
  val1[1] =  1.0*CX;

  val2[0] =  1.0*C2X;
  val2[1] = -2.0*C2X;
  val2[2] =  1.0*C2X;

  for (i = 1; i < M - 1; i++) { /*i = row, col1, col2 = column */
      col1[0] = i-1;
      col1[1] = i+1;
      //ierr   = MatSetValues(DX,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr); 
      ierr   = MatSetValues(D,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr); 
      col2[0] = i-1;
      col2[1] = i;
      col2[2] = i+1;
      ierr   = MatSetValues(D2X,1,&i,3,col2,val2,INSERT_VALUES);CHKERRQ(ierr);
  }  

  i = 0; 
  col1[0] = 0;
  col1[1] = 1; 
  val1[0] = 0*CX; 
  val1[1] = 1*CX;
  //ierr = MatSetValues(DX,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(D,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr);
  i = M-1;
  col1[0] = M - 2;
  col1[1] = M - 1 ;
  val1[0] = -1*CX; 
  val1[1] =  0*CX;
  //ierr = MatSetValues(DX,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(D,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr);
  
  i = 0;
  col1[0] = 0;
  col1[1] = 1 ;
  val1[0] = -2; 
  val1[1] =  1;
  ierr = MatSetValues(D2X,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr);


  i = M-1;
  col1[0] = M - 2;
  col1[1] = M - 1 ;
  val1[0] =  1; 
  val1[1] = -2;
  ierr = MatSetValues(D2X,1,&i,2,col1,val1,INSERT_VALUES);CHKERRQ(ierr);


/* Fix Boundary Errors (default to Dirichlet) */
    for (ny = 1; ny <=  Ny*Nx; ny++) {
        mxlo = (ny - 1)*Nx + 1; // row
        mxhi = (ny - 1)*Nx + Nx; // column

        if (mxlo > 1) {
            //ierr = MatSetValue(DX,  mxlo-1, mxlo-2,  0*CX, INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(D,  mxlo-1, mxlo-2,  0*CX, INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(D2X, mxlo-1, mxlo-2, 0*C2X, INSERT_VALUES);CHKERRQ(ierr); 
        }
        if (mxhi < M) {
            //ierr = MatSetValue(DX,  mxhi-1, mxhi,  0*CX, INSERT_VALUES);CHKERRQ(ierr); 
            ierr = MatSetValue(D,  mxhi-1, mxhi,  0*CX, INSERT_VALUES);CHKERRQ(ierr); 
            ierr = MatSetValue(D2X, mxhi-1, mxhi, 0*C2X, INSERT_VALUES);CHKERRQ(ierr);
        }
    }
/* Incorporate Boundary Conditions */
    switch (BC[0]) {
        case -1: //Periodic
            for (ny = 1; ny <= Ny*Nx; ny++) {
                n1 = (ny - 1)*Nx + 1;
                n2 = (ny - 1)*Nx + Nx;
                //ierr = MatSetValue(DX,  n1-1, n2-1, -1*CX, INSERT_VALUES);CHKERRQ(ierr);
                //ierr = MatSetValue(DX,  n2-1, n1-1,  1*CX, INSERT_VALUES);CHKERRQ(ierr); 
                ierr = MatSetValue(D,  n1-1, n2-1, -1*CX, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D,  n2-1, n1-1,  1*CX, INSERT_VALUES);CHKERRQ(ierr); 
                ierr = MatSetValue(D2X, n1-1, n2-1, 1*C2X, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2X, n2-1, n1-1, 1*C2X, INSERT_VALUES);CHKERRQ(ierr);
            }
            break;
        case +1: //Neumann
            for (ny = 1; ny <= Ny*Nx; ny++) {
                n1 = (ny - 1)*Nx + 1;
                n2 = (ny - 1)*Nx + Nx;
                //ierr = MatSetValue(DX,  n1-1, n1-1, -2*CX, INSERT_VALUES);CHKERRQ(ierr); 
                //ierr = MatSetValue(DX,  n1-1, n1  ,  2*CX, INSERT_VALUES); CHKERRQ(ierr);
                //ierr = MatSetValue(DX,  n2-1, n2-1,  2*CX, INSERT_VALUES);CHKERRQ(ierr);
                //ierr = MatSetValue(DX,  n2-1, n2-2, -2*CX, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D,  n1-1, n1-1, -2*CX, INSERT_VALUES);CHKERRQ(ierr); 
                ierr = MatSetValue(D,  n1-1, n1  ,  2*CX, INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(D,  n2-1, n2-1,  2*CX, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D,  n2-1, n2-2, -2*CX, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2X, n1-1, n1-1, 0*C2X, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2X, n1-1, n1  , 0*C2X, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2X, n2-1, n2-1, 0*C2X, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2X, n2-1, n2-2, 0*C2X, INSERT_VALUES);CHKERRQ(ierr);
            }
            break;
    }// end swicth

 } //end if (Nx > 1)

/***************************************************/
/*                Assemble matrix                  */
/*                   DY and D2Y                    */
/***************************************************/
// CREATE MATRIX IF DIMENSION EXISTS
if (Ny > 1) {
   for (k = 0; k < Nz; k++) {
        for (i = 0; i < Nx*Ny - Nx; i++) { 
            //ierr = MatSetValue(DY,  i +  k * Nx*Ny     , i + Nx +  k * Nx*Ny, 1*CY, INSERT_VALUES);CHKERRQ(ierr);
            //ierr = MatSetValue(DY,  i + Nx +  k * Nx*Ny, i +  k * Nx*Ny     ,-1*CY, INSERT_VALUES);CHKERRQ(ierr);
            ierr = MatSetValue(D,  i +  k * Nx*Ny + M     , i + Nx +  k * Nx*Ny, 1*CY, INSERT_VALUES);CHKERRQ(ierr);
            ierr = MatSetValue(D,  i + Nx +  k * Nx*Ny + M, i +  k * Nx*Ny     ,-1*CY, INSERT_VALUES);CHKERRQ(ierr);

            ierr = MatSetValue(D2Y, i +  k * Nx*Ny     , i + Nx +  k * Nx*Ny, 1*C2Y, INSERT_VALUES);CHKERRQ(ierr);
            ierr = MatSetValue(D2Y, i + Nx +  k * Nx*Ny, i +  k * Nx*Ny     , 1*C2Y, INSERT_VALUES);CHKERRQ(ierr);
         }
   }
  for (i = 0; i < M; i++) {
      ierr = MatSetValue(D2Y, i, i, -2*C2Y, INSERT_VALUES);CHKERRQ(ierr);
  } 

// Incorporate Boundary Conditions
    switch (BC[1]){
        case -1: //periodic
            for (k = 0; k < Nz; k++) {
                for (ny = 1; ny <= Nx; ny++) {
                    mylo =  ny - 1 + k * Nx*Ny;
                    myhi = (Ny - 1)*Nx + ny - 1 + k * Nx*Ny;
                    //ierr = MatSetValue(DY,  mylo, myhi, -1*CY, INSERT_VALUES);CHKERRQ(ierr);
                    //ierr = MatSetValue(DY,  myhi, mylo,  1*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D,  mylo + M, myhi, -1*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D,  myhi + M, mylo,  1*CY, INSERT_VALUES);CHKERRQ(ierr);

                    ierr = MatSetValue(D2Y, mylo, myhi, 1*C2Y, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D2Y, myhi, mylo, 1*C2Y, INSERT_VALUES);CHKERRQ(ierr);
                 }
            }
          break;
        case +1: //extended
            for (k = 0; k < Nz; k++) {
                for (ny = 1; ny <= Nx; ny++) {
                    mylo =  ny - 1 + k * Nx*Ny;
                    myhi = (Ny - 1)*Nx + ny - 1 + k * Nx*Ny;
                    /*ierr = MatSetValue(DY,  mylo, mylo   , -2*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(DY,  mylo, mylo+Nx,  2*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(DY,  myhi, myhi   ,  2*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(DY,  myhi, myhi-Nx, -2*CY, INSERT_VALUES);CHKERRQ(ierr);
*/
                    ierr = MatSetValue(D,  mylo + M, mylo   , -2*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D,  mylo + M, mylo+Nx,  2*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D,  myhi + M, myhi   ,  2*CY, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D,  myhi + M, myhi-Nx, -2*CY, INSERT_VALUES);CHKERRQ(ierr);

                    ierr = MatSetValue(D2Y, mylo, mylo   , 0*C2Y, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D2Y, mylo, mylo+Nx, 0*C2Y, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D2Y, myhi, myhi   , 0*C2Y, INSERT_VALUES);CHKERRQ(ierr);
                    ierr = MatSetValue(D2Y, myhi, myhi-Nx, 0*C2Y, INSERT_VALUES);CHKERRQ(ierr);
                }
            }
          break;
    }// end switch
} //end if (Ny > 1)


/***************************************************/
/*                Assemble matrix                  */
/*                   DZ and D2Z                    */
/***************************************************/
// CREATE MATRIX IF DIMENSION EXISTS
if (Nz > 1) {
  for (i = 0; i < M - Nx*Ny; i++) { 
     //ierr = MatSetValue(DZ,  i      , i+Nx*Ny , 1*CZ, INSERT_VALUES);CHKERRQ(ierr);
     //ierr = MatSetValue(DZ,  i+Nx*Ny, i,-1*CZ , INSERT_VALUES);CHKERRQ(ierr);
     ierr = MatSetValue(D,  i + 2*M      , i+Nx*Ny , 1*CZ, INSERT_VALUES);CHKERRQ(ierr);
     ierr = MatSetValue(D,  i+Nx*Ny + 2*M, i,-1*CZ , INSERT_VALUES);CHKERRQ(ierr);

     ierr = MatSetValue(D2Z, i      , i+Nx*Ny , 1*C2Z, INSERT_VALUES);CHKERRQ(ierr);
     ierr = MatSetValue(D2Z, i+Nx*Ny, i, 1*C2Z, INSERT_VALUES);CHKERRQ(ierr);
  }
  for (i = 0; i < M; i++) {
      ierr = MatSetValue(D2Z, i, i, -2*C2Z, INSERT_VALUES);CHKERRQ(ierr);
  } 

// Incorporate Boundary Conditions
    switch (BC[1]){
        case -1: //periodic
            for (nz = 1; nz <= Nx*Ny; nz++) {
                mzlo =  nz - 1;
                mzhi = (Ny - 1)*Nx*Ny + nz - 1;
                //ierr = MatSetValue(DZ , mzlo, mzhi,-1*CZ, INSERT_VALUES);CHKERRQ(ierr);
                //ierr = MatSetValue(DZ , mzhi, mzlo, 1*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D , mzlo + 2*M, mzhi,-1*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D , mzhi + 2*M, mzlo, 1*CZ, INSERT_VALUES);CHKERRQ(ierr);

                ierr = MatSetValue(D2Z, mzlo, mzhi, 1*C2Z, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2Z, mzhi, mzlo, 1*C2Z, INSERT_VALUES);CHKERRQ(ierr);
            }
          break;
        case +1: //extended
            for (nz = 1; nz <= Nx*Ny; nz++) {
                mzlo = nz - 1;
                mzhi = (Ny - 1)*Nx*Ny + nz - 1;
                //ierr = MatSetValue(DZ , mzlo, mzlo       ,- 2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                //ierr = MatSetValue(DZ , mzlo, mzlo+Nx*Ny ,  2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                //ierr = MatSetValue(DZ , mzhi, mzhi       ,  2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                //ierr = MatSetValue(DZ , mzhi, mzhi-Nx*Ny , -2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D , mzlo + 2*M, mzlo       ,- 2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D , mzlo + 2*M, mzlo+Nx*Ny ,  2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D , mzhi + 2*M, mzhi       ,  2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D , mzhi + 2*M, mzhi-Nx*Ny , -2*CZ, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2Z, mzlo, mzlo       , 0*C2Z, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2Z, mzlo, mzlo+Nx*Ny , 0*C2Z, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2Z, mzhi, mzhi       , 0*C2Z, INSERT_VALUES);CHKERRQ(ierr);
                ierr = MatSetValue(D2Z, mzhi, mzhi-Nx*Ny , 0*C2Z, INSERT_VALUES);CHKERRQ(ierr);
            }
          break;
    }// end switch
}//end if (Nz > 1)

/***************************************************/
/*                 ASSEMBLY MATRIX                 */
/***************************************************/
/* MatAssemblyBegin : Begins assembling the matrix. This routine should be called after completing
                      all calls to MatSetValues(). 
   MatAssemblyEnd : Completes assembling the matrix. This routine should be called after MatAssemblyBegin(). 
*/

/* Assembly  DX & D2X */
//  ierr = MatAssemblyBegin(DX,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
//  ierr = MatAssemblyEnd(DX,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
  ierr = MatAssemblyBegin(D2X,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D2X,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
/* Assembly  DY & D2Y */
//  ierr = MatAssemblyBegin(DY,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(DY,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(D2Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D2Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
/* Assembly  DZ & D2Z */
//  ierr = MatAssemblyBegin(DZ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(DZ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(D2Z,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D2Z,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

/* Assembly D */ 
  ierr = MatAssemblyBegin(D,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

/* MatScale: Scales all elements of a matrix by a given number. */
  /*MatScale(DX, (double) 1/(2*RES[0]));
  MatScale(D2X,(double) 1/(RES[0]*RES[0]));
  MatScale(DY, (double) 1/(2*RES[1]));
  MatScale(D2Y,(double) 1/(RES[1]*RES[1]));
  MatScale(DZ, (double) 1/(2*RES[2]));
  MatScale(D2Z,(double) 1/(RES[2]*RES[2]));
*/

/***************************************************/
/*             SHOW MATLAB FILE ON SCREEN          */
/***************************************************/
/* MatView : Visualizes a matrix object. view the matrix - Visualizes a matrix object. PETSC_VIEWER_STDOUT_SELF - standard output (default)  */

/* Show the operator DX DY DZ D2X D2Y D2Z on screen */
/* ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix DX\n---------------\n");
  ierr = MatView(DX,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D2X\n----------------\n");
  ierr = MatView(D2X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix DY\n---------------\n");
  ierr = MatView(DY,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D2Y\n----------------\n");
  ierr = MatView(D2Y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix DZ\n---------------\n");
  ierr = MatView(DZ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D2Z\n----------------\n");
  ierr = MatView(D2Z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D\n----------------\n");
  ierr = MatView(D,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
*/
/* Show the operator D into screen */
   //ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D FDDER\n----------------\n");
   //ierr = MatView(D,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

/***************************************************/
/*                SAVE INTO MATLAB FILE            */
/***************************************************/
/* Write D operator into a file */
/* PetscViewer : Abstract PETSc object that helps view (in ASCII, binary, graphically etc) other PETSc objects 
   PetscViewerASCIIOpen : Opens an ASCII file as a PetscViewer. 
   PetscViewerSetFormat: Sets the format for PetscViewers.  */

 /*  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"D.m", &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = MatView(DX,viewer);CHKERRQ(ierr);
  ierr = MatView(D2X,viewer);CHKERRQ(ierr);
  ierr = MatView(DY,viewer);CHKERRQ(ierr);
  ierr = MatView(D2Y,viewer);CHKERRQ(ierr);
  ierr = MatView(DZ,viewer);CHKERRQ(ierr);
  ierr = MatView(D2Z,viewer);CHKERRQ(ierr);
  ierr = MatView(D,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);*/


/************************************************************************/
/*                           View  Matrix                               */
/************************************************************************/
     /*PetscInt          m, n;
     ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D on FDDER\n----------------\n");
     ierr = MatView(D,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
     ierr = MatGetSize(D,&m,&n);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_WORLD,"\nD[cols, rows] = [%d %d]\n", m, n);
     ierr = PetscPrintf(PETSC_COMM_WORLD,"\nM = %d\n", M);*/

/* Free work space.  All PETSc objects should be destroyed when they are no longer needed. */
 // ierr = MatDestroy(&DX);CHKERRQ(ierr);
  ierr = MatDestroy(&D2X);CHKERRQ(ierr);
 // ierr = MatDestroy(&DY);CHKERRQ(ierr);
  ierr = MatDestroy(&D2Y);CHKERRQ(ierr);
//  ierr = MatDestroy(&DZ);CHKERRQ(ierr);
  ierr = MatDestroy(&D2Z);CHKERRQ(ierr);
  //ierr = MatDestroy(&D);CHKERRQ(ierr); /* Uncomment this line will destroy Mat D before go out of the subroutine - keep line on Commented */

/* Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary). */
  PetscFunctionReturn(0); 
}
