// UNIT CELL
//extern PetscErrorCode SVL_3D_UNIT_CELL_ARRAY(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal, PetscReal, PetscReal, PetscReal); //GROUP LEAD PROGRAM
//extern PetscErrorCode SVL_3D_ZERO_CELL(PetscInt, PetscInt, PetscInt, PetscReal*);
extern PetscErrorCode SVL_3D_UNIT_CELL(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscReal, PetscReal, PetscReal, PetscReal);
// FFTW & SWAP METHOD AND TRANSPOSE
extern PetscErrorCode SVL_3D_FFTW_SWAP_METHOD_1(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
extern PetscErrorCode SVL_3D_FFTW_SWAP_METHOD_2(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
//extern PetscErrorCode SVL_3D_FFTW_SWAP_METHOD_3(PetscInt, PetscInt, PetscInt, PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_3D_TRANSPOSE_1_COMPLEX(PetscInt, PetscInt, PetscInt, PetscScalar*);
extern PetscErrorCode SVL_3D_TRANSPOSE_2_COMPLEX(PetscInt, PetscInt, PetscInt, PetscScalar*);
extern PetscErrorCode SVL_3D_TRANSPOSE_3_COMPLEX(PetscInt, PetscInt, PetscInt, PetscScalar*);
extern PetscErrorCode SVL_3D_TRANSPOSE_4_COMPLEX(PetscInt, PetscInt, PetscInt, PetscScalar*);
// FFTW 1D AND 3D, IFFWT AND SWAP
extern PetscErrorCode SVL_1D_FFTW(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
extern PetscErrorCode SVL_3D_FFTW(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
extern PetscErrorCode SVL_3D_IFFTW(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
extern PetscErrorCode SVL_3D_SWAP_QUADRANTS(PetscInt, PetscInt, PetscInt, PetscScalar*);
// TRUNCATE 
extern PetscErrorCode SVL_3D_TRUNCATE_FFTW_ARRAY(PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
// SAVE
extern PetscErrorCode SAVE_1D_TO_3D_ARRAY_REAL(const char*, PetscInt, PetscInt, PetscInt, PetscReal*);
extern PetscErrorCode SAVE_1D_TO_3D_ARRAY_COMPLEX(const char*, const char*, PetscInt, PetscInt, PetscInt, PetscScalar*);

// GRADING
//extern PetscErrorCode SVL_3D_GRADING_VECTOR(PetscInt, PetscInt, PetscInt, PetscReal, PetscReal, PetscReal, PetscReal*, PetscReal*, PetscReal*);
// IMPROVEMENTS
//extern PetscErrorCode SVL_3D_IMPLEMENT_IMPROVEMENTS(PetscInt, PetscInt, PetscInt, PetscInt*, PetscReal*, PetscReal*, PetscReal*, PetscScalar*);  //GROUP LEAD PROGRAM
//extern PetscErrorCode SVL_3D_ELIMINATE_GRATING_ACCORD_THEIR_AMPLITUD(PetscInt, PetscInt, PetscInt, PetscScalar*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_IDENTIFIED_COLLINEAR_PLANAR_GRATING(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_IMPROVEMENTS(PetscInt, PetscInt,  PetscInt, PetscInt*, PetscScalar*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*); 
// SPHERICAL SPATIAL_VARIANT
//extern PetscErrorCode SVL_3D_SPHERICAL_SPATIAL_VARIANT(PetscInt, PetscInt, PetscInt, PetscReal, PetscReal, PetscReal, PetscReal*, PetscReal*, PetscReal*, PetscReal*,  PetscReal*, PetscReal*, PetscReal*);  //GROUP LEAD PROGRAM
//extern PetscErrorCode SVL_3D_SPHERICAL_ORIENTATION_FUNCTION(PetscInt, PetscInt, PetscInt, PetscReal, PetscReal, PetscReal, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_SPHERICAL_LATTICE_SPACING_FUNCTION(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*);
// CYLIDRICAL SPATIAL_VARIANT
//extern PetscErrorCode SVL_3D_CYLINDRIC_SPATIAL_VARIANT(PetscInt, PetscInt, PetscInt, PetscReal, PetscReal, PetscReal, PetscReal*, PetscReal*, PetscReal*, PetscReal*,  PetscReal*, PetscReal*, PetscReal*);  //GROUP LEAD PROGRAM
//extern PetscErrorCode SVL_3D_CYLINDRICAL_ORIENTATION_FUNCTION(PetscInt, PetscInt, PetscInt, PetscReal, PetscReal, PetscReal, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_CYLINDRICAL_LATTICE_SPACING_FUNCTION(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*);
// FDDER
//extern PetscErrorCode SVL_3D_FDDER(PetscInt*, PetscInt*, PetscReal*, Mat);
// LOOP
//extern PetscErrorCode SVL_3D_LOOP(char, PetscInt, PetscInt,  PetscInt, PetscInt, PetscInt,  PetscScalar*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, Mat); //GROUP LEAD PROGRAM
//extern PetscErrorCode SVL_3D_ORIENTATION_VECTOR(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal, PetscReal, PetscReal);
//SPHERICAL
//extern PetscErrorCode SVL_3D_SPHERICAL_TRANSLATION(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*); 
//extern PetscErrorCode SVL_3D_CART2SPHER(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_SPHER_ROTATION(PetscInt, PetscInt,  PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_SPHER_SPACING(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_SPHER2CART(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//CYLINCRICAL
//extern PetscErrorCode SVL_3D_CYLINDRICAL_TRANSLATION(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*); 
//extern PetscErrorCode SVL_3D_CART2CYLIN(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_CYLIN_ROTATION(PetscInt, PetscInt,  PetscInt, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_CYLIN_SPACING(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_3D_CYLIN2CART(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//SOLVE LINEAR SYSTEM
//extern PetscErrorCode SVL_3D_RHS(PetscInt, PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal* ); 

