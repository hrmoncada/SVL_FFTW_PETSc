// UNIT CELL
//extern PetscErrorCode SVL_2D_UNIT_CELL(PetscInt, PetscInt, PetscScalar*);
//1extern PetscErrorCode SVL_2D_ZERO_CELL(PetscInt, PetscInt, PetscScalar*);
extern PetscErrorCode SVL_2D_TRIANGLE(PetscInt, PetscInt, PetscReal*);
// FFTW AND SWAP METHOD 1, 2 3  - 1D FFTW
extern PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_1(PetscInt, PetscInt, PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_2(PetscInt, PetscInt, PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_3(PetscInt, PetscInt, PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_1D_FFTW(PetscInt, PetscInt,  PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_2D_TRANSPOSE_REAL(PetscInt, PetscInt, PetscReal*);
extern PetscErrorCode SVL_2D_TRANSPOSE_COMPLEX(PetscInt, PetscInt, PetscScalar*);
// FFTW AND SWAP METHOD 4 - 2D FFTW
extern PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_4(PetscInt, PetscInt, PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_2D_FFTW(PetscInt, PetscInt, PetscReal*, PetscScalar*);
extern PetscErrorCode SVL_2D_IFFTW(PetscInt, PetscInt, PetscScalar*, PetscScalar*);
extern PetscErrorCode SVL_2D_SWAP_QUADRANTS(PetscInt, PetscInt, PetscScalar*);
// TRUNCATE
//extern PetscErrorCode SVL_2D_TRUNCATE_FFTW_SPATIAL_HARMONICS(PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
extern PetscErrorCode SVL_2D_TRUNCATE_FFTW_ARRAY(PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar*, PetscScalar*);
// SAVE
extern PetscErrorCode SAVE_1D_To_2D_ARRAY_REAL(const char*, PetscInt, PetscInt, PetscReal*);
extern PetscErrorCode SAVE_1D_To_2D_ARRAY_COMPLEX(const char*, const char*, PetscInt, PetscInt, PetscScalar*);

// IMPROVEMENTS
//extern PetscErrorCode SVL_2D_GRATING_AND_IMPROVEMENTS(PetscInt, PetscInt, PetscReal, PetscReal, PetscScalar*, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_2D_GRADING_VECTOR(PetscInt, PetscInt, PetscReal, PetscReal, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_2D_ELIMINATED_GRATING_ACCORD_THEIR_AMPLITUD(PetscInt, PetscInt, PetscScalar*, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_2D_IDENTIFIED_COLLINEAR_PLANAR_GRATING(PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*);
//extern PetscErrorCode SVL_2D_IMPLEMENT_IMPROVEMENTS(PetscInt, PetscInt, PetscInt*, PetscScalar*, PetscReal*, PetscReal*, PetscInt*, PetscInt*, PetscReal*, PetscReal*);

/*
extern PetscErrorCode SVL_ORIENTATION_FUNCTION(PetscInt, PetscInt, PetscReal, PetscReal, PetscReal*, PetscReal*);

extern PetscErrorCode SVL_ORIENTATION_VECTOR(PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal, PetscReal);
extern PetscErrorCode SVL_LATTICE_SPACING_FUNCTION (PetscInt, PetscInt, PetscReal*, PetscReal*);
extern PetscErrorCode SVL_Cart2Pol(PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
extern PetscErrorCode SVL_Rotation(PetscInt, PetscInt, PetscReal*, PetscReal*);
extern PetscErrorCode SVL_Spacing(PetscInt, PetscInt, PetscReal*, PetscReal*);
extern PetscErrorCode SVL_Pol2Cart(PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*, PetscReal*);
extern PetscErrorCode SVL_RHS(PetscInt, PetscInt, PetscReal*, PetscReal*, PetscReal*); 
extern PetscErrorCode SVL_fdder(PetscInt*, PetscInt*, PetscReal*, Mat);
extern PetscErrorCode SVL_LOOP(PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, Mat);
extern PetscErrorCode SVL_Print_Real_Array(const char*, PetscInt, PetscInt, PetscReal*);
extern PetscErrorCode SVL_Print_Complex_Array(PetscInt, PetscInt, PetscScalar*); 
*/
