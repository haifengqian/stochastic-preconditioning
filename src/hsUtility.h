
#if !defined(hsUtility_P)
#define hsUtility_P

#include <iostream>
#include <set>
#include <vector>
#include <assert.h>
#include "hsSymRmatCore.h"

#define DEFAULTSTEPDELTA 0.4
#define VERYLITTLE 0.001
#define VERYVERYLITTLE 0.000001 

namespace hybridSolver {

enum MatrixStatus{
    hscNonExistent              = 0,
    hscSymRmatRaw               = 1, // R-matrices defined in Qian ICCAD05
    hscSymRmatPreconditioned    = 2,
    hscSymSpdRaw                = 3, // Not R-matrices, but sym. posi. definite
    hscSymSpdPreconditioned     = 4,
    hscSymIndefRaw              = 5, // sym. indefinite matrices
    hscSymIndefPreconditioned   = 6,
    hscAsymPdRaw                = 7, // asymmetric posi. definite
    hscAsymPdPreconditioned     = 8,
    hscAsymIndefRaw             = 9, // asymmetric indefinite
    hscAsymIndefPreconditioned  = 10,
};

enum ErrorCode{
    hscCannotMarkNonsinkPath = 0,
    hscCannotFindPathEnd     = 1,
    hscCannotCutEdge         = 2,
    hscCannotRecoversinkPath = 3,
    hscCannotRecoverPath     = 4,
    hscZeroDiagonal          = 5,
    hscFailedGetOrdering     = 6,
};

class Error {
 public:
    Error( ErrorCode code ){
        _code = code;
    }
    ErrorCode getCode(){
        return _code;
    }
    void showMessage(){
        switch( _code ){
            case hscZeroDiagonal:
                std::cout << " Numerical error. Zero digonal component encountered. ";
                return;
            case hscCannotMarkNonsinkPath:
            case hscCannotFindPathEnd:
            case hscCannotCutEdge:
            case hscCannotRecoversinkPath:
            case hscCannotRecoverPath:
            case hscFailedGetOrdering:
                std::cout << " Internal error. Code " << _code << ". Please send bug report to Haifeng Qian. ";
                return;
            default:
                std::cout << " Unrecognized internal error. Please send bug report to Haifeng Qian. ";
        }
    }
 private:
    ErrorCode _code;
};

//******************************************************************************
// create a symmetric matrix object.
// return a positive ID number to be used to access this matrix.
// return 0 if failed.
unsigned createSymMatrix( std::vector<unsigned>& rowIndex,
                          std::vector<unsigned>& colIndex,
                          std::vector<double>&   value,
                          unsigned               dimension );

//******************************************************************************
// create a matrix object.
// return a positive ID number to be used to access this matrix.
// return 0 if failed.
unsigned createAsymMatrix( std::vector<unsigned>& rowIndex,
                           std::vector<unsigned>& colIndex,
                           std::vector<double>&   value,
                           unsigned               dimension );

//******************************************************************************
// free memory
// return 0 if sucessful, return 1 if failed
bool freeMatrix( unsigned matrixID );

//******************************************************************************
// free memory
void freeSymRawMatrix( struct global* globalpointer );

//******************************************************************************
// preconditioning
// return 0 if sucessful, return 1 if failed
bool precond( unsigned matrixID,
              double   quality = DEFAULTSTEPDELTA );

//******************************************************************************
// solve
bool getSolution( unsigned             matrixID,
                  std::vector<double>& rhs,
                  double               tolerance,
                  std::vector<double>& solution );

//******************************************************************************
// write a RHS into the global data structure
void writeRHS( struct global* pointer,
               std::vector<double>& rhs );

//******************************************************************************
// copy soluton from the global data structure to a vector
// free memory related to this solve in the global data structure
void getSolutionFreeMemory( struct global* pointer,
                            std::vector<double>& solution );

//******************************************************************************
// cleanup when something goes wrong in the preconditioning stage
// delete matrix
void showErrorCleanUp( Error& error, unsigned matrixID );

//******************************************************************************
// cleanup when something goes wrong in the solving stage
// delete matrix
void showErrorCleanUp2( Error& error, unsigned matrixID );

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
// return 0 if sucessful, return 1 if failed
bool generic( unsigned             matrixID,
              std::vector<double>& rhs,
              double               delta,
              std::vector<double>& solution );

//******************************************************************************
// return the "entryIndex"-th entry of the solution vector in "solution"
// accuracy: P[|error| < "delta"] > 99%
// return 0 if sucessful, return 1 if failed
bool singleEntry( unsigned             matrixID,
                  std::vector<double>& rhs,
                  double               delta,
                  unsigned             entryIndex,
                  double&              solution );

//******************************************************************************
// return the current preconditioner for a symmetric matrix "matrixID".
// for numerical study purpose, not needed in normal usage.
// if no preconditioner available, precondition(matrixID) is called internally.
// "ordering" definition: original_index = ordering[ new_index ].
// "rowIndex" "colIndex" "value" describe an incomplete cholesky factor
// as an upper-triangular matrix in the new (i.e. after ordering) indices.
// return 0 if sucessful, return 1 if failed.
bool writePreconditioner( unsigned               matrixID,
                          std::vector<unsigned>& ordering,
                          std::vector<unsigned>& rowIndex, 
                          std::vector<unsigned>& colIndex,
                          std::vector<double>&   value );

}

#endif
