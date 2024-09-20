
// Hybrid Solver Release 1.2
// Interface declarations
// Copyright (c) 2005-2007 Haifeng Qian

#if !defined(HybridSolver_P)
#define HybridSolver_P

#include <vector>

namespace hybridSolver {

//******************************************************************************
// create a matrix object.
// return a positive ID number to be used to access this matrix.
// return 0 if failed.
unsigned createMatrix( std::vector<unsigned>& rowIndex,
                       std::vector<unsigned>& colIndex,
                       std::vector<double>&   value,
                       bool                   isSymmetric = true );

//******************************************************************************
// free memory
// return 0 if sucessful, return 1 if failed
bool deleteMatrix( unsigned matrixID );

//******************************************************************************
// perform preconditioning with default quality
// overwrite previous preconditioner if any
// return 0 if sucessful, return 1 if failed
bool precondition( unsigned matrixID );

//******************************************************************************
// perform preconditioning with user-specified quality
// recommended "quality" value range (0.2,0.6)
// smaller value -> larger size and better accuracy
// overwrite previous preconditioner if any
// return 0 if sucessful, return 1 if failed
bool precondition( unsigned matrixID,
                   double   quality );

//******************************************************************************
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// convergence condition: norm(b-Ax)/norm(b) < "tolerance"
// if no preconditioner available, precondition(matrixID) is called internally
// return 0 if sucessful, return 1 if failed
bool solve( unsigned             matrixID,
            std::vector<double>& rhs,
            double               tolerance,
            std::vector<double>& solution );

//******************************************************************************
// return the current preconditioner for a symmetric matrix "matrixID".
// for numerical study purpose, not needed in normal usage.
// if no preconditioner available, precondition(matrixID) is called internally.
// "ordering" definition: original_index = ordering[ new_index ].
// "rowIndex" "colIndex" "value" describe an incomplete cholesky factor
// as an upper-triangular matrix in the new (i.e. after ordering) indices.
// return 0 if sucessful, return 1 if failed.
bool getPreconditioner( unsigned               matrixID,
                        std::vector<unsigned>& ordering,
                        std::vector<unsigned>& rowIndex, 
                        std::vector<unsigned>& colIndex,
                        std::vector<double>&   value );

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// return the "entryIndex"-th entry of the solution vector in "solution"
// accuracy: P[|error| < "delta"] > 99%
// donot call precondition(...) on the matrix prior to calling this function
// return 0 if sucessful, return 1 if failed
bool solveSingleEntry( unsigned             matrixID,
                       std::vector<double>& rhs,
                       double               delta,
                       unsigned             entryIndex,
                       double&              solution );

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// donot call precondition(...) on the matrix prior to calling this function
// return 0 if sucessful, return 1 if failed
bool genericRandomWalk( unsigned             matrixID,
                        std::vector<double>& rhs,
                        double               delta,
                        std::vector<double>& solution );

}

#endif
