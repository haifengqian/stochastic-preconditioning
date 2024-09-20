
#include "HybridSolver.h"
#include "hsUtility.h"

namespace hybridSolver {

//******************************************************************************
// create a matrix object.
// return a positive ID number to be used to access this matrix.
// return 0 if failed.
unsigned
createMatrix( std::vector<unsigned>& rowIndex,
              std::vector<unsigned>& colIndex,
              std::vector<double>&   value,
              bool                   isSymmetric )
{
    std::cout << "[Hybrid Solver]: Creating a ";
    if( isSymmetric ){
        std::cout << "symmetric";
    }
    else{
        std::cout << "asymmetric";
    }
    std::cout << " matrix.\n";

    if( ! isSymmetric ){
        std::cout << "[Hybrid Solver]: This version does not support "
                  << "asymmetric matrices.\n"
                  << "[Hybrid Solver]: Failed to create a matrix.\n";
        return 0;
    }
    unsigned size = rowIndex.size();
    if(( colIndex.size() != size )||( value.size() != size )){
        std::cout << "[Hybrid Solver]: Invalid input. Row/Col/Value vectrors "
                  << "must have the same size\n"
                  << "[Hybrid Solver]: Failed to create a matrix.\n";
        return 0;
    }
    if( !size ){
        std::cout << "[Hybrid Solver]: Empty matrix.\n"
                  << "[Hybrid Solver]: Failed to create a matrix.\n";
        return 0;
    }

    unsigned dimension1=0, dimension2=0;
    for( unsigned i=0; i<size; i++ ){
        if( dimension1 < rowIndex[i] ) dimension1 = rowIndex[i];
        if( dimension2 < colIndex[i] ) dimension2 = colIndex[i];
    }
    if( dimension1 != dimension2 ){
        std::cout << "[Hybrid Solver]: Invalid input. Not a square matrix.\n"
                  << "[Hybrid Solver]: Failed to create a matrix.\n";
        return 0;
    }
    dimension1++; 

    if( isSymmetric ){
        return createSymMatrix( rowIndex, colIndex, value, dimension1 );
    }
    else{
        return createAsymMatrix( rowIndex, colIndex, value, dimension1 );
    }
}

//******************************************************************************
// free memory
// return 0 if sucessful, return 1 if failed
bool
deleteMatrix( unsigned matrixID )
{
    return freeMatrix( matrixID );
}

//******************************************************************************
// perform preconditioning with default quality
// overwrite previous preconditioner if any
// return 0 if sucessful, return 1 if failed
bool
precondition( unsigned matrixID )
{
    return precond( matrixID );
}

//******************************************************************************
// perform preconditioning with user-specified quality
// recommended "quality" value range (0.2,0.6)
// smaller value -> larger size and better accuracy
// overwrite previous preconditioner if any
// return 0 if sucessful, return 1 if failed
bool
precondition( unsigned matrixID,
              double   quality )
{
    return precond( matrixID, quality );
}

//******************************************************************************
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// convergence condition: norm(b-Ax)/norm(b) < "tolerance"
// if no conditioner available, precondition(matrixID) is called internally
// return 0 if sucessful, return 1 if failed
bool
solve( unsigned             matrixID,
       std::vector<double>& rhs,
       double               tolerance,
       std::vector<double>& solution )
{
    return getSolution( matrixID, rhs, tolerance, solution );
}

//******************************************************************************
// return the "entryIndex"-th entry of the solution vector in "solution"
// accuracy: P[|error| < "delta"] > 99%
// return 0 if sucessful, return 1 if failed
bool
solveSingleEntry( unsigned             matrixID,
                  std::vector<double>& rhs,
                  double               delta,
                  unsigned             entryIndex,
                  double&              solution )
{
    return singleEntry( matrixID, rhs, delta, entryIndex, solution );
}

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// return 0 if sucessful, return 1 if failed
bool
genericRandomWalk( unsigned             matrixID,
                   std::vector<double>& rhs,
                   double               delta,
                   std::vector<double>& solution )
{
    return generic( matrixID, rhs, delta, solution );
}

//******************************************************************************
// return the current preconditioner for a symmetric matrix "matrixID".
// for numerical study purpose, not needed in normal usage.
// if no preconditioner available, precondition(matrixID) is called internally.
// "ordering" definition: original_index = ordering[ new_index ].
// "rowIndex" "colIndex" "value" describe an incomplete cholesky factor
// as an upper-triangular matrix in the new (i.e. after ordering) indices.
// return 0 if sucessful, return 1 if failed.
bool
getPreconditioner( unsigned               matrixID,
                   std::vector<unsigned>& ordering,
                   std::vector<unsigned>& rowIndex, 
                   std::vector<unsigned>& colIndex,
                   std::vector<double>&   value )
{
    return writePreconditioner(matrixID,ordering,rowIndex,colIndex,value);
}

}
