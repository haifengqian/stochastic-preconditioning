
#include <iostream>
#include "HybridSolver.h"

int
main()
{
    std::vector<unsigned> rowIndex, colIndex;
    std::vector<double> value, rhs, solution;
    rowIndex.clear();
    colIndex.clear();
    value.clear();
    
    rowIndex.push_back( 0 );
    colIndex.push_back( 0 );
    value.push_back( 1.5 );

    rowIndex.push_back( 0 );
    colIndex.push_back( 2 );
    value.push_back( -1 );
    
    rowIndex.push_back( 1 );
    colIndex.push_back( 1 );
    value.push_back( 2 );
    
    rowIndex.push_back( 1 );
    colIndex.push_back( 2 );
    value.push_back( -1 );
    
    rowIndex.push_back( 2 );
    colIndex.push_back( 2 );
    value.push_back( 2.25 );
    
    rowIndex.push_back( 2 );
    colIndex.push_back( 3 );
    value.push_back( -0.25 );

    rowIndex.push_back( 3 );
    colIndex.push_back( 3 );
    value.push_back( 1.25 );
    
    unsigned matrixID = hybridSolver::createMatrix( rowIndex,
                                                    colIndex,
                                                    value,
                                                    true );

    if( ! matrixID ) return 1;

    bool returnFlag = hybridSolver::precondition( matrixID );
    
    if( returnFlag ) return 1;

    rhs.resize( 4 );
    rhs[ 0 ] = 0.2;
    rhs[ 1 ] = 0.9;
    rhs[ 2 ] = -0.05;
    rhs[ 3 ] = 0.95;

    double tolerance = 0.000001;

    returnFlag = hybridSolver::solve( matrixID,
                                      rhs,
                                      tolerance,
                                      solution );

    if( returnFlag ) return 1;
           
    std::cout << "Soluton#1: " << solution[0] << " " << solution[1] << " "
              << solution[2] << " " << solution[3] << "\n";
              
    rhs[ 0 ] = -0.5;
    rhs[ 1 ] = 0;
    rhs[ 2 ] = 2.25;
    rhs[ 3 ] = 0.75;
    
    returnFlag = hybridSolver::solve( matrixID,
                                      rhs,
                                      tolerance,
                                      solution );

    if( returnFlag ) return 1;
           
    std::cout << "Soluton#2: " << solution[0] << " " << solution[1] << " "
              << solution[2] << " " << solution[3] << "\n";

    hybridSolver::deleteMatrix( matrixID );

    return 0;
}
