
#include <stdio.h>
#include <iostream>
#include <string>
#include "HybridSolver.h"

void
generateMatrix( unsigned dimension, std::string filename )
{
    std::vector<std::vector<std::vector<unsigned> > > grid;
    grid.resize(10);
    for(unsigned i=0;i<10;i++){
        grid[i].resize(dimension);
        for(unsigned ii=0;ii<dimension;ii++){
            grid[i][ii].resize(dimension);
            for(unsigned iii=0;iii<dimension;iii++){
                grid[i][ii][iii] = i * dimension * dimension + ii * dimension + iii;
            }
        }
    }
    FILE *inFile = fopen (filename.c_str(), "w");
    fprintf( inFile, "%d\n", 10 * dimension * dimension );
    fprintf( inFile, "%d\n", 40 * (dimension-1) * dimension + 18 * dimension * dimension + 10 * dimension * dimension);
    for(unsigned i=0;i<10;i++){
        for(unsigned ii=0;ii<dimension;ii++){
            for(unsigned iii=0;iii<dimension;iii++){
                unsigned counter = 0;
                if( i > 0 ){
                    fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i-1][ii][iii], -1.0 );
                    counter++;
                }
                if( i < 9 ){
                    fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i+1][ii][iii], -1.0 );
                    counter++;
                }
                if( ii > 0 ){
                    fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i][ii-1][iii], -1.0 );
                    counter++;
                }
                if( ii < dimension-1 ){
                    fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i][ii+1][iii], -1.0 );
                    counter++;
                }
                if( iii > 0 ){
                    fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i][ii][iii-1], -1.0 );
                    counter++;
                }
                if( iii < dimension-1 ){
                    fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i][ii][iii+1], -1.0 );
                    counter++;
                }
                if(( i == 0 )&&( (ii/10)*10 == ii )&&( (iii/10)*10 == iii )){
                    counter++;
                }
                fprintf( inFile, "%d %d %f\n", grid[i][ii][iii], grid[i][ii][iii], (float)counter );
            }
        }
    }
    for(unsigned i=0;i<10*dimension*dimension;i++) fprintf( inFile, "1.0\n" ); 
    fclose(inFile);
}

void
generateMatrixFiles()
{
    // testcase 1
    // 10x10x10 3D uniform grid, 1 corner node diagonally dominant, unit rhs
    generateMatrix( 10, "test1.matrix" );
    std::cout << "Test case 1 written in file test1.matrix\n\n";
    
    // testcase 2
    // 50x50x10 3D uniform grid, 25 top-layer nodes diagonally dominant, unit rhs
    generateMatrix( 50, "test2.matrix" );
    std::cout << "Test case 2 written in file test2.matrix\n\n";

    // testcase 3
    // 100x100x10 3D uniform grid, 100 top-layer nodes diagonally dominant, unit rhs
    generateMatrix( 100, "test3.matrix" );
    std::cout << "Test case 3 written in file test3.matrix\n\n";
}

void
solveTestcase( std::string filename )
{
    std::vector<unsigned> rowIndex, colIndex;
    std::vector<double> value, rhs, solution;
    rowIndex.clear();
    colIndex.clear();
    value.clear();

    // read matrix file
	FILE *inFile;
	char inputbuf[256];
	int i,dimension, entrynumber, fromnodeindex,tonodeindex;
	float devicevalue;
	inFile = fopen(filename.c_str(), "r");
	fgets(inputbuf,255,inFile);
	sscanf(inputbuf,"%d",&dimension);
	fgets(inputbuf,255,inFile);
	sscanf(inputbuf,"%d",&entrynumber);
	for(i=0;i<entrynumber;i++){
		fgets(inputbuf,255,inFile);
		sscanf(inputbuf,"%d %d %f",&fromnodeindex,&tonodeindex,&devicevalue);
		if(fromnodeindex <= tonodeindex){
            rowIndex.push_back( fromnodeindex );
            colIndex.push_back( tonodeindex );
            value.push_back( devicevalue );
		}
	}
	rhs.resize( dimension );
	for(i=0;i<dimension;i++){
		fgets(inputbuf,255,inFile);
		sscanf(inputbuf,"%f",&devicevalue);
		rhs[i] = devicevalue;
	}
	fclose(inFile);
    // end of read matrix file

    // solver package routines
    unsigned matrixID = hybridSolver::createMatrix( rowIndex,
                                                     colIndex,
                                                     value,
                                                     true );
    if( ! matrixID ) return;
    bool returnFlag = hybridSolver::precondition( matrixID );
    if( returnFlag ) return;
    double tolerance = 0.000001;
    returnFlag = hybridSolver::solve( matrixID,
                                      rhs,
                                      tolerance,
                                      solution );    
    if( returnFlag ) return;
    // end of solver package routines

    // write solution into output file
    std::string name = filename + ".solution";
    inFile = fopen (name.c_str(), "w");
    for( i=0;i<dimension;i++ ){
        fprintf(inFile, "%d %g\n", i, solution[i]);
    }
    fclose(inFile);
    // end of write solution into output file

    std::cout << "\nTest case solved, result written in file " << name <<"\n\n";

    // free memory
    returnFlag = hybridSolver::deleteMatrix( matrixID );
}

int
main()
{
    std::cout << "This is a testing program that generates three matrices and then solve them.\n\n";
    generateMatrixFiles();
    solveTestcase( "test1.matrix" );
    solveTestcase( "test2.matrix" );
    solveTestcase( "test3.matrix" );
}
