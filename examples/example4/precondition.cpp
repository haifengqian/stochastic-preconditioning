
#include <stdio.h>
#include <iostream>
#include <string>
#include "HybridSolver.h"

unsigned
readMatrixFile( char* filename )
{
    std::vector<unsigned> rowIndex, colIndex;
    std::vector<double> value;
    rowIndex.clear();
    colIndex.clear();
    value.clear();

    // read matrix file
	FILE *inFile;
	char inputbuf[256];
	int i,dimension, entrynumber, fromnodeindex,tonodeindex;
	float devicevalue;
	inFile = fopen(filename, "r");
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
	fclose(inFile);
    // end of read matrix file

    return hybridSolver::createMatrix( rowIndex,
                                       colIndex,
                                       value,
                                       true );
}

void
writeFile( std::vector<unsigned>& values,
		  std::string& filename )
{
	FILE* transfer = fopen( filename.c_str(), "wb" );
	unsigned size = values.size();
	fwrite( &size, sizeof(unsigned), 1, transfer );
	for(unsigned i=0; i<size; i++){
		unsigned number = values[i];
		fwrite( &number, sizeof(unsigned), 1, transfer );
	}
	fclose( transfer );
}

void
writeFile( unsigned dimension,
		   std::vector<unsigned>& rowIndex,
		   std::vector<unsigned>& colIndex,
		   std::vector<double>& values,
		   std::string& filename )
{
	FILE* transfer = fopen( filename.c_str(), "wb" );
	unsigned numnonzero = rowIndex.size();
	fwrite( &dimension, sizeof(unsigned), 1, transfer );
	fwrite( &numnonzero, sizeof(unsigned), 1, transfer );
	unsigned i;
	for(i=0; i<numnonzero; i++){
		unsigned row = rowIndex[i];
		fwrite( &row, sizeof(unsigned), 1, transfer );
	}
	for(i=0; i<numnonzero; i++){
		unsigned col = colIndex[i];
		fwrite( &col, sizeof(unsigned), 1, transfer );
	}
	for(i=0; i<numnonzero; i++){
		double value = values[i];
		fwrite( &value, sizeof(double), 1, transfer );
	}
	fclose( transfer );
}

void
writePreconditioner( char* filename, double delta )
{
	unsigned matrixID = readMatrixFile( filename );
	if( !matrixID ) return;

    bool returnFlag = hybridSolver::precondition( matrixID, delta );
    if( returnFlag ) return;

	std::vector<unsigned> ordering, rowIndex, colIndex;
	std::vector<double>   values;
    returnFlag = hybridSolver::getPreconditioner( matrixID,
												  ordering,
												  rowIndex, 
												  colIndex,
												  values );
    if( returnFlag ) return;

    // write preconditioner into output file
    std::string name( filename );
	unsigned start = 0;
	for(unsigned i=0; i<name.size(); i++){
		if( name[i] == '\\' ) start = i+1;
	}
	name = name.substr( start, name.size()-start );
	std::string orderingFile = name + ".ordering";
	writeFile( ordering,
			   orderingFile );
	std::string choleskyFile = name + ".cholesky";
	writeFile( ordering.size(),
		       rowIndex,
               colIndex,
			   values,
			   choleskyFile );

    // free memory
    returnFlag = hybridSolver::deleteMatrix( matrixID );
}

int
main(int argv, char *argc[])
{
    if (argv!=3) {
		printf ("ERROR: Wrong number of arguments!\nUsage:\n\t  %s <filename> <delta>\n", argc[0]);
		return 1;
	}
	float temp;
	unsigned returnvalue = sscanf( argc[2], "%f", &temp);
	if( returnvalue != 1 ) return 1;
    double delta = temp;
	writePreconditioner( argc[1], delta );
	return 0;
}
