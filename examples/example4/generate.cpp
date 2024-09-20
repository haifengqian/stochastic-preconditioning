
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

void
generateMatrix( unsigned dimension, char* filename )
{
    std::vector<std::vector<std::vector<unsigned> > > grid;
    grid.resize(dimension);
	unsigned i;
    for(i=0;i<dimension;i++){
        grid[i].resize(dimension);
        for(unsigned ii=0;ii<dimension;ii++){
            grid[i][ii].resize(dimension);
            for(unsigned iii=0;iii<dimension;iii++){
                grid[i][ii][iii] = i * dimension * dimension + ii * dimension + iii;
            }
        }
    }
    FILE *inFile = fopen (filename, "w");
    fprintf( inFile, "%d\n", dimension * dimension * dimension );
    fprintf( inFile, "%d\n", 6 * (dimension-1) * dimension * dimension + dimension * dimension * dimension);
    for(i=0;i<dimension;i++){
        for(unsigned ii=0;ii<dimension;ii++){
            for(unsigned iii=0;iii<dimension;iii++){
                if( i > 0 ){
                    fprintf( inFile, "%d %d -1\n", grid[i][ii][iii], grid[i-1][ii][iii]);
                }
                if( i < dimension-1 ){
                    fprintf( inFile, "%d %d -1\n", grid[i][ii][iii], grid[i+1][ii][iii]);
                }
                if( ii > 0 ){
                    fprintf( inFile, "%d %d -1\n", grid[i][ii][iii], grid[i][ii-1][iii]);
                }
                if( ii < dimension-1 ){
                    fprintf( inFile, "%d %d -1\n", grid[i][ii][iii], grid[i][ii+1][iii]);
                }
                if( iii > 0 ){
                    fprintf( inFile, "%d %d -1\n", grid[i][ii][iii], grid[i][ii][iii-1]);
                }
                if( iii < dimension-1 ){
                    fprintf( inFile, "%d %d -1\n", grid[i][ii][iii], grid[i][ii][iii+1]);
                }
                fprintf( inFile, "%d %d 6\n", grid[i][ii][iii], grid[i][ii][iii] );
            }
        }
    }
    for(i=0;i<dimension*dimension*dimension;i++) fprintf( inFile, "1.0\n" ); 
    fclose(inFile);
}


int
main(int argv, char *argc[])
{
    if (argv!=3) {
		printf ("ERROR: Wrong number of arguments!\nUsage:\n\t  %s <dimension> <FileName>\n", argc[0]);
		return 1;
	}
	unsigned dimension;
	sscanf( argc[1], "%d", &dimension );
    generateMatrix( dimension, argc[2] );
	return 0;
}
