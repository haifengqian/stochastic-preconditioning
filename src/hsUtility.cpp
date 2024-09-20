
#include "hsUtility.h"
#include "hsGenericRW.h"

namespace hybridSolver {

static std::vector<void*>        matrices;
static std::vector<MatrixStatus> status;

//******************************************************************************
// create a symmetric matrix object.
// return a positive ID number to be used to access this matrix.
// return 0 if failed.
unsigned
createSymMatrix( std::vector<unsigned>& rowIndex,
                 std::vector<unsigned>& colIndex,
                 std::vector<double>&   value,
                 unsigned               dimension )
{
    if( matrices.empty() ){ 
        matrices.resize( 1 );
        status.resize( 1 );
        matrices[ 0 ] = NULL;
        status[ 0 ] = hscNonExistent;
    }
    
	struct edge** adjmatrix=(struct edge **)calloc(dimension,sizeof(struct edge *));
	int* nodedegree=(int *)calloc(dimension,sizeof(int));
	double* capacitance=(double *)calloc(dimension,sizeof(double));
	struct edge* homeregistry=NULL;

	unsigned i; 
	for( i=0; i<rowIndex.size(); i++ ){
        unsigned fromnodeindex = rowIndex[ i ];
        unsigned tonodeindex = colIndex[ i ];
        double devicevalue = value[ i ];
		if(fromnodeindex==tonodeindex){
            if( devicevalue <=0 ){  
                std::cout << "[Hybrid Solver]: This version only supports "
                          << "matrices with positive diagonal entries.\n"
    	                  << "[Hybrid Solver]: Failed to create a matrix.\n";
    	        for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
	    		free( adjmatrix );
		    	free( nodedegree );
			    free( capacitance );
			    return 0;
            }
		    if( capacitance[fromnodeindex] != 0 ){  
		        std::cout << "[Hybrid Solver]: Duplicate diagonal entries.\n"
    	                  << "[Hybrid Solver]: Failed to create a matrix.\n";
    	        for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
	    		free( adjmatrix );
		    	free( nodedegree );
			    free( capacitance );
			    return 0;
		    }
		    capacitance[fromnodeindex]=devicevalue; 
		}
		else if( fromnodeindex > tonodeindex ){  
		    std::cout << "[Hybrid Solver]: When describing a symmetric matrix, "
                      << "only describe the upper triangular part, and "
                      << "row index cannot be greater than column index.\n"
                      << "[Hybrid Solver]: Failed to create a matrix.\n";
			for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
			free( adjmatrix );
			free( nodedegree );
			free( capacitance );
			return 0;
		}
		else{                         
			if(devicevalue > 0){      
                std::cout << "[Hybrid Solver]: This version does not support "
                          << "positive off-diagonal entries.\n"
                          << "[Hybrid Solver]: Failed to create a matrix.\n";
				for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
				free( adjmatrix );
				free( nodedegree );
				free( capacitance );
				return 0;
			}
			if( devicevalue == 0 ){  
			    continue;
			}
			
			struct edge* edgepointer=adjmatrix[fromnodeindex];
			adjmatrix[fromnodeindex]=(struct edge *)malloc(1*sizeof(struct edge));
			adjmatrix[fromnodeindex]->nodeindex=tonodeindex;
			adjmatrix[fromnodeindex]->conductance=-devicevalue;
			adjmatrix[fromnodeindex]->next=edgepointer;
			nodedegree[fromnodeindex]++;
			
			edgepointer=adjmatrix[tonodeindex];
			adjmatrix[tonodeindex]=(struct edge *)malloc(1*sizeof(struct edge));
			adjmatrix[tonodeindex]->nodeindex=fromnodeindex;
			adjmatrix[tonodeindex]->conductance=-devicevalue;
			adjmatrix[tonodeindex]->next=edgepointer;
			nodedegree[tonodeindex]++;
		}
	}

	unsigned realnodenum = dimension;
	for( i=0; i<dimension; i++ ){
	    if( capacitance[i] <= 0 ){
	        if(( capacitance[i] == 0 )&&( nodedegree[i] == 0 )){
	            std::cout << "[Hybrid Solver]: Empty row.\n";
	        }
	        else{
	            std::cout << "[Hybrid Solver]: This version only supports ";
                std::cout << "matrices with positive diagonal entries.\n";
	        }
	        std::cout << "[Hybrid Solver]: Failed to create a matrix.\n";
	        for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
			free( adjmatrix );
			free( nodedegree );
			free( capacitance );
			freeedgepointer( homeregistry );
			return 0;
	    }
        double temp2 = 0;
		struct edge* edgepointer=adjmatrix[i];
		std::set<int> checkFlag;
		checkFlag.clear();
		for(int j=0;j<nodedegree[i];j++){
		    if( checkFlag.find( edgepointer->nodeindex ) != checkFlag.end() ){ 
		        std::cout << "[Hybrid Solver]: Duplicate off-diagonal entries.\n"
                          << "[Hybrid Solver]: Failed to create a matrix.\n";
			    for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
			    free( adjmatrix );
			    free( nodedegree );
			    free( capacitance );
			    freeedgepointer( homeregistry );
			    return 0;
		    }
            checkFlag.insert( edgepointer->nodeindex );
			temp2=temp2+edgepointer->conductance;
			edgepointer=edgepointer->next;
		}
		if ((capacitance[i]-temp2)/capacitance[i] < -VERYLITTLE){ 
            std::cout << "[Hybrid Solver]: This version only supports "
                      << "matrices that are irreducibly diagonally dominant.\n"
                      << "[Hybrid Solver]: Failed to create a matrix.\n";
			for(unsigned ii=0;ii<dimension;ii++) freeedgepointer( adjmatrix[ii] );
			free( adjmatrix );
			free( nodedegree );
			free( capacitance );
			freeedgepointer( homeregistry );
			return 0;
		}
		if ((capacitance[i]-temp2)/capacitance[i] > VERYLITTLE){ 
			edgepointer=adjmatrix[i];
			adjmatrix[i]=(struct edge *)malloc(1*sizeof(struct edge));
			adjmatrix[i]->nodeindex=realnodenum;
			adjmatrix[i]->conductance=capacitance[i]-temp2;
			adjmatrix[i]->next=edgepointer;
			nodedegree[i]++;
			edgepointer = homeregistry;
			homeregistry = (struct edge *)malloc(1*sizeof(struct edge));
			homeregistry->nodeindex = i;
			homeregistry->probability = realnodenum; 
			homeregistry->conductance=capacitance[i]-temp2;
			homeregistry->next=edgepointer;
			realnodenum++;
		}
	}

    struct global* globalpointer = (struct global *)calloc(1,sizeof(struct global));
	globalpointer->dimension = dimension;
	globalpointer->realnodenum = realnodenum;
	globalpointer->adjmatrix = adjmatrix;
	globalpointer->nodedegree = nodedegree;
	globalpointer->diagonal = capacitance;
	globalpointer->homeregistry = homeregistry;
    
    matrices.push_back( static_cast<void*>( globalpointer ) );
    status.push_back( hscSymRmatRaw );
    
    std::cout << "[Hybrid Solver]: A symmetric matrix created with dimension "
              << dimension << " and " << rowIndex.size()*2-dimension
              << " entries. ID " << matrices.size() - 1 << ".\n";
    return matrices.size() - 1;
}

//******************************************************************************
// create a matrix object.
// return a positive ID number to be used to access this matrix.
// return 0 if failed.
unsigned
createAsymMatrix( std::vector<unsigned>& rowIndex,
                  std::vector<unsigned>& colIndex,
                  std::vector<double>&   value,
                  unsigned               dimension )
{
    std::cout << "[Hybrid Solver]: This version does not support ";
    std::cout << "asymmetric matrices.\n";
    std::cout << "[Hybrid Solver]: Failed to create a matrix.\n";
    return 0;
}

//******************************************************************************
// free memory
// return 0 if sucessful, return 1 if failed
bool
freeMatrix( unsigned matrixID )
{
    if( matrixID >= matrices.size() ){
        std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
        return 1;
    }

    switch( status[ matrixID ] ){
        case hscNonExistent:
            std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
            return 1;
        case hscSymRmatRaw:
            freeSymRawMatrix( static_cast<struct global*>( matrices[ matrixID ] ) );
            status[ matrixID ] = hscNonExistent;
            std::cout << "[Hybrid Solver]: Matrix " << matrixID << " deleted.\n";
            return 0;
        case hscSymRmatPreconditioned:
            freespace( static_cast<struct global*>( matrices[ matrixID ] ) );
            status[ matrixID ] = hscNonExistent;
            std::cout << "[Hybrid Solver]: Matrix " << matrixID << " deleted.\n";
            return 0;
        case hscSymSpdRaw:
        case hscSymSpdPreconditioned:
        case hscSymIndefRaw:
        case hscSymIndefPreconditioned:
        case hscAsymPdRaw:
        case hscAsymPdPreconditioned:
        case hscAsymIndefRaw:
        case hscAsymIndefPreconditioned:
        default:
            std::cout << "[Hybrid Solver]: Internal error.\n";
            std::cout << "[Hybrid Solver]: Failed to delete matrix.\n";
            return 1;
    }
}

//******************************************************************************
// free memory
void
freeSymRawMatrix( struct global* globalpointer )
{
	for(int i=0;i<globalpointer->dimension;i++) freeedgepointer(globalpointer->adjmatrix[i]);
	free(globalpointer->adjmatrix);
	free(globalpointer->diagonal);
	freeedgepointer(globalpointer->homeregistry);
	free(globalpointer->nodedegree);
	free(globalpointer);
}


//******************************************************************************
// preconditioning
// return 0 if sucessful, return 1 if failed
bool
precond( unsigned matrixID,
         double   quality )
{
    std::cout << "[Hybrid Solver]: Preconditioning matrix " << matrixID << ".\n";
    double startime= seconds();
    
    if( matrixID >= matrices.size() ){
        std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
        return 1;
    }

    try{
        switch( status[ matrixID ] ){
            case hscNonExistent:
                std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
                return 1;
            case hscSymRmatRaw:
                shrinkGraph( static_cast<struct global*>( matrices[ matrixID ] ) );
                randomwalk(  static_cast<struct global*>( matrices[ matrixID ] ),
                             quality );
                status[ matrixID ] = hscSymRmatPreconditioned;
                break;
            case hscSymRmatPreconditioned:
                std::cout << "[Hybrid Solver]: This version does not support "
                          << "preconditioning more than once on a matrix. "
                          << "Please create another matrix and call again.\n"
                          << "[Hybrid Solver]: Preconditioning aborted.\n";
                return 1;
            case hscSymSpdRaw:
            case hscSymSpdPreconditioned:
            case hscSymIndefRaw:
            case hscSymIndefPreconditioned:
            case hscAsymPdRaw:
            case hscAsymPdPreconditioned:
            case hscAsymIndefRaw:
            case hscAsymIndefPreconditioned:
            default:
                std::cout << "[Hybrid Solver]: Internal error.\n"
                          << "[Hybrid Solver]: Failed to run preconditioning.\n";
                return 1;
        }
    }
    catch( Error& error ){
        showErrorCleanUp( error, matrixID );
        return 1;
    }

    std::cout << "[Hybrid Solver]: Done preconditioning matrix " << matrixID << ".\n"
              << "[Hybrid Solver]: Preconditioning runtime " << seconds()-startime << " seconds.\n";
    return 0;
}

//******************************************************************************
// cleanup when something goes wrong in the preconditioning stage
// delete matrix
void
showErrorCleanUp( Error& error, unsigned matrixID )
{
    std::cout << "[Hybrid Solver]:";
    error.showMessage();
    std::cout << "\n"
              << "[Hybrid Solver]: Failed to precondition matrix " << matrixID << ".\n"
              << "[Hybrid Solver]: Attempting to delete matrix " << matrixID << ".\n";

    struct global* pointer;
    switch( status[ matrixID ] ){
        case hscSymRmatRaw:
            pointer = static_cast<struct global*>( matrices[ matrixID ] );
            switch( error.getCode() ){
                case hscCannotMarkNonsinkPath:
                case hscCannotFindPathEnd:
                case hscCannotCutEdge:
                    freeSymRawMatrix( pointer );
                    status[ matrixID ] = hscNonExistent;
                    std::cout << "[Hybrid Solver]: Matrix " << matrixID << " deleted.\n";
                    return;
                case hscCannotRecoversinkPath:
                case hscCannotRecoverPath:
                case hscZeroDiagonal:
                default:
                    std::cout << "[Hybrid Solver]: Internal error.\n"
                              << "[Hybrid Solver]: Failed to free memory from matrix " << matrixID << ".\n";
                    status[ matrixID ] = hscNonExistent;
                    return;
            }
        case hscSymSpdRaw:
        case hscNonExistent:
        case hscSymRmatPreconditioned:
        case hscSymSpdPreconditioned:
        case hscSymIndefRaw:
        case hscSymIndefPreconditioned:
        case hscAsymPdRaw:
        case hscAsymPdPreconditioned:
        case hscAsymIndefRaw:
        case hscAsymIndefPreconditioned:
        default:
            std::cout << "[Hybrid Solver]: Internal error.\n"
                      << "[Hybrid Solver]: Failed to free memory from matrix " << matrixID << ".\n";
            status[ matrixID ] = hscNonExistent;
            return;
    }
}

//******************************************************************************
// solve
// return 0 if sucessful, return 1 if failed
bool
getSolution( unsigned             matrixID,
             std::vector<double>& rhs,
             double               tolerance,
             std::vector<double>& solution )
{
    std::cout << "[Hybrid Solver]: Solving matrix " << matrixID << " with a given RHS.\n";
    double startime= seconds();
    
    if( matrixID >= matrices.size() ){
        std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
        return 1;
    }

    struct global* pointer;

    try{
        switch( status[ matrixID ] ){
            case hscNonExistent:
                std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
                return 1;
            case hscSymRmatRaw:
                if( precond( matrixID ) ){
                    return 1;
                }
            case hscSymRmatPreconditioned:
                pointer = static_cast<struct global*>( matrices[ matrixID ] );
                if( rhs.size() != static_cast<unsigned>( pointer->dimension ) ){
                    std::cout << "[Hybrid Solver]: Invalid RHS. Solving aborted.\n";
                    return 1;
                }
                writeRHS( pointer, rhs );
                solve( pointer, tolerance );
                getSolutionFreeMemory( pointer, solution );
                break;
            case hscSymSpdRaw:
            case hscSymSpdPreconditioned:
            case hscSymIndefRaw:
            case hscSymIndefPreconditioned:
            case hscAsymPdRaw:
            case hscAsymPdPreconditioned:
            case hscAsymIndefRaw:
            case hscAsymIndefPreconditioned:
            default:
                std::cout << "[Hybrid Solver]: Internal error.\n"
                          << "[Hybrid Solver]: Failed solving.\n";
                return 1;
        }
    }
    catch( Error& error ){
        showErrorCleanUp2( error, matrixID );
        return 1;
    }

    std::cout << "[Hybrid Solver]: Done solving matrix " << matrixID << " with the given RHS.\n"
              << "[Hybrid Solver]: Solving runtime " << seconds()-startime << " seconds.\n";
    return 0;
}

//******************************************************************************
// cleanup when something goes wrong in the solving stage
// delete matrix
void
showErrorCleanUp2( Error& error, unsigned matrixID )
{
    std::cout << "[Hybrid Solver]:";
    error.showMessage();
    std::cout << "\n"
              << "[Hybrid Solver]: Failed to solve matrix " << matrixID << " with the given RHS.\n"
              << "[Hybrid Solver]: Attempting to delete matrix " << matrixID << ".\n";

    struct global* pointer;
    switch( status[ matrixID ] ){
        case hscSymRmatPreconditioned:
            pointer = static_cast<struct global*>( matrices[ matrixID ] );
            switch( error.getCode() ){
                case hscCannotRecoversinkPath:
                case hscCannotRecoverPath:
                    freespace( pointer );
                    status[ matrixID ] = hscNonExistent;
                    std::cout << "[Hybrid Solver]: Matrix " << matrixID << " deleted.\n";
                    return;
                case hscCannotMarkNonsinkPath:
                case hscCannotFindPathEnd:
                case hscCannotCutEdge:
                case hscZeroDiagonal:
                default:
                    std::cout << "[Hybrid Solver]: Internal error.\n"
                              << "[Hybrid Solver]: Failed to free memory from matrix " << matrixID << ".\n";
                    status[ matrixID ] = hscNonExistent;
                    return;
            }
        case hscSymSpdPreconditioned:
        case hscSymRmatRaw:
        case hscNonExistent:
        case hscSymSpdRaw:
        case hscSymIndefRaw:
        case hscSymIndefPreconditioned:
        case hscAsymPdRaw:
        case hscAsymPdPreconditioned:
        case hscAsymIndefRaw:
        case hscAsymIndefPreconditioned:
        default:
            std::cout << "[Hybrid Solver]: Internal error.\n"
                      << "[Hybrid Solver]: Failed to free memory from matrix " << matrixID << ".\n";
            status[ matrixID ] = hscNonExistent;
            return;
    }
}

//******************************************************************************
// write a RHS into the global data structure
void
writeRHS( struct global* pointer,
          std::vector<double>& rhs )
{
    unsigned dimension = pointer->dimension;
    double* storage = (double *)calloc(dimension,sizeof(double));
	for(unsigned i=0;i<dimension;i++) storage[i] = rhs[i];
	pointer->rhs = storage;
}

//******************************************************************************
// copy soluton from the global data structure to a vector
// free memory related to this solve in the global data structure
void
getSolutionFreeMemory( struct global* globalpointer,
                       std::vector<double>& solution )
{
    unsigned dimension = globalpointer->dimension;
    solution.resize( dimension );
    for(unsigned i=0;i<dimension;i++) solution[i] = globalpointer->solution[i];
	free(globalpointer->rhs);
	globalpointer->rhs = NULL;
	free(globalpointer->newrhs);
	globalpointer->newrhs = NULL;
	free(globalpointer->solution);
	globalpointer->solution = NULL;
}


//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
// return 0 if sucessful, return 1 if failed
bool
generic( unsigned             matrixID,
         std::vector<double>& rhs,
         double               delta,
         std::vector<double>& solution )
{
    std::cout << "[Hybrid Solver]: Using generic random walks to solve matrix " << matrixID << " with a given RHS.\n";
    double startime= seconds();
    
    if( matrixID >= matrices.size() ){
        std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
        return 1;
    }

    struct global* pointer;

    try{
        switch( status[ matrixID ] ){
            case hscNonExistent:
                std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
                return 1;
            case hscSymRmatRaw:
                pointer = static_cast<struct global*>( matrices[ matrixID ] );
                if( rhs.size() != static_cast<unsigned>( pointer->dimension ) ){
                    std::cout << "[Hybrid Solver]: Invalid RHS. Solving aborted.\n";
                    return 1;
                }
                genericOnSymRaw( pointer, rhs, delta, solution );
                break;
            case hscSymRmatPreconditioned:
                std::cout << "[Hybrid Solver]: This version does not support "
                          << "generic random walks on a preconditioned matrix. "
                          << "Please create another matrix (copy of the original) and call again.\n"
                          << "[Hybrid Solver]: Generic random walks aborted.\n";
                return 1;
            case hscSymSpdPreconditioned:
            case hscSymSpdRaw:
            case hscSymIndefRaw:
            case hscSymIndefPreconditioned:
            case hscAsymPdRaw:
            case hscAsymPdPreconditioned:
            case hscAsymIndefRaw:
            case hscAsymIndefPreconditioned:
            default:
                std::cout << "[Hybrid Solver]: Internal error.\n"
                          << "[Hybrid Solver]: Failed solving.\n";
                return 1;
        }
    }
    catch( Error& error ){
        std::cout << "[Hybrid Solver]: Generic random-walk solve failed due to";
        error.showMessage();
        std::cout << "\n";
        return 1;
    }

    std::cout << "[Hybrid Solver]: Done solving matrix " << matrixID << " with the given RHS by generic random walks.\n"
              << "[Hybrid Solver]: Generic random walks solving runtime " << seconds()-startime << " seconds.\n";
    return 0;
}

//******************************************************************************
// return the "entryIndex"-th entry of the solution vector in "solution"
// accuracy: P[|error| < "delta"] > 99%
// return 0 if sucessful, return 1 if failed
bool
singleEntry( unsigned             matrixID,
             std::vector<double>& rhs,
             double               delta,
             unsigned             entryIndex,
             double&              solution )
{
    std::cout << "[Hybrid Solver]: Solving the " << entryIndex
              <<"th solution entry for matrix " << matrixID << " with a given RHS.\n";
    double startime= seconds();
    
    if( matrixID >= matrices.size() ){
        std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
        return 1;
    }

    struct global* pointer;

    try{
        switch( status[ matrixID ] ){
            case hscNonExistent:
                std::cout << "[Hybrid Solver]: Invalid matrix ID.\n";
                return 1;
            case hscSymRmatRaw:
                pointer = static_cast<struct global*>( matrices[ matrixID ] );
                if( rhs.size() != static_cast<unsigned>( pointer->dimension ) ){
                    std::cout << "[Hybrid Solver]: Invalid RHS. Solving aborted.\n";
                    return 1;
                }
                solution = singleEntry( pointer, rhs, delta, entryIndex );
                break;
            case hscSymRmatPreconditioned:
                std::cout << "[Hybrid Solver]: This version does not support "
                          << "solving single-entry for a preconditioned matrix. "
                          << "Please create another matrix (copy of the original) and call again.\n"
                          << "[Hybrid Solver]: Single-entry solving aborted.\n";
                return 1;
            case hscSymSpdRaw:
            case hscSymSpdPreconditioned:
            case hscSymIndefRaw:
            case hscSymIndefPreconditioned:
            case hscAsymPdRaw:
            case hscAsymPdPreconditioned:
            case hscAsymIndefRaw:
            case hscAsymIndefPreconditioned:
            default:
                std::cout << "[Hybrid Solver]: Internal error.\n"
                          << "[Hybrid Solver]: Failed solving.\n";
                return 1;
        }
    }
    catch( Error& error ){
        std::cout << "[Hybrid Solver]: Solving single-entry failed due to";
        error.showMessage();
        std::cout << "\n";
        return 1;
    }

    std::cout << "[Hybrid Solver]: Done solving the " << entryIndex
              <<"th solution entry for matrix " << matrixID << " with the given RHS.\n"
              << "[Hybrid Solver]: Single-entry solving runtime " << seconds()-startime << " seconds.\n";
    return 0;
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
writePreconditioner( unsigned               matrixID,
                     std::vector<unsigned>& ordering,
                     std::vector<unsigned>& rowIndex, 
                     std::vector<unsigned>& colIndex,
                     std::vector<double>&   value )
{
    if( matrixID >= matrices.size() ){
        std::cout << "[Hybrid Solver]: Invalid matrix ID. "
                  << "getPreconditioner() aborted. \n";
        return 1;
    }

    try{
        switch( status[ matrixID ] ){
            case hscNonExistent:
                std::cout << "[Hybrid Solver]: Invalid matrix ID. "
                          << "getPreconditioner() aborted. \n";
                return 1;
            case hscSymRmatRaw:
                if( precond( matrixID ) ){
                    return 1;
                }
            case hscSymRmatPreconditioned:
                writePreconditioner( static_cast<struct global*>( matrices[ matrixID ] ),
                                     ordering,
                                     rowIndex, 
                                     colIndex,
                                     value );
                break;
            case hscSymSpdRaw:
            case hscSymSpdPreconditioned:
            case hscSymIndefRaw:
            case hscSymIndefPreconditioned:
            case hscAsymPdRaw:
            case hscAsymPdPreconditioned:
            case hscAsymIndefRaw:
            case hscAsymIndefPreconditioned:
            default:
                std::cout << "[Hybrid Solver]: Internal error.\n"
                          << "[Hybrid Solver]: Failed to write out preconditioner.\n";
                return 1;
        }
    }
    catch( Error& error ){
        std::cout << "[Hybrid Solver]: Failed to write out preconditioner due to";
        error.showMessage();
        std::cout << "\n";
        return 1;
    }
    return 0;
}

}

