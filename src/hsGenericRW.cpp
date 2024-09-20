
#include "hsGenericRW.h"

namespace hybridSolver {

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
void
genericOnSymRaw( struct global*       globalpointer, 
                 std::vector<double>& rhs,
                 double               delta,
                 std::vector<double>& solution )
{
    int i,ii,iii,j,k,dimension,realnodenum,*doneflag,*ordering,temp,*nodedegree,probaccumulator;
    int randomseed,kk,randomnumber,currentnodeindex=0,validwalknumber,start,end,middle;
	double temp2,temp4,*currentload,DELTA,sum,var,money;
	float divider = (float)1/(float)IQ, temp2temp;
	struct edge **adjmatrix,*edgepointer;
	struct edge2 **adjmatrix2,*adjhead,*edge2pointer;
	
	DELTA = delta*delta/6.6049;
	if(RANDOMSEED < 0) randomseed=(unsigned)time( NULL );
	else randomseed = RANDOMSEED;
	
	dimension = globalpointer->dimension;
	realnodenum = globalpointer->realnodenum;
	nodedegree = globalpointer->nodedegree;
	adjmatrix = globalpointer->adjmatrix;
    
    currentload=(double *)malloc(dimension*sizeof(double));
	ordering=(int *)malloc(dimension*sizeof(int));
	doneflag = (int *)calloc(realnodenum,sizeof(int)); 
	for(i=0;i<dimension;i++){
	    ordering[i]=i;
	    currentload[i] = rhs[i];
	}
	disorder(ordering, dimension);
	for(i=dimension;i<realnodenum;i++) doneflag[i]=1;

	solution.resize( realnodenum );
	edgepointer = globalpointer->homeregistry;
	while(edgepointer){
		solution[edgepointer->probability] = currentload[edgepointer->nodeindex]/edgepointer->conductance; 
		currentload[edgepointer->nodeindex] = 0;
		edgepointer = edgepointer->next;
	}

    temp=0;
	for(i=0;i<dimension;i++){
	    solution[i] = 0;  
		temp2 = 0;
		edgepointer=adjmatrix[i];
		for(j=0;j<nodedegree[i];j++){
			temp2 += edgepointer->conductance;
			edgepointer=edgepointer->next;
		}
		temp=temp+nodedegree[i];
		edgepointer=adjmatrix[i];
		for(j=0;j<nodedegree[i];j++){
			edgepointer->probability=(int)floor(IM*(edgepointer->conductance)/temp2);
			edgepointer=edgepointer->next;
		}
		currentload[i] = currentload[i] / temp2;
	}

	adjmatrix2 = (struct edge2 **)malloc(dimension*sizeof(struct edge2 *));
	adjhead = (struct edge2 *)malloc(temp*sizeof(struct edge2));
	temp=0;
	for(i=0;i<dimension;i++){
		adjmatrix2[i]=adjhead+temp;
		temp += nodedegree[i];
		probaccumulator=0;
		edgepointer=adjmatrix[i];
		for(j=0;j<nodedegree[i];j++){
			probaccumulator=probaccumulator+edgepointer->probability;
			if(probaccumulator<0) probaccumulator=IM;
			adjmatrix2[i][j].nodeindex = edgepointer->nodeindex;
			adjmatrix2[i][j].probability = probaccumulator;
			adjmatrix2[i][j].prob = (double)(edgepointer->probability)/(double)IM;
			edgepointer=edgepointer->next;
		}
		adjmatrix2[i][nodedegree[i]-1].probability = IM;
	}

    for(ii=0;ii<dimension;ii++){
		i=ordering[ii];
		edge2pointer = (struct edge2 *)malloc(nodedegree[i]*sizeof(struct edge2));
		randomnumber = 0;
		probaccumulator = 0;
		iii = 0;
		for(k=0;k<nodedegree[i];k++){
			currentnodeindex = adjmatrix2[i][k].nodeindex;
			if(!doneflag[currentnodeindex]){
				probaccumulator += adjmatrix2[i][k].probability - randomnumber;
				edge2pointer[iii].nodeindex = currentnodeindex;
				edge2pointer[iii].probability = probaccumulator;
				iii++;
			}
			else{
				solution[i] += adjmatrix2[i][k].prob * solution[ currentnodeindex ];
			}
			randomnumber = adjmatrix2[i][k].probability;
		}
		if(!iii){  
			doneflag[i] = 1;
			free(edge2pointer);
			continue;
		}
		for(k=0;k<iii-1;k++){
			randomnumber = (int)floor((double)edge2pointer[k].probability * (double)IM / (double)probaccumulator);
			if(randomnumber<0) randomnumber=IM;
			edge2pointer[k].probability = randomnumber;
		}
		edge2pointer[iii-1].probability = IM;
		
		validwalknumber=0;
		sum = 0;
		var = 0;
		while(1){
			for(j=0;j<10;j++){
				randomseed ^= MASK;
				temp2temp  = randomseed*divider;
				kk = (int)temp2temp;
				randomseed=IA*(randomseed-kk*IQ)-IR*kk;
				if (randomseed < 0) randomseed += IM;
				randomnumber=randomseed;
				randomseed ^= MASK;
				if(iii < 5)
        			for(k=0;k<iii;k++){
        				if(randomnumber <= edge2pointer[k].probability){
        					currentnodeindex=edge2pointer[k].nodeindex;
        					break;
        				}
        			}
				else{
					start = -1;
					end = iii-1;
					while(end-start > 1){
						middle = (start+end) >> 1;
						if(randomnumber <= edge2pointer[middle].probability) end = middle;
						else start = middle;
					}
					currentnodeindex=edge2pointer[end].nodeindex;
				}
				money = 0;
				for(temp=1;temp<WALKSTEPLIMIT;temp++){
					if(doneflag[currentnodeindex]){
					    money += solution[ currentnodeindex ];
						break;
					}
                    money += currentload[ currentnodeindex ];
					randomseed ^= MASK;
					temp2temp  = randomseed*divider;
					kk = (int)temp2temp;
					randomseed=IA*(randomseed-kk*IQ)-IR*kk;
					if (randomseed < 0) randomseed += IM;
					randomnumber=randomseed;
					randomseed ^= MASK;
    				if(nodedegree[currentnodeindex] < 5)
    					for(k=0;k<nodedegree[currentnodeindex];k++){
    						if(randomnumber <= adjmatrix2[currentnodeindex][k].probability){
    							currentnodeindex=adjmatrix2[currentnodeindex][k].nodeindex;
    							break;
    						}
    					}
    				else{
    					start = -1;
    					end = nodedegree[currentnodeindex]-1;
    					while(end-start > 1){
    						middle = (start+end) >> 1;
    						if(randomnumber <= adjmatrix2[currentnodeindex][middle].probability) end = middle;
    						else start = middle;
    					}
    					currentnodeindex = adjmatrix2[currentnodeindex][end].nodeindex;
    				}
				}
				sum += money;
				var += money*money;
			}
			validwalknumber += 10;
			temp4 = sum*sum;
			temp2 = var - temp4/(double)validwalknumber;
			if( temp2/((double)validwalknumber*(double)validwalknumber) < DELTA ) break;
		}
		temp2 = (double)probaccumulator / ((double)IM * (double)validwalknumber);
        solution[i] += temp2 * sum + currentload[ i ];
		doneflag[i]=1;		
		free(edge2pointer);
	}
    
	free(ordering);
	free(currentload);
	free(adjmatrix2);
	free(adjhead);
	free(doneflag);
	solution.resize( dimension );
}

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
void
genericOnSymPreconditioned( struct global*       globalpointer, 
                            std::vector<double>& rhs,
                            double               delta,
                            std::vector<double>& solution )
{
}

//******************************************************************************
// return the "entryIndex"-th entry of the solution vector
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
double
singleEntry( struct global*       globalpointer,
             std::vector<double>& rhs,
             double               delta,
             unsigned             i )
{
    struct edge2 *edge2pointer;
    struct edge **adjmatrix,*edgepointer;
    int randomseed,kk,randomnumber,probaccumulator,iii,j,k,*nodedegree,dimension;
    int currentnodeindex=0,validwalknumber,start,end,middle,temp;
    double solution,diagonal,sum,var,DELTA,money,*diagonals,tempdiag,tempaccu,probability,temp2,temp4;
    float divider = (float)1/(float)IQ, temp2temp;

	DELTA = delta*delta/6.6049;
	if(RANDOMSEED < 0) randomseed=(unsigned)time( NULL );
	else randomseed = RANDOMSEED;

	dimension = globalpointer->dimension;
	nodedegree = globalpointer->nodedegree;
	adjmatrix = globalpointer->adjmatrix;
	diagonals = globalpointer->diagonal;  

    diagonal = diagonals[i];
    solution = rhs[i] / diagonal;

    edge2pointer = (struct edge2 *)malloc(nodedegree[i]*sizeof(struct edge2));
    probaccumulator = 0;
	iii = 0;
	edgepointer = adjmatrix[i];
	for(k=0;k<nodedegree[i];k++){
	    currentnodeindex = edgepointer->nodeindex;
	    if( currentnodeindex < dimension ){
	        probaccumulator += (int)floor((double)IM*(edgepointer->conductance)/diagonal);
	        if(probaccumulator<0) probaccumulator=IM;
            edge2pointer[iii].nodeindex = currentnodeindex;
			edge2pointer[iii].probability = probaccumulator;
			iii++;
	    }
	    edgepointer = edgepointer->next;
	}
	if(!iii){  
        return solution;
	}
    for(k=0;k<iii-1;k++){
		randomnumber = (int)floor((double)edge2pointer[k].probability * (double)IM / (double)probaccumulator);
		if(randomnumber<0) randomnumber=IM;
		edge2pointer[k].probability = randomnumber;
	}
	edge2pointer[iii-1].probability = IM;
	
	validwalknumber=0;
	sum = 0;
	var = 0;
	while(1){
		for(j=0;j<10;j++){
			randomseed ^= MASK;
			temp2temp  = randomseed*divider;
			kk = (int)temp2temp;
			randomseed=IA*(randomseed-kk*IQ)-IR*kk;
			if (randomseed < 0) randomseed += IM;
			randomnumber=randomseed;
			randomseed ^= MASK;
			if(iii < 5)
    			for(k=0;k<iii;k++){
    				if(randomnumber <= edge2pointer[k].probability){
    					currentnodeindex=edge2pointer[k].nodeindex;
    					break;
    				}
    			}
			else{
				start = -1;
				end = iii-1;
				while(end-start > 1){
					middle = (start+end) >> 1;
					if(randomnumber <= edge2pointer[middle].probability) end = middle;
					else start = middle;
				}
				currentnodeindex=edge2pointer[end].nodeindex;
			}
			money = 0;
			for(temp=1;temp<WALKSTEPLIMIT;temp++){
				if(currentnodeindex >= dimension){ 
					break;
				}
				tempdiag = diagonals[ currentnodeindex ];
                money += rhs[ currentnodeindex ] / tempdiag;

				randomseed ^= MASK;
				temp2temp  = randomseed*divider;
				kk = (int)temp2temp;
				randomseed=IA*(randomseed-kk*IQ)-IR*kk;
				if (randomseed < 0) randomseed += IM;
                probability = ((double)randomseed) / (double)IM; 
				randomseed ^= MASK;

				edgepointer = adjmatrix[ currentnodeindex ];
				tempaccu = 0;
				for(k=0;k<nodedegree[currentnodeindex];k++){
				    tempaccu += edgepointer->conductance / tempdiag;
					if(( probability <= tempaccu )||( k == nodedegree[currentnodeindex]-1 )){
						currentnodeindex = edgepointer->nodeindex;
    					break;
					}
					edgepointer = edgepointer->next;
				}
			}
			sum += money;
			var += money*money;
		}
		validwalknumber += 10;
		temp4 = sum*sum;
		temp2 = var - temp4/(double)validwalknumber;
		if( temp2/((double)validwalknumber*(double)validwalknumber) < DELTA ) break; 
	}
	temp2 = (double)probaccumulator / ((double)IM * (double)validwalknumber);
    solution += temp2 * sum;

	free(edge2pointer);
    return solution;
}

}
