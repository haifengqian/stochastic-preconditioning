
#include "hsSymRmatCore.h"
#include "hsUtility.h"

namespace hybridSolver {

void
shrinkGraph( struct global *globalpointer )
{
    int i,dimension,*shrinkflag,*nodedegree,*ordering,arraysize=0;
    struct edge **adjmatrix;
    struct shrink *shrinkSequence=NULL,*shrinkpointer,*shrinkpointer2;
    struct shrinkEntry *shrinkSequencearray;
	
    dimension = globalpointer->dimension;
    nodedegree = globalpointer->nodedegree;
    adjmatrix = globalpointer->adjmatrix;

    shrinkflag = (int *)calloc(dimension,sizeof(int));
    ordering=(int *)malloc(dimension*sizeof(int));
    for(i=0;i<dimension;i++) ordering[i]=i;
    disorder(ordering, dimension);

    shrinkSequence = shrinkCall1(adjmatrix,nodedegree,shrinkSequence,dimension,shrinkflag);
    shrinkSequence = shrinkCall2(adjmatrix,nodedegree,shrinkSequence,dimension,shrinkflag);
    shrinkSequence = shrinkCall3(ordering,adjmatrix,nodedegree,shrinkSequence,dimension,shrinkflag);

    shrinkpointer = shrinkSequence;
    while(shrinkpointer){
        arraysize++;
        shrinkpointer = shrinkpointer->next;
    }
    shrinkSequencearray = (struct shrinkEntry *)calloc(arraysize,sizeof(struct shrinkEntry));
    shrinkpointer = shrinkSequence;
    i=0;
    while(shrinkpointer){
        shrinkSequencearray[i].type = shrinkpointer->type;
        shrinkSequencearray[i].pointer = shrinkpointer->pointer;
        i++;
        shrinkpointer2 = shrinkpointer->next;
        free(shrinkpointer);
        shrinkpointer = shrinkpointer2;
    }	

    globalpointer->shrinkflag = shrinkflag;
    globalpointer->shrinkSequence = shrinkSequencearray;
    globalpointer->shrinkEntrynumber = arraysize;
    globalpointer->ordering = ordering;
}

struct shrink*
shrinkCall1( struct edge  **adjmatrix,
             int           *nodedegree,
             struct shrink *shrinkSequence,
             int            dimension,
             int           *shrinkflag )
{
    for(int i=0;i<dimension;i++){
        if(shrinkflag[i]) continue;
        if(nodedegree[i]==1){
            shrinkSequence = processSinkPath( i, adjmatrix, nodedegree, shrinkSequence, shrinkflag, dimension );
        }
    }
    return shrinkSequence;
}

struct shrink*
processSinkPath( int            sinkIdx,
                 struct edge  **adjmatrix,
                 int           *nodedegree,
                 struct shrink *shrinkSequence,
                 int           *shrinkflag,
                 int            dimension )
{
    struct sinkPath* pathPointer = (struct sinkPath *)calloc(1, sizeof(struct sinkPath));
    marksinkPath( adjmatrix,nodedegree,shrinkflag,&(pathPointer->inEndIdx),
                  &(pathPointer->exEndIdx),sinkIdx,dimension );
    pathPointer->sinkIdx = sinkIdx;
    shrinksinkPath( adjmatrix,nodedegree,pathPointer,dimension );
    struct shrink* shrinkpointer = shrinkSequence;
    shrinkSequence = (struct shrink *)malloc(1*sizeof(struct shrink));
    shrinkSequence->type = 0;
    shrinkSequence->pointer = (void*)pathPointer;
    shrinkSequence->next = shrinkpointer;    
    return shrinkSequence;
}

void
preprocessing( struct global *globalpointer )
{
    int i,j,dimension,*shrinkflag,*ordering,temp,*nodedegree,probcounter,actualsize=0;
    double temp2,*capacitance;
    struct edge **adjmatrix,*edgepointer;
    struct edge2 **adjmatrix2,*adjhead;
	
    dimension = globalpointer->dimension;
    nodedegree = globalpointer->nodedegree;
    adjmatrix = globalpointer->adjmatrix;
    capacitance = globalpointer->diagonal;

    if(globalpointer->shrinkflag==NULL){   
        shrinkflag = (int *)calloc(dimension,sizeof(int));
        globalpointer->shrinkflag = shrinkflag;
    }
    else shrinkflag = globalpointer->shrinkflag;

    if(globalpointer->ordering==NULL){   
        ordering=(int *)malloc(dimension*sizeof(int));
        for(i=0;i<dimension;i++) ordering[i]=i;
        disorder(ordering, dimension);
        globalpointer->ordering = ordering;
    }
    else ordering = globalpointer->ordering;	

    temp=0;
    for(i=0;i<dimension;i++){
        if(shrinkflag[i]) continue;
        temp2 = 0;
        edgepointer=adjmatrix[i];
        for(j=0;j<nodedegree[i];j++){
            temp2 += edgepointer->conductance;
            edgepointer=edgepointer->next;
        }
        capacitance[i] = temp2; 
        temp=temp+nodedegree[i];
        edgepointer=adjmatrix[i];
        for(j=0;j<nodedegree[i];j++){
            edgepointer->probability=(int)floor(IM*(edgepointer->conductance)/temp2);
            edgepointer=edgepointer->next;
        }
    }

    adjmatrix2 = (struct edge2 **)malloc(dimension*sizeof(struct edge2 *));
    adjhead = (struct edge2 *)malloc(temp*sizeof(struct edge2));

    temp=0;
    for(i=0;i<dimension;i++){
        if(shrinkflag[i]) continue;
        actualsize++;
        adjmatrix2[i]=adjhead+temp;
        temp += nodedegree[i];
        probcounter=0;
        edgepointer=adjmatrix[i];
        for(j=0;j<nodedegree[i];j++){
            probcounter=probcounter+edgepointer->probability;
            if(probcounter<0) probcounter=IM;
            adjmatrix2[i][j].nodeindex = edgepointer->nodeindex;
            adjmatrix2[i][j].probability = probcounter;
            adjmatrix2[i][j].prob = (double)(edgepointer->probability)/(double)IM;
            edgepointer=edgepointer->next;
        }
        adjmatrix2[i][nodedegree[i]-1].probability = IM;
        freeedgepointer(adjmatrix[i]);
        adjmatrix[i] = NULL;
    }
    globalpointer->adjmatrix2 = adjmatrix2;
    globalpointer->actualsize = actualsize;
}

void
randomwalk( struct global *globalpointer,
            double STEPDELTA )
{
    float divider,temp2temp;
    double DELTA,*partiallog,stepsum,stepvar,temp2,temp4,*diagonalscale,*capacitance;
    int randomseed,kk,*appearsequence,*homeappeartimes,*hashnext,*appearsizes,dimension,realnodenum;
    int i,ii,iii,j,k,temp,*ordering,*shrinkflag,probcounter,randomnumber,*hashtable,appearsize;
    int currentnodeindex=0,*doneflag,validwalknumber,*nodedegree;
    int *terminalmnappearsizes,fixterminalcounter;
    int returncounter;
    int hashindex,currentnode;
    struct record **mn,**terminalmn;
    struct edge2 *adjhead,**adjmatrix2;
    int start, end, middle;

    preprocessing( globalpointer );

    dimension = globalpointer->dimension;
    realnodenum = globalpointer->realnodenum;
    ordering = globalpointer->ordering;
    shrinkflag = globalpointer->shrinkflag;
    nodedegree = globalpointer->nodedegree;
    adjmatrix2 = globalpointer->adjmatrix2;
    capacitance = globalpointer->diagonal;
	
    appearsequence = (int *)malloc((realnodenum+1)*sizeof(int));
    homeappeartimes = (int *)malloc((realnodenum+1)*sizeof(int));
    hashnext = (int *)malloc((realnodenum+1)*sizeof(int)); 
    mn = (struct record **)malloc(dimension*sizeof(struct record *));
    appearsizes = (int *)malloc(dimension*sizeof(int));
    terminalmn = (struct record **)malloc(dimension*sizeof(struct record *));
    terminalmnappearsizes = (int *)calloc(dimension,sizeof(int));
    diagonalscale = (double *)calloc(dimension,sizeof(double));
    doneflag = (int *)calloc(realnodenum,sizeof(int)); 
    for(i=dimension;i<realnodenum;i++) doneflag[i]=1;

    DELTA = STEPDELTA*STEPDELTA/6.6049;
    divider = (float)1/(float)IQ;

    if(RANDOMSEED < 0) randomseed=(unsigned)time( NULL );
    else randomseed = RANDOMSEED;

    for(ii=0;ii<dimension;ii++){
        i=ordering[ii];
        if(doneflag[i])	continue;
        if(shrinkflag[i]) continue;
		
        adjhead = (struct edge2 *)malloc(nodedegree[i]*sizeof(struct edge2));
        partiallog = (double *)malloc(nodedegree[i]*sizeof(double));
        randomnumber = 0;
        probcounter = 0;
        iii = 0;
        hashtable=(int *)calloc(1024,sizeof(int));
        appearsize=0;
        for(k=0;k<nodedegree[i];k++){
            currentnodeindex = adjmatrix2[i][k].nodeindex;
            if(!doneflag[currentnodeindex]){
                probcounter += adjmatrix2[i][k].probability - randomnumber;
                adjhead[iii].nodeindex = currentnodeindex;
                adjhead[iii].probability = probcounter;
                iii++;
            }
            else{
                partiallog[k-iii] = adjmatrix2[i][k].prob;

                currentnode = currentnodeindex >> 10;
                hashindex = currentnodeindex - (currentnode << 10);

                currentnode=hashtable[hashindex];
                while(currentnode!=0){
                    if(appearsequence[currentnode]==currentnodeindex){
                        homeappeartimes[currentnode]++;
                        break;
                    }
                    currentnode=hashnext[currentnode];
                }
                if (currentnode == 0){
                    appearsize++;
                    appearsequence[appearsize]=currentnodeindex;
                    homeappeartimes[appearsize]=0;     
                    hashnext[appearsize]=hashtable[hashindex];
                    hashtable[hashindex]=appearsize;
                }

            }
            randomnumber = adjmatrix2[i][k].probability;
        }
        if(!iii){
            fixterminalcounter = 0;
            for(k=0;k<nodedegree[i];k++) if(adjmatrix2[i][k].nodeindex >= dimension) fixterminalcounter++;
            appearsizes[i] = nodedegree[i] - fixterminalcounter;
            if( appearsizes[i] > 0 ) mn[i] = (struct record *)malloc(appearsizes[i]*sizeof(struct record));
            else mn[i] = NULL;
            if(fixterminalcounter>0) terminalmn[i] = (struct record *)malloc(fixterminalcounter*sizeof(struct record));			
            else terminalmn[i] = NULL;
            terminalmnappearsizes[i] = fixterminalcounter;
            fixterminalcounter = 0;
            for(k=0;k<nodedegree[i];k++){
                currentnodeindex = adjmatrix2[i][k].nodeindex;
                if(currentnodeindex < dimension){
                    mn[i][k-fixterminalcounter].nodeindex = currentnodeindex;
                    mn[i][k-fixterminalcounter].value = partiallog[k];
                }
                else{
                    terminalmn[i][fixterminalcounter].nodeindex = currentnodeindex;
                    terminalmn[i][fixterminalcounter].value = partiallog[k];
                    fixterminalcounter++;
                }
            }
            diagonalscale[i] = 1/capacitance[i];
            doneflag[i] = 1;
            free(partiallog);
            free(hashtable);
            free(adjhead);
            continue;
        }
        for(k=0;k<iii-1;k++){
            randomnumber = (int)floor((double)adjhead[k].probability * (double)IM / (double)probcounter);
            if(randomnumber<0) randomnumber=IM;
            adjhead[k].probability = randomnumber;
        }
        adjhead[iii-1].probability = IM;
		
        validwalknumber=0;
        stepsum = 0;
        stepvar = 0;
        returncounter = 0;
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
                        if(randomnumber <= adjhead[k].probability){
                            currentnodeindex=adjhead[k].nodeindex;
                            break;
                        }
				}
                else{
                    start = -1;
                    end = iii-1;
                    while(end-start > 1){
                        middle = (start+end) >> 1;
                        if(randomnumber <= adjhead[middle].probability) end = middle;
                        else start = middle;
                    }
                    currentnodeindex=adjhead[end].nodeindex;
                }
                for(temp=1;temp<WALKSTEPLIMIT;temp++){
                    if(doneflag[currentnodeindex]){

                        currentnode = currentnodeindex >> 10;
                        hashindex = currentnodeindex - (currentnode << 10);

                        currentnode=hashtable[hashindex];
                        while(currentnode!=0){
                            if(appearsequence[currentnode]==currentnodeindex){
                                homeappeartimes[currentnode]++;
                                break;
                            }
                            currentnode=hashnext[currentnode];
                        }
                        if (currentnode == 0){
                            appearsize++;
                            appearsequence[appearsize]=currentnodeindex;
                            homeappeartimes[appearsize]=1;
                            hashnext[appearsize]=hashtable[hashindex];
                            hashtable[hashindex]=appearsize;
                        }

                        break;
                    }
                    if(currentnodeindex == i) returncounter++;
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

                stepsum += temp;
                stepvar += (double)temp*(double)temp;
            }
            validwalknumber += 10;
            temp4 = stepsum*stepsum;
            temp2 = stepvar-temp4/(double)validwalknumber;
            if(temp2 < DELTA * temp4) break;
            if(validwalknumber < 60){
                if(temp2*(double)validwalknumber/(DELTA*temp4) > SWITCHTHRESHOLD){
                    i = UNKNOWN;
                    break;
                }
            }
        }
        if(i==UNKNOWN){
            temp = ii + (int)((double)randomnumber * (double)(dimension-1-ii) * AM);
            if(temp <= ii) temp = ii+1;
            if(temp >= dimension) temp = dimension - 1;
            i=ordering[ii];
            ordering[ii] = ordering[temp];
            ordering[temp] = i;
            ii--;
            free(adjhead);
            free(hashtable);
            free(partiallog);
            continue;
        }
        temp2 = (double)probcounter / ((double)IM * (double)validwalknumber);
        fixterminalcounter = 0;
        for(k=0;k<appearsize;k++) if(appearsequence[k+1] >= dimension) fixterminalcounter++;
        appearsizes[i] = appearsize - fixterminalcounter;
        if(appearsizes[i]>0) mn[i] = (struct record *)malloc(appearsizes[i]*sizeof(struct record));
        else mn[i] = NULL;
        if(fixterminalcounter) terminalmn[i] = (struct record *)malloc(fixterminalcounter*sizeof(struct record));
        else terminalmn[i] = NULL;
        terminalmnappearsizes[i] = fixterminalcounter;
        fixterminalcounter = 0;
        for(k=0;k<appearsize;k++){
            currentnodeindex = appearsequence[k+1];
            if(currentnodeindex < dimension){
                mn[i][k-fixterminalcounter].nodeindex = currentnodeindex;
                mn[i][k-fixterminalcounter].value = (double)homeappeartimes[k+1] * temp2;
            }
            else{
                terminalmn[i][fixterminalcounter].nodeindex = currentnodeindex;
                terminalmn[i][fixterminalcounter].value = (double)homeappeartimes[k+1] * temp2;
                fixterminalcounter++;
            }
        }
        fixterminalcounter = 0;
        for(k=0;k<nodedegree[i] - iii;k++){
            currentnodeindex = appearsequence[k+1];
            if(currentnodeindex < dimension){
                mn[i][k-fixterminalcounter].value += partiallog[k];
            }
            else{
                terminalmn[i][fixterminalcounter].value += partiallog[k];
                fixterminalcounter++;
            }
        }
        diagonalscale[i] = ((double)returncounter * temp2 + 1) / capacitance[i]; 
        doneflag[i]=1;		
        free(adjhead);
        free(hashtable);
        free(partiallog);
    }

    free(appearsequence);
    free(homeappeartimes);
    free(hashnext);
    free(doneflag);

    globalpointer->mn = mn;
    globalpointer->appearsizes = appearsizes;
    globalpointer->terminalmn = terminalmn;
    globalpointer->terminalmnappearsizes = terminalmnappearsizes;
    globalpointer->diagonalscale = diagonalscale;

    postprocessing( globalpointer );
}

void
postprocessing( struct global *globalpointer )
{
    int *inverseordering,temp,i,j,ii,dimension,*ordering,*shrinkflag,actualsize,*newnodedegree,appearsize;
    int *nodedegree,*newappearsizes,*appearsizes,*newterminalmnappearsizes,*terminalmnappearsizes;
    int *nmappearsizes;
    double *newcapacitance,*capacitance,*newdiagonalscale,*diagonalscale;
    struct record **newmn,**mn,**newterminalmn,**terminalmn,**nm;
    struct edge2 **newadjmatrix2,**adjmatrix2;
	
    dimension = globalpointer->dimension;
    ordering = globalpointer->ordering;
    shrinkflag = globalpointer->shrinkflag;
    actualsize = globalpointer->actualsize;
    nodedegree = globalpointer->nodedegree;
    adjmatrix2 = globalpointer->adjmatrix2;
    capacitance = globalpointer->diagonal;
    appearsizes = globalpointer->appearsizes;
    mn = globalpointer->mn;
    terminalmn = globalpointer->terminalmn;
    terminalmnappearsizes = globalpointer->terminalmnappearsizes;
    diagonalscale = globalpointer->diagonalscale;
	
    inverseordering = (int *)calloc(dimension,sizeof(int));
    temp = 0;
    for(ii=0;ii<dimension;ii++){
        i=ordering[ii];
        if(!shrinkflag[i]){
            ordering[temp] = i;
            inverseordering[i] = temp;
            temp++;
        }
    }
    ordering = (int *)realloc(ordering,actualsize*sizeof(int));

    newnodedegree = (int *)calloc(actualsize,sizeof(int));
    newadjmatrix2 = (struct edge2 **)calloc(actualsize,sizeof(struct edge2 *));
    newcapacitance = (double *)calloc(actualsize,sizeof(double));
    newappearsizes = (int *)calloc(actualsize,sizeof(int));
    newmn = (struct record **)calloc(actualsize,sizeof(struct record *));
    newterminalmn = (struct record **)malloc(actualsize*sizeof(struct record *));
    newterminalmnappearsizes = (int *)calloc(actualsize,sizeof(int));
    newdiagonalscale = (double *)calloc(actualsize,sizeof(double));
    for(ii=0;ii<actualsize;ii++){
        i = ordering[ii];
        newnodedegree[ii] = nodedegree[i];
        newadjmatrix2[ii] = adjmatrix2[i];
        for(j=0;j<newnodedegree[ii];j++)
            if(newadjmatrix2[ii][j].nodeindex < dimension)
                newadjmatrix2[ii][j].nodeindex = inverseordering[newadjmatrix2[ii][j].nodeindex];
        newcapacitance[ii] = capacitance[i];
        newappearsizes[ii] = appearsizes[i];
        newmn[ii] = mn[i];
        for(j=0;j<newappearsizes[ii];j++) newmn[ii][j].nodeindex = inverseordering[newmn[ii][j].nodeindex];
        newterminalmn[ii] = terminalmn[i];
        newterminalmnappearsizes[ii] = terminalmnappearsizes[i];
        newdiagonalscale[ii] = diagonalscale[i];
    }

    nm = (struct record **)malloc(actualsize*sizeof(struct record *));
    nmappearsizes = (int *)calloc(actualsize,sizeof(int));
    for(i=0;i<actualsize;i++){
        appearsize = newappearsizes[i];
        for(j=0;j<appearsize;j++) nmappearsizes[newmn[i][j].nodeindex]++;
    }
    for(i=0;i<actualsize;i++)
        if(nmappearsizes[i])
            nm[i] = (struct record *)calloc(nmappearsizes[i],sizeof(struct record));
    free(nmappearsizes);
    nmappearsizes = (int *)calloc(actualsize,sizeof(int));
    for(i=0;i<actualsize;i++){
        appearsize = newappearsizes[i];
        for(j=0;j<appearsize;j++){
            temp = newmn[i][j].nodeindex;
            nm[temp][nmappearsizes[temp]].nodeindex = i;
            nm[temp][nmappearsizes[temp]].value = newmn[i][j].value;
            nmappearsizes[temp]++;
        }
    }

    free(inverseordering);

    free(nodedegree);
    globalpointer->nodedegree = newnodedegree;
    free(adjmatrix2);
    globalpointer->adjmatrix2 = newadjmatrix2;
    free(capacitance);
    globalpointer->diagonal = newcapacitance;
    free(appearsizes);
    globalpointer->appearsizes = newappearsizes;
    free(mn);
    globalpointer->mn = newmn;
    free(terminalmn);
    globalpointer->terminalmn = newterminalmn;
    free(terminalmnappearsizes);
    globalpointer->terminalmnappearsizes = newterminalmnappearsizes;
    free(diagonalscale);
    globalpointer->diagonalscale = newdiagonalscale;

    globalpointer->ordering = ordering; 
    globalpointer->nm = nm;
    globalpointer->nmappearsizes = nmappearsizes;
}

void
processrhs( struct global *globalpointer,
            double TOLERANCE )
{
    int i,dimension,shrinkEntrynumber,actualsize,*ordering;
    double *currentload,*voltage,rhsnorm=0,*capacitance,*newcurrentload;
    struct edge *edgepointer;
    struct shrinkEntry *shrinkSequence;
    struct sinkPath *pathPointer;
    struct nonsinkPath *nonsinkPathPointer;
	
    dimension = globalpointer->dimension;
    actualsize = globalpointer->actualsize;
    currentload = globalpointer->rhs;
    edgepointer = globalpointer->homeregistry;
    shrinkSequence = globalpointer->shrinkSequence;
    shrinkEntrynumber = globalpointer->shrinkEntrynumber;
    capacitance = globalpointer->diagonal;
    ordering = globalpointer->ordering;
	
    for(i=0;i<dimension;i++){
        rhsnorm += currentload[i]*currentload[i];
        currentload[i] = -currentload[i];
    }
    rhsnorm = sqrt(rhsnorm)*TOLERANCE;

    voltage = (double *)calloc(globalpointer->realnodenum,sizeof(double));	
    while(edgepointer){
        voltage[edgepointer->probability] = -currentload[edgepointer->nodeindex]/edgepointer->conductance; 
        currentload[edgepointer->nodeindex] = 0;
        edgepointer = edgepointer->next;
    }
	
    for(i=shrinkEntrynumber-1;i>=0;i--){
        switch(shrinkSequence[i].type){
        case 1:
            nonsinkPathPointer = (struct nonsinkPath *)shrinkSequence[i].pointer;
            rhsshrinknonsinkPath(globalpointer,nonsinkPathPointer);
            break;
        case 0:
            pathPointer = (struct sinkPath *)shrinkSequence[i].pointer;
            rhsshrinksinkPath(globalpointer,pathPointer);
            break;
        case 2:
            rhsStarShrink(globalpointer,(int)(shrinkSequence[i].pointer));
        }
    }

    newcurrentload = (double *)calloc(actualsize,sizeof(double));
    for(i=0;i<actualsize;i++) newcurrentload[i] = currentload[ordering[i]]/capacitance[i];
	
    globalpointer->newrhs = newcurrentload;
    globalpointer->targetnorm = rhsnorm;	
    globalpointer->solution = voltage;
}

void
solve( struct global *globalpointer,
       double TOLERANCE )
{
    double rhsnorm,*voltage,*currentload,temp2,*d,*z,*r,*p,*q,*capacitance;
    double *diagonalscale,*newvoltage,inner_rz,newinner_rz,inner_app,alpha,beta;
    int i,j,*ordering,*nodedegree,dimension,realnodenum;
    int iter,appearsize,*appearsizes,shrinkEntrynumber,*nmappearsizes,*terminalmnappearsizes;
    int actualsize;
    struct edge2 **adjmatrix2;
    struct edge **adjmatrix;
    struct record **mn,**nm,**terminalmn;
    struct shrinkEntry *shrinkSequence;
    struct nonsinkPath *nonsinkPathPointer;
    struct sinkPath *pathPointer;

    processrhs(globalpointer, TOLERANCE);
	
    rhsnorm = globalpointer->targetnorm;
    dimension = globalpointer->dimension;
    ordering = globalpointer->ordering;
    nodedegree = globalpointer->nodedegree;
    adjmatrix2 = globalpointer->adjmatrix2;
    adjmatrix = globalpointer->adjmatrix;
    currentload = globalpointer->newrhs;
    voltage = globalpointer->solution;
    realnodenum = globalpointer->realnodenum;
    capacitance = globalpointer->diagonal;
    mn = globalpointer->mn;
    appearsizes = globalpointer->appearsizes;
    nm = globalpointer->nm;
    nmappearsizes = globalpointer->nmappearsizes;
    terminalmn = globalpointer->terminalmn;
    terminalmnappearsizes = globalpointer->terminalmnappearsizes;
    diagonalscale = globalpointer->diagonalscale;
    shrinkSequence = globalpointer->shrinkSequence;
    shrinkEntrynumber = globalpointer->shrinkEntrynumber;
    actualsize = globalpointer->actualsize;

    d = (double *)calloc(actualsize,sizeof(double));
    z = (double *)calloc(realnodenum,sizeof(double));
    r = (double *)calloc(actualsize,sizeof(double));
    p = (double *)calloc(realnodenum,sizeof(double));
    q = (double *)calloc(actualsize,sizeof(double));
	
    for(i=actualsize-1;i>=0;i--){
        temp2 = currentload[i]*capacitance[i];
        appearsize = nmappearsizes[i];
        for(j=0;j<appearsize;j++) temp2 += d[nm[i][j].nodeindex] * nm[i][j].value;
        d[i] = temp2;
    }
    for(i=0;i<actualsize;i++){
        temp2 = d[i]*diagonalscale[i];
        appearsize = appearsizes[i];
        for(j=0;j<appearsize;j++) temp2 += d[mn[i][j].nodeindex] * mn[i][j].value;
        appearsize = terminalmnappearsizes[i];
        for(j=0;j<appearsize;j++) temp2 -= voltage[terminalmn[i][j].nodeindex] * terminalmn[i][j].value;
        d[i] = temp2;
        voltage[i] = -temp2;
    }
	
    for(i=0;i<actualsize;i++){
        temp2 = voltage[i] + currentload[i];
        for(j=0;j<nodedegree[i];j++) temp2 -= adjmatrix2[i][j].prob * voltage[adjmatrix2[i][j].nodeindex];
        voltage[i] -= temp2; 
    }
	
    for(i=0;i<actualsize;i++){
        temp2 = voltage[i] + currentload[i];
        for(j=0;j<nodedegree[i];j++) temp2 -= adjmatrix2[i][j].prob * voltage[adjmatrix2[i][j].nodeindex];
        r[i] = -temp2*capacitance[i];
    }
	
    for(i=actualsize-1;i>=0;i--){
        temp2 = -r[i];
        appearsize = nmappearsizes[i];
        for(j=0;j<appearsize;j++) temp2 += d[nm[i][j].nodeindex] * nm[i][j].value;
        d[i] = temp2;
    }
    inner_rz = 0;
    for(i=0;i<actualsize;i++){
        temp2 = d[i]*diagonalscale[i];
        appearsize = appearsizes[i];
        for(j=0;j<appearsize;j++) temp2 += d[mn[i][j].nodeindex] * mn[i][j].value;
        d[i] = temp2;
        z[i] = -temp2;
        p[i] = -temp2;
        inner_rz -= r[i]*temp2;
    }
	
    for(iter=2;1;iter++){
        inner_app = 0;
        for(i=0;i<actualsize;i++){
            temp2 = p[i];
            for(j=0;j<nodedegree[i];j++) temp2 -= adjmatrix2[i][j].prob * p[adjmatrix2[i][j].nodeindex];
            q[i] = temp2*capacitance[i];
            inner_app += p[i]*q[i];
        }
        alpha = inner_rz / inner_app;
        for(i=0;i<actualsize;i++){
            voltage[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        temp2 = sqrt(inner_rz);

        if(temp2 < rhsnorm) break;

        for(i=actualsize-1;i>=0;i--){
            temp2 = -r[i];
            appearsize = nmappearsizes[i];
            for(j=0;j<appearsize;j++) temp2 += d[nm[i][j].nodeindex] * nm[i][j].value;
            d[i] = temp2;
        }
        newinner_rz = 0;
        for(i=0;i<actualsize;i++){
            temp2 = d[i]*diagonalscale[i];
            appearsize = appearsizes[i];
            for(j=0;j<appearsize;j++) temp2 += d[mn[i][j].nodeindex] * mn[i][j].value;
            d[i] = temp2;
            z[i] = -temp2;
            newinner_rz -= r[i]*temp2;
        }
        beta = newinner_rz / inner_rz;
        inner_rz = newinner_rz;
        for(i=0;i<actualsize;i++) p[i] = p[i]*beta + z[i];
    }

    newvoltage = (double *)calloc(realnodenum,sizeof(double));
    for(i=0;i<actualsize;i++) newvoltage[ordering[i]] = voltage[i];
    for(i=dimension;i<realnodenum;i++) newvoltage[i] = voltage[i];

    currentload = globalpointer->rhs;
    for(i=0;i<shrinkEntrynumber;i++){
        switch(shrinkSequence[i].type){
        case 1:
            nonsinkPathPointer = (struct nonsinkPath *)shrinkSequence[i].pointer;
            recovernonsinkPath(adjmatrix,currentload,newvoltage,nonsinkPathPointer);
            break;
        case 0:
            pathPointer = (struct sinkPath *)shrinkSequence[i].pointer;
            recoversinkPath(adjmatrix,currentload,newvoltage,pathPointer);
            break;
        case 2:
            starRecover((int)(shrinkSequence[i].pointer),adjmatrix,currentload,newvoltage);
        }
    }

    free(d);
    free(z);
    free(r);
    free(p);
    free(q);

    free(voltage);
    globalpointer->solution = newvoltage;
}

void
freeedgepointer( struct edge  *edgepointer )
{
    struct edge  *edgepointer2;
    while(edgepointer!=NULL){
        edgepointer2=edgepointer->next;
        free(edgepointer);
        edgepointer=edgepointer2;
    }
}

void
disorder( int *ports2,
          int  portnumber2 )
{
    int i, num1,num2,exch;
    int randomseed, kk, randomnumber;
    float divider,temp2temp;

    randomseed=5;
    divider = (float)1/(float)IQ;
    for (i=0;i<portnumber2;i++){

        randomseed ^= MASK;
        temp2temp  = randomseed*divider;
        kk = (int)temp2temp;
        randomseed=IA*(randomseed-kk*IQ)-IR*kk;
        if (randomseed < 0) randomseed += IM;
        randomnumber=randomseed;
        randomseed ^= MASK;

        num1=(int)floor((double)randomnumber*(double)portnumber2/(double)IM);
        if(num1<0)num1=0;
        if(num1>=portnumber2)num1=portnumber2-1;

        randomseed ^= MASK;
        temp2temp  = randomseed*divider;
        kk = (int)temp2temp;
        randomseed=IA*(randomseed-kk*IQ)-IR*kk;
        if (randomseed < 0) randomseed += IM;
        randomnumber=randomseed;
        randomseed ^= MASK;

        num2=(int)floor((double)randomnumber*(double)portnumber2/(double)IM);
        if(num2<0)num2=0;
        if(num2>=portnumber2)num2=portnumber2-1;

        exch=ports2[num1];
        ports2[num1]=ports2[num2];
        ports2[num2]=exch;
    }
}

void marksinkPath( struct edge **adjmatrix,
                   int *nodedegree,
                   int *shrinkflag,
                   int *portaddrs,
                   int *portexaddrs,
                   int currentnode,
                   int dimension )
{
    while(1){
        shrinkflag[currentnode] = 1;
        struct edge *edgepointer = adjmatrix[currentnode];
        int nextnode = edgepointer->nodeindex;     
        if((nextnode >= dimension)||(nodedegree[nextnode] >= 3)) {  
            *portaddrs = currentnode;
            *portexaddrs = nextnode;
            return;
        }
        else if( ! shrinkflag[ nextnode ] ){   
            currentnode = nextnode;
            continue;
        }
        if( nodedegree[ currentnode ] != 2 ){  
            throw Error( hscCannotFindPathEnd );
        }
        nextnode = edgepointer->next->nodeindex; 
        if((nextnode >= dimension)||(nodedegree[nextnode] >= 3)) { 
            *portaddrs = currentnode;
            *portexaddrs = nextnode;
            return;
        }
        else{
            if( shrinkflag[ nextnode ] ){     
                throw Error( hscCannotFindPathEnd );
            }
            currentnode = nextnode;
            continue;
        }
    }        
}

void
marknonsinkPath( struct edge       **adjmatrix,
                 int                *nodedegree,
                 int                *shrinkflag,
                 struct nonsinkPath *pathPointer,
                 int                 i,
                 int                 dimension )
{
    shrinkflag[i] = 1;
    struct edge *edgepointer = adjmatrix[i];
    int nodeindex = edgepointer->nodeindex;
    
    if((nodeindex >= dimension)||(nodedegree[nodeindex]>=3)){
        pathPointer->exEndIdx = nodeindex;
        pathPointer->inEndIdx = i;
    }
    else if(!shrinkflag[nodeindex]){
        marksinkPath( adjmatrix, nodedegree, shrinkflag, &(pathPointer->inEndIdx),
                      &(pathPointer->exEndIdx),nodeindex,dimension);
    }
    else throw Error( hscCannotMarkNonsinkPath );
    
    edgepointer = edgepointer->next;
    
    nodeindex = edgepointer->nodeindex;
    if((nodeindex >= dimension)||(nodedegree[nodeindex]>=3)){
        pathPointer->exEndIdx2 = nodeindex;
        pathPointer->inEndIdx2 = i;
    }
    else if(!shrinkflag[nodeindex]){
        marksinkPath( adjmatrix, nodedegree, shrinkflag, &(pathPointer->inEndIdx2),
                      &(pathPointer->exEndIdx2),nodeindex,dimension);
    }
    else throw Error( hscCannotMarkNonsinkPath );
}

void
shrinksinkPath( struct edge    **adjmatrix,
                int             *nodedegree,
                struct sinkPath *pathPointer,
                int              dimension )
{
    int node = pathPointer->exEndIdx;
    if( node < dimension ) cutonedirectionedge( node,
                                                pathPointer->inEndIdx,
                                                adjmatrix,nodedegree );
}

void
rhsshrinksinkPath( struct global   *globalpointer,
                   struct sinkPath *pathPointer )
{
    struct edge** adjmatrix = globalpointer->adjmatrix;
    double* currentload = globalpointer->rhs;
    int dimension = globalpointer->dimension;
    int port = pathPointer->exEndIdx;
    int currentnode  = pathPointer->sinkIdx;
    int oldnode = UNKNOWN;
    int nextnode = dimension;
    double load = 0;
    while(1){
        load += currentload[currentnode];
        struct edge *edgepointer = adjmatrix[currentnode];
        if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
        nextnode = edgepointer->nodeindex;
        if(nextnode == port) break;
        oldnode = currentnode;
        currentnode = nextnode;
    }
    if(nextnode < dimension) currentload[nextnode] += load;
    pathPointer->load = load;
}

void
shrinknonsinkPath( struct edge       **adjmatrix,
                   int                *nodedegree,
                   struct nonsinkPath *pathPointer,
                   int                 dimension )
{
    int inEnd1 = pathPointer->inEndIdx;
    int exEnd1 = pathPointer->exEndIdx;
    int inEnd2 = pathPointer->inEndIdx2;
    int exEnd2 = pathPointer->exEndIdx2;
    int currentnode = inEnd1;
    int oldnode = exEnd1;
    if( exEnd1 < dimension ) cutonedirectionedge(exEnd1,inEnd1,adjmatrix,nodedegree);
    if( exEnd2 < dimension ) cutonedirectionedge(exEnd2,inEnd2,adjmatrix,nodedegree);
    double counter;
    if( adjmatrix[inEnd1]->nodeindex == exEnd1 ) counter = 1/(adjmatrix[inEnd1]->conductance);
    else counter = 1/(adjmatrix[inEnd1]->next->conductance);
    while(1){
        struct edge *edgepointer = adjmatrix[currentnode];
        if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
        int nextnode = edgepointer->nodeindex;
        counter += 1/(edgepointer->conductance);
        if(nextnode == exEnd2) break;
        oldnode = currentnode;
        currentnode = nextnode;
    }
    pathPointer->res = counter;
    if(exEnd1 != exEnd2) addedge(exEnd1,exEnd2,1/counter,adjmatrix,nodedegree,dimension);
}

void
rhsshrinknonsinkPath(struct global      *globalpointer,
                     struct nonsinkPath *nonsinkPathPointer)
{
    int inEnd1 = nonsinkPathPointer->inEndIdx;
    int exEnd1 = nonsinkPathPointer->exEndIdx;
    int exEnd2 = nonsinkPathPointer->exEndIdx2;
    int currentnode = inEnd1;
    int oldnode = exEnd1;
    double counter = nonsinkPathPointer->res;
    struct edge** adjmatrix = globalpointer->adjmatrix;
    double* currentload = globalpointer->rhs;
    int dimension = globalpointer->dimension;
    double counter2,accu1=0,accu2=0;
    if( adjmatrix[inEnd1]->nodeindex == exEnd1 ) counter2 = 1/(adjmatrix[inEnd1]->conductance);
    else counter2 = 1/(adjmatrix[inEnd1]->next->conductance);
    while(1){
        accu1 += currentload[currentnode] * (counter - counter2);
        accu2 += currentload[currentnode] * counter2;
        struct edge* edgepointer = adjmatrix[currentnode];
        if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
        int nextnode = edgepointer->nodeindex;
        if(nextnode == exEnd2) break;
        counter2 += 1/(edgepointer->conductance);
        oldnode = currentnode;
        currentnode = nextnode;
    }
    if(exEnd1 < dimension) currentload[exEnd1] += accu1/counter;
    if(exEnd2 < dimension) currentload[exEnd2] += accu2/counter;
    nonsinkPathPointer->load1 = accu1/counter;
    nonsinkPathPointer->load2 = accu2/counter;
}

void
recoversinkPath( struct edge    **adjmatrix,
                 double          *currentload,
                 double          *voltage,
                 struct sinkPath *pathPointer )
{
    int oldnode = pathPointer->exEndIdx;
    int currentnode = pathPointer->inEndIdx;
    int sink = pathPointer->sinkIdx;
    double load = pathPointer->load;
    double nodevoltage = voltage[ oldnode ];
    if(load == 0 ){
        while(1){
            voltage[ currentnode ] = nodevoltage;
            if( currentnode==sink ) return;
            struct edge* edgepointer = adjmatrix[currentnode];
            if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
            int nextnode = edgepointer->nodeindex;
            oldnode = currentnode;
            currentnode = nextnode;
        }
        throw Error( hscCannotRecoversinkPath );
    }
    while(1){
        struct edge* edgepointer = adjmatrix[ currentnode ];
        if( currentnode == sink ){
            nodevoltage -= load / ( edgepointer->conductance );
            voltage[ currentnode ] = nodevoltage;
            return;
        }
        int nextnode=0;
        for(unsigned i=0; i<2; i++){
            if( edgepointer->nodeindex == oldnode ){
                nodevoltage -= load / ( edgepointer->conductance );
                voltage[ currentnode ] = nodevoltage;
            }
            else nextnode = edgepointer->nodeindex;
            edgepointer = edgepointer->next;
        }
        load -= currentload[ currentnode ];
        oldnode = currentnode;
        currentnode = nextnode;
    }
    throw Error( hscCannotRecoversinkPath );
}

void
recovernonsinkPath( struct edge       **adjmatrix,
                    double             *currentload,
                    double             *voltage,
                    struct nonsinkPath *pathPointer )
{
    int inEnd1 = pathPointer->inEndIdx;
    int exEnd1 = pathPointer->exEndIdx;
    int inEnd2 = pathPointer->inEndIdx2;
    int exEnd2 = pathPointer->exEndIdx2;
    double nodevoltage = voltage[exEnd1];
    double load = ( voltage[exEnd1] - voltage[exEnd2] ) / pathPointer->res + pathPointer->load1;
    int currentnode = inEnd1;
    int oldnode = exEnd1;
    while(1){
        struct edge* edgepointer = adjmatrix[currentnode];
        int nextnode=0;
        for(unsigned i=0; i<2; i++){
            if( edgepointer->nodeindex == oldnode){
                nodevoltage -= load / ( edgepointer->conductance );
                voltage[ currentnode ] = nodevoltage;
            }
            else nextnode = edgepointer->nodeindex;
            edgepointer = edgepointer->next;
        }
        if( currentnode == inEnd2 ) return;
        load -= currentload[ currentnode ];
        oldnode = currentnode;
        currentnode = nextnode;
    }
    throw Error( hscCannotRecoverPath );
}

void
addedge( int           fromnodeindex,
         int           tonodeindex,
         double        conductance,
         struct edge **adjmatrix,
         int          *nodedegree,
         int           dimension )
{
    struct edge *edgepointer;
    int i;
    if(fromnodeindex < dimension){
        edgepointer=adjmatrix[fromnodeindex];
        for(i=0;i<nodedegree[fromnodeindex];i++){
            if(edgepointer->nodeindex == tonodeindex){
                edgepointer->conductance += conductance;
                break;
            }
            edgepointer = edgepointer->next;
        }
        if (i==nodedegree[fromnodeindex]){
            edgepointer=adjmatrix[fromnodeindex];
            adjmatrix[fromnodeindex]=(struct edge *)malloc(1*sizeof(struct edge));
            adjmatrix[fromnodeindex]->nodeindex=tonodeindex;
            adjmatrix[fromnodeindex]->conductance=conductance;
            adjmatrix[fromnodeindex]->next=edgepointer;
            nodedegree[fromnodeindex]++;
        }
    }
    if(tonodeindex < dimension){
        edgepointer=adjmatrix[tonodeindex];
        for(i=0;i<nodedegree[tonodeindex];i++){
            if(edgepointer->nodeindex == fromnodeindex){
                edgepointer->conductance += conductance;
                break;
            }
            edgepointer = edgepointer->next;
        }
        if (i==nodedegree[tonodeindex]){
            edgepointer=adjmatrix[tonodeindex];
            adjmatrix[tonodeindex]=(struct edge *)malloc(1*sizeof(struct edge));
            adjmatrix[tonodeindex]->nodeindex=fromnodeindex;
            adjmatrix[tonodeindex]->conductance=conductance;
            adjmatrix[tonodeindex]->next=edgepointer;
            nodedegree[tonodeindex]++;
        }
    }
}

struct shrink*
shrinkCall2( struct edge  **adjmatrix,
             int           *nodedegree,
             struct shrink *shrinkSequence,
             int            dimension,
             int           *shrinkflag )
{
    int* stack = (int *)malloc(2*dimension*sizeof(int)); 
    int stackpointer=0;
    for(int ii=0; ii<dimension; ii++){
        if( shrinkflag[ii] ) continue;
        if( nodedegree[ii] != 2 ) continue;
        stack[ stackpointer ] = ii;
        stackpointer++;
        int node1,node2;
        while( stackpointer ){
            stackpointer--;
            int i = stack[ stackpointer ];
            if( shrinkflag[i] ) continue;
            switch( nodedegree[i] ){
            case 1:
                shrinkSequence = processSinkPath( i, adjmatrix, nodedegree, shrinkSequence, shrinkflag, dimension );
                node1 = ((struct sinkPath *)(shrinkSequence->pointer))->exEndIdx;
                if( ( node1 < dimension )&&( nodedegree[node1] <= 2 ) ){
                    stack[stackpointer] = node1;
                    stackpointer++;
                }
                break;
            case 2:
                shrinkSequence = processNonSinkPath( i, adjmatrix, nodedegree, shrinkSequence, shrinkflag, dimension );
                node1 = ((struct nonsinkPath *)(shrinkSequence->pointer))->exEndIdx;
                node2 = ((struct nonsinkPath *)(shrinkSequence->pointer))->exEndIdx2;
                if((node2 < dimension)&&(nodedegree[node2]<=2)){
                    stack[stackpointer] = node2;
                    stackpointer++;
                }
                if((node1 < dimension)&&(nodedegree[node1]<=2)){
                    stack[stackpointer] = node1;
                    stackpointer++;
                }
            }
        }
    }
    free(stack);
    return shrinkSequence;
}

struct shrink*
processNonSinkPath( int            nodeindex,
                    struct edge  **adjmatrix,
                    int           *nodedegree,
                    struct shrink *shrinkSequence,
                    int           *shrinkflag,
                    int            dimension )
{
    struct nonsinkPath* nonsinkPathPointer = (struct nonsinkPath*)calloc(1, sizeof(struct nonsinkPath));
    marknonsinkPath(adjmatrix,nodedegree,shrinkflag,nonsinkPathPointer,nodeindex,dimension);
    shrinknonsinkPath(adjmatrix,nodedegree,nonsinkPathPointer,dimension);
    struct shrink* shrinkpointer = shrinkSequence;
    shrinkSequence = (struct shrink*)malloc(1*sizeof(struct shrink));
    shrinkSequence->type = 1;
    shrinkSequence->pointer = (void*)nonsinkPathPointer;
    shrinkSequence->next = shrinkpointer;
    return shrinkSequence;
}

struct shrink*
shrinkCall3( int           *ordering,
             struct edge  **adjmatrix,
             int           *nodedegree,
             struct shrink *shrinkSequence,
             int            dimension,
             int           *shrinkflag )
{
    int* stack = (int *)malloc(2*dimension*sizeof(int)); 
    int stackpointer=0;
    for(int ii=0; ii<dimension; ii++){
        int i = ordering[ ii ];
        if( shrinkflag[i] ) continue;
        if( nodedegree[i] != 3 ) continue;
        stack[ stackpointer ] = i;
        stackpointer++;
        while( stackpointer ){
            stackpointer--;
            i = stack[ stackpointer ];
            if( shrinkflag[i] ) continue;
            int node1,node2,node3;
            switch( nodedegree[i] ){
            case 1:
                shrinkSequence = processSinkPath( i, adjmatrix, nodedegree, shrinkSequence, shrinkflag, dimension );
                node1 = ((struct sinkPath*)(shrinkSequence->pointer))->exEndIdx;
                if( ( node1 < dimension )&&( nodedegree[node1] <= 3 ) ){
                    stack[ stackpointer ] = node1;
                    stackpointer++;
                }
                break;
            case 2:
                shrinkSequence = processNonSinkPath( i, adjmatrix, nodedegree, shrinkSequence, shrinkflag, dimension );
                node1 = ((struct nonsinkPath*)(shrinkSequence->pointer))->exEndIdx;
                node2 = ((struct nonsinkPath*)(shrinkSequence->pointer))->exEndIdx2;
                if( ( node1 < dimension )&&( nodedegree[node1] == 3 ) ){
                    stack[ stackpointer ] = node1;
                    stackpointer++;
                }
                if( ( node2 < dimension )&&( nodedegree[node2] == 3 )&&( node2 != node1 ) ){
                    stack[ stackpointer ] = node2;
                    stackpointer++;
                }
                if( ( node1 < dimension )&&( nodedegree[node1] == 2 ) ){
                    stack[ stackpointer ] = node1;
                    stackpointer++;
                }
                if( ( node2 < dimension )&&( nodedegree[node2] == 2 )&&( node2 != node1 ) ){
                    stack[ stackpointer ] = node2;
                    stackpointer++;
                }
                if( ( node1 < dimension )&&( nodedegree[node1] == 1 ) ){ 
                    stack[ stackpointer ] = node1;
                    stackpointer++;
                }
                break;
            case 3:
                shrinkflag[i] = 1;
                struct shrink* shrinkpointer = shrinkSequence;
                shrinkSequence = (struct shrink*)malloc(1*sizeof(struct shrink));
                shrinkSequence->type = 2;
                shrinkSequence->pointer = (void*)i;
                shrinkSequence->next = shrinkpointer;
                struct edge* edgepointer = adjmatrix[i];
                node1 = edgepointer->nodeindex;
                double cond1 = edgepointer->conductance;
                edgepointer = edgepointer->next;
                node2 = edgepointer->nodeindex;
                double cond2 = edgepointer->conductance;
                edgepointer = edgepointer->next;
                node3 = edgepointer->nodeindex;
                double cond3 = edgepointer->conductance;
                double sum = cond1 + cond2 + cond3;
                double newcond = cond1 * cond2 / sum;
                addedge(node1,node2,newcond,adjmatrix,nodedegree,dimension);
                newcond = cond1 * cond3 / sum;
                addedge(node1,node3,newcond,adjmatrix,nodedegree,dimension);
                newcond = cond3 * cond2 / sum;
                addedge(node3,node2,newcond,adjmatrix,nodedegree,dimension);
                if( node1 < dimension ) cutonedirectionedge(node1,i,adjmatrix,nodedegree);
                if( node2 < dimension ) cutonedirectionedge(node2,i,adjmatrix,nodedegree);
                if( node3 < dimension ) cutonedirectionedge(node3,i,adjmatrix,nodedegree);
                if( ( node1 < dimension )&&( nodedegree[node1] == 3 ) ){
                    stack[ stackpointer ] = node1;
                    stackpointer++;
                }
                if( ( node2 < dimension )&&( nodedegree[node2] == 3 ) ){
                    stack[ stackpointer ] = node2;
                    stackpointer++;
                }
                if( ( node3 < dimension )&&( nodedegree[node3] == 3 ) ){
                    stack[ stackpointer ] = node3;
                    stackpointer++;
                }
                if( ( node1 < dimension )&&( nodedegree[node1] == 2 ) ){
                    stack[ stackpointer ] = node1;
                    stackpointer++;
                }
                if( ( node2 < dimension )&&( nodedegree[node2] == 2 ) ){
                    stack[ stackpointer ] = node2;
                    stackpointer++;
                }
                if( ( node3 < dimension )&&( nodedegree[node3] == 2 ) ){
                    stack[ stackpointer ] = node3;
                    stackpointer++;
                }
            }
        }
    }
    free(stack);
    return shrinkSequence;
}

void
rhsStarShrink( struct global *globalpointer,
               int            nodeindex )
{
	struct edge* edgepointer = globalpointer->adjmatrix[nodeindex];
	double* currentload = globalpointer->rhs;
	int dimension = globalpointer->dimension;
	int node1 = edgepointer->nodeindex;
	double cond1 = edgepointer->conductance;
	edgepointer = edgepointer->next;
	int node2 = edgepointer->nodeindex;
	double cond2 = edgepointer->conductance;
	edgepointer = edgepointer->next;
	int node3 = edgepointer->nodeindex;
	double cond3 = edgepointer->conductance;
	double sum = cond1 + cond2 + cond3;
	double load = currentload[nodeindex];
	if(node1 < dimension) currentload[node1] += load * cond1 / sum;
	if(node2 < dimension) currentload[node2] += load * cond2 / sum;
	if(node3 < dimension) currentload[node3] += load * cond3 / sum;
}

void
starRecover( int           nodeindex,
             struct edge **adjmatrix,
             double       *currentload,
             double       *voltage )
{
    struct edge *edgepointer = adjmatrix[nodeindex];
    int node1 = edgepointer->nodeindex;
    double cond1 = edgepointer->conductance;
    edgepointer = edgepointer->next;
    int node2 = edgepointer->nodeindex;
    double cond2 = edgepointer->conductance;
    edgepointer = edgepointer->next;
    int node3 = edgepointer->nodeindex;
    double cond3 = edgepointer->conductance;
    double sum = cond1 + cond2 + cond3;
    voltage[ nodeindex ] = ( voltage[node1] * cond1 + voltage[node2] * cond2
                             + voltage[node3]* cond3 - currentload[nodeindex] ) / sum;
}

void
cutonedirectionedge( int           fromnodeindex,
                     int           tonodeindex,
                     struct edge **adjmatrix,
                     int          *nodedegree)
{
    struct edge * edgepointer, *edgepointer2;
    int i;
    edgepointer = adjmatrix[fromnodeindex];
    if(edgepointer->nodeindex == tonodeindex) {
        adjmatrix[fromnodeindex] = edgepointer->next;
        nodedegree[fromnodeindex]--;
        free(edgepointer);
        return;
    }
    else for(i=0;i<nodedegree[fromnodeindex];i++){
            edgepointer2 = edgepointer->next;
            if(edgepointer2->nodeindex == tonodeindex){
                edgepointer->next = edgepointer2->next;
                nodedegree[fromnodeindex]--;
                free(edgepointer2);
                return;
            }
            edgepointer = edgepointer->next;
        }
    throw Error( hscCannotCutEdge );
}

double
seconds()
{
    return (double)(clock())/CLOCKS_PER_SEC;
}

void
freespace( struct global* globalpointer )
{
    int i,dimension,minvalue,minindex,actualsize;
	
    dimension = globalpointer->dimension;
    actualsize = globalpointer->actualsize;

    for(i=0;i<actualsize;i++){
        if(globalpointer->terminalmnappearsizes[i]) free(globalpointer->terminalmn[i]);
    }
    free(globalpointer->terminalmn);

    for(i=0;i<actualsize;i++){
        if(globalpointer->nmappearsizes[i]) free(globalpointer->nm[i]);
    }
    free(globalpointer->nm);

    for(i=0;i<actualsize;i++){
        if(globalpointer->appearsizes[i]) free(globalpointer->mn[i]);
    }
    free(globalpointer->mn);
	
    minvalue = dimension + 1;
    minindex = UNKNOWN;
    for(i=0;i<actualsize;i++){
        if(minvalue > globalpointer->ordering[i]){
            minvalue = globalpointer->ordering[i];
            minindex = i;
        }
    }
    if(minindex != UNKNOWN) free(globalpointer->adjmatrix2[minindex]);
    free(globalpointer->adjmatrix2);

    for(i=0;i<dimension;i++) freeedgepointer(globalpointer->adjmatrix[i]);
    free(globalpointer->adjmatrix);

    if(globalpointer->solution != NULL) free(globalpointer->solution);
    free(globalpointer->diagonalscale);
    free(globalpointer->terminalmnappearsizes);
    free(globalpointer->nmappearsizes);
    free(globalpointer->appearsizes);
    free(globalpointer->ordering);
    for( i=0; i<globalpointer->shrinkEntrynumber; i++ ){
        if( globalpointer->shrinkSequence[i].type < 2){
            free( globalpointer->shrinkSequence[i].pointer );
        }
    }
    free(globalpointer->shrinkSequence);
    free(globalpointer->shrinkflag);
    freeedgepointer(globalpointer->homeregistry);
    free(globalpointer->diagonal);
    if(globalpointer->newrhs != NULL) free(globalpointer->newrhs);
    if(globalpointer->rhs != NULL) free(globalpointer->rhs);
    free(globalpointer->nodedegree);
    free(globalpointer);
}


static void
getOrdering( struct nonsinkPath    *pathPointer,
             struct edge          **adjmatrix,
             std::vector<unsigned>& storage )
{
    storage.clear();
    int inEnd1 = pathPointer->inEndIdx;
    int exEnd1 = pathPointer->exEndIdx;
    int inEnd2 = pathPointer->inEndIdx2;
    int currentnode = inEnd1;
    int oldnode = exEnd1;
    while(1){
        storage.push_back( currentnode );
        struct edge* edgepointer = adjmatrix[currentnode];
        if( currentnode == inEnd2 ) return;
        if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
        oldnode = currentnode;
        currentnode = edgepointer->nodeindex;
    }
}


static void
getOrdering( struct sinkPath       *pathPointer,
             struct edge          **adjmatrix,
             std::vector<unsigned>& storage )
{
    std::vector<int> buffer;
    buffer.clear();
    int oldnode = pathPointer->exEndIdx;
    int currentnode = pathPointer->inEndIdx;
    int sink = pathPointer->sinkIdx;
    while(1){
        buffer.push_back( currentnode );
        if( currentnode == sink ) break;
        struct edge* edgepointer = adjmatrix[currentnode];
        if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
        oldnode = currentnode;
        currentnode = edgepointer->nodeindex;
    }

    unsigned size = buffer.size();
    storage.resize( size );
    for(unsigned i=0;i<size;i++){
        assert( buffer[size-i-1] >= 0 );
        storage[i] = buffer[size-i-1];
    }
}


static void
getOrdering( struct global*         globalpointer,
             std::vector<unsigned>& ordering )
{
    ordering.clear();

    struct shrinkEntry* shrinkSequence = globalpointer->shrinkSequence;
    unsigned shrinkEntrynumber = globalpointer->shrinkEntrynumber;
    std::vector<unsigned> storage;
    unsigned ii;
    for(int i=shrinkEntrynumber-1; i>=0; i--){ 
        switch(shrinkSequence[i].type){
            case 1:
                getOrdering( (struct nonsinkPath *)shrinkSequence[i].pointer,
                             globalpointer->adjmatrix,
                             storage );
                for(ii=0;ii<storage.size();ii++){
                    ordering.push_back( storage[ii] );
                }
                break;
            case 0:
                getOrdering( (struct sinkPath *)shrinkSequence[i].pointer,
                             globalpointer->adjmatrix,
                             storage );
                for(ii=0;ii<storage.size();ii++){
                    ordering.push_back( storage[ii] );
                }
                break;
            case 2:
                ordering.push_back( (unsigned)(shrinkSequence[i].pointer) );
        }
	}

    int actualsize = globalpointer->actualsize;
    int* mnOrdering = globalpointer->ordering;
    for(int iii = actualsize-1; iii >= 0; iii--){
        assert( mnOrdering[iii] >= 0 );
        ordering.push_back( mnOrdering[iii] );
    }

    if( ordering.size() != globalpointer->dimension ){
        throw Error( hscFailedGetOrdering );
    }
}


static void
getEntries( struct nonsinkPath    *pathPointer,
            struct global*         globalpointer,
            std::vector<unsigned>& rows,
            std::vector<unsigned>& cols,
            std::vector<double>&   values )
{
    rows.clear();
    cols.clear();
    values.clear();

    int inEnd1 = pathPointer->inEndIdx;
    int exEnd1 = pathPointer->exEndIdx;
    int exEnd2 = pathPointer->exEndIdx2;
    int currentnode = inEnd1;
    int oldnode = exEnd1;
    double counter2;
    struct edge** adjmatrix = globalpointer->adjmatrix;
    if( adjmatrix[inEnd1]->nodeindex == exEnd1 ){
        counter2 = 1/(adjmatrix[inEnd1]->conductance);
    }
    else counter2 = 1/(adjmatrix[inEnd1]->next->conductance);
    while(1){
        struct edge* edgepointer = adjmatrix[ currentnode ];
        if( edgepointer->nodeindex == oldnode ) edgepointer = edgepointer->next;
        int nextnode = edgepointer->nodeindex;
        
        double diagonal = sqrt( edgepointer->conductance + 1/counter2 );
        rows.push_back( currentnode );
		cols.push_back( currentnode );
		values.push_back( diagonal );
        
        if( exEnd1 < globalpointer->dimension ){
            rows.push_back( currentnode );
            cols.push_back( exEnd1 );
            values.push_back(  -1 / ( counter2 * diagonal ) );
        }
        if( nextnode < globalpointer->dimension ){
            rows.push_back( currentnode );
            cols.push_back( nextnode );
            values.push_back( - edgepointer->conductance / diagonal );
        }
        
        if(nextnode == exEnd2) break;
        counter2 += 1/(edgepointer->conductance);
        oldnode = currentnode;
        currentnode = nextnode;
    }
}


static void
getEntries( struct sinkPath       *pathPointer,
            struct global*         globalpointer,
            std::vector<unsigned>& rows,
            std::vector<unsigned>& cols,
            std::vector<double>&   values )
{
    rows.clear();
    cols.clear();
    values.clear();
    int oldnode = pathPointer->exEndIdx;
    int currentnode = pathPointer->inEndIdx;
    int sink = pathPointer->sinkIdx;
    while(1){
        struct edge* edgepointer = globalpointer->adjmatrix[ currentnode ];
        if( currentnode == sink ){
            assert( edgepointer->conductance > 0 );
            double value = sqrt( edgepointer->conductance );

            rows.push_back( sink );
            cols.push_back( sink );
            values.push_back( value );

            if( edgepointer->nodeindex >= globalpointer->dimension ) return;
            rows.push_back( sink );
            cols.push_back( edgepointer->nodeindex );
            values.push_back( - value );
            return;
        }
        int nextnode=0;
        for(unsigned i=0; i<2; i++){
            if( edgepointer->nodeindex == oldnode ){
                assert( edgepointer->conductance > 0 );
                double value = sqrt( edgepointer->conductance );

                rows.push_back( currentnode );
                cols.push_back( currentnode );
                values.push_back( value );

                if( edgepointer->nodeindex < globalpointer->dimension ){
                    rows.push_back( currentnode );
                    cols.push_back( edgepointer->nodeindex );
                    values.push_back( - value );
                }
            }
            else nextnode = edgepointer->nodeindex;
            edgepointer = edgepointer->next;
        }
        oldnode = currentnode;
        currentnode = nextnode;
    }
}


static void
getEntries( int                    degree3node,
            struct global*         globalpointer,
            std::vector<unsigned>& rows,
            std::vector<unsigned>& cols,
            std::vector<double>&   values )
{
    rows.clear();
    cols.clear();
    values.clear();

    std::vector<unsigned> node(3);
    std::vector<double> cond(3);
    struct edge *edgepointer = globalpointer->adjmatrix[ degree3node ];
    node[0] = edgepointer->nodeindex;
    cond[0] = edgepointer->conductance;
    edgepointer = edgepointer->next;
    node[1] = edgepointer->nodeindex;
    cond[1] = edgepointer->conductance;
    edgepointer = edgepointer->next;
    node[2] = edgepointer->nodeindex;
    cond[2] = edgepointer->conductance;
    double sum = cond[0] + cond[1] + cond[2];
    assert( sum > 0 );
    
    rows.push_back( degree3node );
    cols.push_back( degree3node );
    values.push_back( sqrt( sum ) );
    
    for(unsigned i=0; i<3; i++){
        if( node[i] >= globalpointer->dimension ) continue; 
        rows.push_back( degree3node );
        cols.push_back( node[i] );
        values.push_back( - cond[i] / sqrt( sum ) );
    }
}


void
writePreconditioner( struct global*         globalpointer,
                     std::vector<unsigned>& ordering,
                     std::vector<unsigned>& rowIdx, 
                     std::vector<unsigned>& colIdx,
                     std::vector<double>&   value )
{
    rowIdx.clear();
    colIdx.clear();
    value.clear();

    getOrdering( globalpointer, ordering );
    unsigned dimension = globalpointer->dimension;
    assert( ordering.size() == dimension );
    std::vector<unsigned> inverseOrdering( dimension );
    for(unsigned newIdx = 0; newIdx < dimension; newIdx++){
        assert( ordering[ newIdx ] < dimension );
        inverseOrdering[ ordering[ newIdx ] ] = newIdx;
    }

    struct shrinkEntry* shrinkSequence = globalpointer->shrinkSequence;
    unsigned shrinkEntrynumber = globalpointer->shrinkEntrynumber;
    for(int i=0; i<shrinkEntrynumber; i++){
        std::vector<unsigned> storageRows, storageCols;
        std::vector<double> storageValues;
		switch(shrinkSequence[i].type){
			case 1:
			    getEntries( (struct nonsinkPath *)shrinkSequence[i].pointer,
			                globalpointer,
			                storageRows,
			                storageCols,
			                storageValues );
				break;
			case 0:
			    getEntries( (struct sinkPath *)shrinkSequence[i].pointer,
			                globalpointer,
			                storageRows,
			                storageCols,
			                storageValues );
				break;
			case 2:
			    getEntries( (int)(shrinkSequence[i].pointer),
			                globalpointer,
			                storageRows,
			                storageCols,
			                storageValues );
		}
		assert( storageRows.size() == storageValues.size() );
		assert( storageCols.size() == storageValues.size() );
		for(unsigned ii=0; ii<storageRows.size(); ii++){
		    assert(    inverseOrdering[ storageCols[ii] ]
		            >= inverseOrdering[ storageRows[ii] ] );
            rowIdx.push_back( inverseOrdering[ storageRows[ii] ] );
            colIdx.push_back( inverseOrdering[ storageCols[ii] ] );
            value.push_back( storageValues[ii] );
        }
	}

    int actualsize = globalpointer->actualsize;
    int* mnOrdering = globalpointer->ordering;
    int* appearsizes = globalpointer->appearsizes;
    double* diagonalscale = globalpointer->diagonalscale;
    struct record** mn = globalpointer->mn;
    for(int iii = 0; iii < actualsize; iii++){
        assert( mnOrdering[iii] >= 0 );
        assert( mnOrdering[iii] < dimension );
        unsigned rowNewIdx = inverseOrdering[ mnOrdering[iii] ];
        assert( diagonalscale[ iii ] > 0 );
        double diagonalEntry = 1 / sqrt( diagonalscale[ iii ] );
        rowIdx.push_back( rowNewIdx );
        colIdx.push_back( rowNewIdx );
        value.push_back( diagonalEntry );
        int appearsize = appearsizes[iii];
        for(int j=0; j < appearsize; j++){
            assert( mn[iii][j].nodeindex < iii );
            int oldIdx = mnOrdering[ mn[iii][j].nodeindex ];
            assert( oldIdx >=0 );
            assert( oldIdx < dimension );
            unsigned colNewIdx = inverseOrdering[ oldIdx ];
            assert( colNewIdx > rowNewIdx );
            double entryValue = - mn[iii][j].value * diagonalEntry;
            assert( entryValue < 0 );
            rowIdx.push_back( rowNewIdx );
            colIdx.push_back( colNewIdx );
            value.push_back( entryValue );
        }
    }
}

}

