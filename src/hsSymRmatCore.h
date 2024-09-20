
#if !defined(hsSymRmatCore_P)
#define hsSymRmatCore_P

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>

#define UNKNOWN -235
#define WALKSTEPLIMIT 50000

#define SWITCHTHRESHOLD 100
#define RANDOMSEED -1       

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

namespace hybridSolver {

struct edge {
    int nodeindex, probability;
    double conductance;
    struct edge *next;
};

struct edge2 {
    int nodeindex, probability;
    double prob;
};

struct sinkPath {
    double load;
    int inEndIdx, exEndIdx, sinkIdx;   
};

struct nonsinkPath {
    double res, load1, load2;    
    int inEndIdx, exEndIdx, inEndIdx2, exEndIdx2;
};

struct shrink {
    void* pointer;
    int   type; 
    struct shrink *next;
};

struct shrinkEntry {
    void* pointer;
    int   type; 
};

struct record {
    int    nodeindex;
    double value;
};

struct global {
    int dimension; 
    int realnodenum; 
    int actualsize; 
    double targetnorm; 
    struct edge **adjmatrix; 
    struct edge2 **adjmatrix2; 
    int *nodedegree;
    double *rhs; 
    double *newrhs; 
    double *diagonal; 
    struct edge *homeregistry; 
    
    int *shrinkflag; 
    struct shrinkEntry *shrinkSequence; 
    int shrinkEntrynumber; 
    int *ordering;
    struct record **mn;
    int *appearsizes;
    struct record **nm;
    int *nmappearsizes;
    struct record **terminalmn;
    int *terminalmnappearsizes;
    double *diagonalscale;
    double *solution;
};

double seconds();
void freeedgepointer( struct edge *edgepointer );
void disorder( int *ports2,
               int portnumber2);
void marksinkPath( struct edge **adjmatrix,
                   int          *nodedegree,
                   int          *shrinkflag,
                   int          *port1addrs,
                   int          *port1exaddrs,
                   int           currentnode,
                   int           dimension);
void marknonsinkPath( struct edge         **adjmatrix,
                      int                  *nodedegree,
                      int                  *shrinkflag,
                      struct nonsinkPath   *pathPointer,
                      int                   i,
                      int                   dimension );
void shrinksinkPath( struct edge      **adjmatrix,
                     int               *nodedegree,
                     struct sinkPath   *pathPointer,
                     int                dimension );
void rhsshrinksinkPath( struct global     *globalpointer,
                        struct sinkPath   *pathPointer );
void shrinknonsinkPath( struct edge         **adjmatrix,
                        int                  *nodedegree,
                        struct nonsinkPath   *pathPointer,
                        int                   dimension );
void rhsshrinknonsinkPath( struct global      *globalpointer,
                           struct nonsinkPath *nonsinkPathPointer );
void recoversinkPath( struct edge      **adjmatrix,
                      double            *currentload,
                      double            *voltage,
                      struct sinkPath   *pathPointer);
void recovernonsinkPath( struct edge         **adjmatrix,
                         double               *currentload,
                         double               *voltage,
                         struct nonsinkPath   *pathPointer);
void addedge( int           fromnodeindex,
              int           tonodeindex,
              double        conductance,
              struct edge **adjmatrix,
              int          *nodedegree,
              int           dimension );
struct shrink *shrinkCall1( struct edge  **adjmatrix,
                            int           *nodedegree,
                            struct shrink *shrinkSequence,
                            int            dimension,
                            int           *shrinkflag);
struct shrink *shrinkCall2( struct edge  **adjmatrix,
                            int           *nodedegree,
                            struct shrink *shrinkSequence,
                            int            dimension,
                            int           *shrinkflag);
struct shrink *shrinkCall3( int           *ordering,
                            struct edge  **adjmatrix,
                            int           *nodedegree,
                            struct shrink *shrinkSequence,
                            int            dimension,
                            int           *shrinkflag );
struct shrink *processSinkPath( int            sinkIdx,
                                struct edge  **adjmatrix,
                                int           *nodedegree,
                                struct shrink *shrinkSequence,
                                int           *shrinkflag,
                                int            dimension );
struct shrink *processNonSinkPath( int            nodeindex,
                                   struct edge  **adjmatrix,
                                   int           *nodedegree,
                                   struct shrink *shrinkSequence,
                                   int           *shrinkflag,
                                   int            dimension );
void rhsStarShrink( struct global *globalpointer,
                    int            nodeindex );
void starRecover( int           nodeindex,
                  struct edge **adjmatrix,
                  double       *currentload,
                  double       *voltage );
void cutonedirectionedge( int           fromnodeindex,
                          int           tonodeindex,
                          struct edge **adjmatrix,
                          int          *nodedegree );
void readinput( char *,
                struct global * );
void shrinkGraph( struct global * );
void preprocessing( struct global * );
void postprocessing( struct global * );
void randomwalk( struct global *,
                 double );
void processrhs( struct global *,
                 double );
void solve( struct global *,
            double );
void freespace( struct global * globalpointer );
void writePreconditioner( struct global*         globalpointer,
                          std::vector<unsigned>& ordering,
                          std::vector<unsigned>& rowIdx, 
                          std::vector<unsigned>& colIdx,
                          std::vector<double>&   value );

}

#endif

