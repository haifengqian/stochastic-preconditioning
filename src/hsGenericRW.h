
#if !defined(hsGenericRW_P)
#define hsGenericRW_P

#include "hsUtility.h"
#include "hsSymRmatCore.h"

namespace hybridSolver {

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
void genericOnSymRaw( struct global*       globalpointer, 
                      std::vector<double>& rhs,
                      double               delta,
                      std::vector<double>& solution );

//******************************************************************************
// a reference implementation of Qian et al. DAC2003
// solve with the given right-hand-side vector b="rhs"
// return solution in "solution"
// accuracy: P[|error| < "delta"] > 99%
// it must not modify any global data
void genericOnSymPreconditioned( struct global*       globalpointer, 
                                 std::vector<double>& rhs,
                                 double               delta,
                                 std::vector<double>& solution );

//******************************************************************************
// return the "entryIndex"-th entry of the solution vector
// accuracy: P[|error| < "delta"] > 99%
double singleEntry( struct global*       globalpointer,
                    std::vector<double>& rhs,
                    double               delta,
                    unsigned             entryIndex );

}

#endif
