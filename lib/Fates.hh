// Fates.hh
// T. M. Kelley
// Oct 16, 2015
// (c) Copyright 2015 LANSLLC, all rights reserved


#ifndef FATES_HH
#define FATES_HH

#include <ostream>

namespace nut
{
    enum struct Fates
    {
        NUCLEON_ABS = 0,
        ESCAPE = 1,
        STEP_END = 2,
        NUM_VALID_FATES,
        NOT_DEAD_YET = 99,
        UNHANDLED_EVENT = 998,    // Couldn't apply apply event
        UNDETERMINED_EVENT = 999  // Couldn't decide event
    };

} // nut

#endif // include guard


// End of file
