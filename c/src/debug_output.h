#pragma once

#include <iostream>

/*
 * Debugging output -- should be used for any information that helps in
 * debugging (it's ok if it slows things down).
 */

#ifdef DEBUG_OUTPUT
#define DEBUG_LOCATION "(" << __FILE__ << ":" << __LINE__ << "):"
#define DMSG(msg) std::cerr << DEBUG_LOCATION << " " << msg << std::endl;
#define DOUT(x) std::cerr << DEBUG_LOCATION << " " << x;
#else
#define DMSG(msg) (void)0
#define DOUT(x) (void)0
#endif
