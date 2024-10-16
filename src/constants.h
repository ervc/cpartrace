#include <math.h>

// Swithces for data
#define RHO    1
#define VPHI   2
#define VR     3
#define VTHETA 4
#define VX     5
#define VY     6
#define VZ     7
#define DRHODX 8
#define DRHODY 9
#define DRHODZ 10

// constants (cgs)
#define MSUN    ( 1.9891e33 )
#define MEARTH  ( 5.97e27 )
#define AU      ( 1.49597871e13 )
#define G       ( 6.674e-8 )
#define YR      ( 3.1557600e7 )

// derived constants
#define R0          ( 5.2 * AU )
#define TIME        ( sqrt(R0*R0*R0/G/MSUN) )

// Read in from FARGO or JUPITER output
#define FARGO_MODEL   1
#define JUPITER_MODEL 2
