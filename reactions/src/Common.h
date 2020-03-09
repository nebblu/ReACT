#ifndef COMMON_H
#define COMMON_H

#ifndef NULL
#define NULL 0
#endif

#include <cmath>

/* Type definitions */
typedef double real;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

/* Forward declarations of generic classes */
class CorrelationFunction;
class Cosmology;
class Datafile;
class PowerSpectrum;
class Spline;
class Timer;
class array;
class pstring;

/* Convenient print functions (that flush the buffer afterwards) */
//   Print to stdout
void info(const char* format, ...);
//   Print to stdout if VERBOSE is defined.
void verbose(const char* format, ...);
//   Print to stdout if DEBUG is defined.
void debug(const char* format, ...);
//   Print to stderr.
void warning(const char* format, ...);
//   Print to stderr and abort.
void error(const char* format, ...);

/* Physical constants in SI units */
namespace Constants {
    const real second = 1;
    const real meter = 1;
    const real km = 1000 * meter;
    const real Mpc = 3.08568025e19 * km;

    const real c = 299792458;           // speed of light
    const real G = 6.67428e-11;         // Newton's constant
    const real hbar = 1.054571628e-34;  // Planck's constant (over 2*pi)
    const real k = 1.3806504e-23;       // Boltzmann constant

    const real t_100 = 3.08568e17;      // 1/(100 km/sec/Mpc) in seconds
}

/***** Math routines *****/

/* Small integer powers */
static inline real pow2(real x) { return x*x; }
static inline real pow3(real x) { return x*x*x; }
static inline real pow4(real x) { return pow2(pow2(x)); }
static inline real pow5(real x) { return x*pow4(x); }
static inline real pow6(real x) { return pow3(pow2(x)); }
static inline real pow7(real x) { return x*pow6(x); }
static inline real pow8(real x) { return pow2(pow4(x)); }

#endif // COMMON_H
