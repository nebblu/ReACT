%{
#include "Common.h"
%}

/* Type definitions */
typedef double real;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

/* Physical constants in SI units */
%nodefaultctor Constants;
%nodefaultdtor Constants;
struct Constants {
    static const real second = 1;
    static const real t_100 = 3.08568e17;      // 1/(100 km/sec/Mpc) in seconds
    static const real meter = 1;
    static const real km = 1000 * meter;
    static const real Mpc = 3.08568025e19 * km;

    static const real c = 299792458;           // speed of light
    static const real G = 6.67428e-11;         // Newton's constant
    static const real hbar = 1.054571628e-34;  // Planck's constant (over 2*pi)
    static const real k = 1.3806504e-23;       // Boltzmann constant
};


/* Allow any float-able type to be treated as a real */
%typemap(in) real {
    PyObject* floatobj = PyNumber_Float($input);
    $1 = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
}
%typemap(typecheck) real {
    PyObject* floatobj = PyNumber_Float($input);
    if(floatobj != NULL) {
        $1 = 1;
        Py_DECREF(floatobj);
    }
    else
        $1 = 0;
}
