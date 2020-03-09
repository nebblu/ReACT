/* TODO:
 * - make PyIntegrate intelligent about what dimension of integration is
 *   desired (to avoid having separate function names for 2- and 3-dimensional
 *   integrals) */

%{
#include "Quadrature.h"
%}

/* Typemaps to allow integration of an arbitrary Python function */
%typemap(in) PyObject* pyfunc {
    if(!PyCallable_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object.");
        return NULL;
    }
    $1 = $input;
}

%typemap(typecheck) PyObject* pyfunc {
    $1 = PyCallable_Check($input) ? 1 : 0;
}

%{
static double PythonCallback(PyObject* f, double x) {
    PyObject* arglist = Py_BuildValue("(d)", x);
    PyObject* result = PyEval_CallObject(f, arglist);
    double fx = PyFloat_AsDouble(result);
    Py_DECREF(arglist);
    Py_XDECREF(result);
    return fx;
}

static double PythonCallback2(PyObject* f, double x1, double x2) {
    PyObject* arglist = Py_BuildValue("(dd)", x1, x2);
    PyObject* result = PyEval_CallObject(f, arglist);
    double fx = PyFloat_AsDouble(result);
    Py_DECREF(arglist);
    Py_XDECREF(result);
    return fx;
}

static double PythonCallback3(PyObject* f, double x1, double x2, double x3) {
    PyObject* arglist = Py_BuildValue("(ddd)", x1, x2, x3);
    PyObject* result = PyEval_CallObject(f, arglist);
    double fx = PyFloat_AsDouble(result);
    Py_DECREF(arglist);
    Py_XDECREF(result);
    return fx;
}
%}

%inline %{
#include <boost/bind.hpp>

double PyIntegrate(PyObject* pyfunc, double a, double b, double epsrel = 1e-5, double epsabs = 1e-10) {
    return Integrate(bind(PythonCallback, pyfunc, _1), a, b, epsrel, epsabs);
}

double PyIntegrate2(PyObject* pyfunc, PyObject* pya, PyObject* pyb, double epsrel = 1e-5, double epsabs = 1e-10) {
    double a[2], b[2];
    for(int i = 0; i < 2; i++) {
        PyObject* ai = PySequence_GetItem(pya, i);
        a[i] = PyFloat_AsDouble(ai);
        Py_DECREF(ai);
        PyObject* bi = PySequence_GetItem(pyb, i);
        b[i] = PyFloat_AsDouble(bi);
        Py_DECREF(bi);
    }
    return Integrate<2>(bind(PythonCallback2, pyfunc, _1, _2), a, b, epsrel, epsabs);
}

double PyIntegrate3(PyObject* pyfunc, PyObject* pya, PyObject* pyb, double epsrel = 1e-5, double epsabs = 1e-10) {
    double a[3], b[3];
    for(int i = 0; i < 3; i++) {
        PyObject* ai = PySequence_GetItem(pya, i);
        a[i] = PyFloat_AsDouble(ai);
        Py_DECREF(ai);
        PyObject* bi = PySequence_GetItem(pyb, i);
        b[i] = PyFloat_AsDouble(bi);
        Py_DECREF(bi);
    }
    return Integrate<3>(bind(PythonCallback3, pyfunc, _1, _2, _3), a, b, epsrel, epsabs);
}

%}
