%module pycopter

%{
#define SWIG_FILE_WITH_INIT
%}

%naturalvar;

%include "numpy.i"
%include "typemaps.i"

%init %{
    import_array();
%}

%include "Common.i"
%include "array.i"
%include "Datafile.i"
%include "Quadrature.i"
%include "Spline.i"
%include "Timer.i"

%include "Cosmology.i"
%include "GrowthFunction.i"
%include "CorrelationFunction.i"
%include "PowerSpectrum.i"
%include "Closure.i"
%include "ConsistentSPT.i"
%include "InterpolatedPS.i"
%include "LagrangianResummation.i"
%include "LargeN.i"
%include "LinearPS.i"
%include "NoWigglePS.i"
%include "RPT.i"
%include "SPT.i"
%include "SimpleRG.i"
