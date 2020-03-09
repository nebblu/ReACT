EXTRA_DIST = \
    Closure.i Common.i ConsistentSPT.i CorrelationFunction.i Cosmology.i Datafile.i \
    GrowthFunction.i InterpolatedPS.i LagrangianResummation.i LargeN.i LinearPS.i \
    NoWigglePS.i PowerSpectrum.i Quadrature.i RPT.i SPT.i SimpleRG.i Spline.i Timer.i \
    array.i numpy.i pycopter.i
EXTRA_PYTHON = pycopter.py
EXTRA_LTLIBRARIES = _pycopter.la

nodist_python_PYTHON = $(PYCOPTER_PY)
pyexec_LTLIBRARIES = $(_PYCOPTER_LA)

nodist__pycopter_la_SOURCES = pycopter_wrap.cpp
_pycopter_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_pycopter_la_LDFLAGS = -shared -module -avoid-version -rpath $(pyexecdir)
_pycopter_la_LIBADD = ../src/libcopter.la

pycopter_wrap.cpp: pycopter.i
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src -o $@ $<

CLEANFILES = $(_PYCOPTER_LA) $(PYCOPTER_PY) pycopter_wrap.cpp
