EXTRA_DIST = 

lib_LTLIBRARIES = libcopter.la
libcopter_la_CPPFLAGS = -DDATADIR='"$(datadir)/Copter"'
libcopter_la_CXXFLAGS = $(OPENMP_CFLAGS)
libcopter_la_LDFLAGS = $(OPENMP_CFLAGS)
libcopter_la_SOURCES = \
    Closure.cpp Common.cpp ConsistentSPT.cpp CorrelationFunction.cpp \
    Cosmology.cpp Datafile.cpp FlowingWithTime.cpp GrowthFunction.cpp \
    HigherSPT.cpp InterpolatedPS.cpp HigherStat.cpp Kaiser.cpp \
    LagrangianResummation.cpp LargeN.cpp \
    LinearPS.cpp MonteCarlo.cpp NoWigglePS.cpp PowerSpectrum.cpp \
    Quadrature.cpp RPT.cpp SPT.cpp SimpleRG.cpp SpecialFunctions.cpp \
    Spline.cpp Timer.cpp \
    array.cpp cfg.c pstring.cpp rng.c RegPT.cpp PMG_Interpolation.cpp EFT.cpp BSPT.cpp BSTPN.cpp CDE.cpp CMBc.cpp SCOL.cpp HALO.cpp

copterincdir = $(includedir)/Copter
copterinc_HEADERS = \
    Closure.h Common.h ConsistentSPT.h CorrelationFunction.h \
    Cosmology.h Datafile.h FlowingWithTime.h GrowthFunction.h \
    HigherSPT.h InterpolatedPS.h HigherStat.h Kaiser.h \
    LagrangianResummation.h LargeN.h \
    LinearPS.h MonteCarlo.h NoWigglePS.h ODE.h ODE.inl PowerSpectrum.h \
    Quadrature.h Quadrature.inl RPT.h SPT.h SimpleRG.h SpecialFunctions.h \
    Spline.h Timer.h \
    array.h cfg.h pstring.h rng.h RegPT.h PMG_Interpolation.h EFT.h BSPT.h BSPTN.h  CDE.h CMBc.h SCOL.h HALO.h
