# ReACT 
# This branch is aimed at implementing massive neutrinos (1909.02561) into the current ReACT code. It is currently undergoing development and validation. 

## Introduction

ReACT is an extension of the software package Copter (0905:0479) and MG-Copter (1606.02520) which allows for 
the calculation of redshift and real space large scale structure 
observables for a wide class of gravity and dark energy models. 

Additions to original Copter code http://mwhite.berkeley.edu/Copter/: 

* Spherical collapse in modified gravity (1812.05594): `reactions/src/SCOL.cpp`

* Halo model power spectrum for general theories (1812.05594):  `reactions/src/HALO.cpp`

* Real and redshift space LSS 2 point statistics for modified gravity and dark energy (1606.02520): `reactions/src/SPT.cpp`

* Numerical perturbation theory kernel solvers (1606.02168): `reactions/src/SpericalFunctions.cpp`

* Real space bispectra in modified gravity (1808.01120): `reactions/src/BSPT.cpp`

* Numerical perturbation theory kernel solver up to 4th order for 1-loop bispectrum (1808.01120): `reactions/src/BSPTN.cpp`

## Requirements
### C++ Compiler and automake
ReACT is written in C++, so you'll need a relatively modern C++ compiler.
It also uses automake.

### GSL
The ReACT extension makes use of a number of [gsl](http://www.gnu.org/software/gsl/) packages such as odeiv2 package which is part of the gsl library.  You will need to have a version of gsl installed (gsl 2. or later). This is needed for the numerical calculation of the perturbation theory kernels.

### SUNDIALS

For Spherical Collapse module (`reactions/src/SCOL.cpp`) you will also need the [SUNDIALS](https://computing.llnl.gov/projects/sundials) package version 4.0.

Get these packages using your package manager of choice (e.g., `homebrew` on mac OS).

One should also make sure that sundials is on your library path, for example 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/install_dir/lib64:${LD_LIBRARY_PATH} 

## Installation
The Python interface can then be installed with
```
$ python setup.py install
```
or 
```
$ pip install .
```
If you want to work on the library, install it in "developer" mode:
```
$ python setup.py develop
```
or
```
$ pip install -e .
```

## Choosing a model of gravity or dark energy

To choose a model of gravity or dark energy within the framework described in `arXiv:1606.02520` for example, open 
the `SpecialFunction.cpp` file in the `reactions/src` directory. Towards the top of the file you will 
find the background Hubble and mu, gamma2, gamma3 functions as well as the modified spherical collapse function `F`.
As examples, the Hu Sawicki, nDGP and GR versions of these functions have been included. wCDM background Hubble rates have also been included as presets. Simply edit in your version of these functions and edit-out the unwanted ones. Then just re-install the package as described above. Note that you may need to run `make clean` in the pyreact or reactions directory before reinstalling. 


### Adding in models

One can add in new models. As a default, a maximum of 3 theory parameters have been allowed for. 
There are instructions within the SpecialFunction.cpp file to include more than 
3 theory parameters. The entire file is heavily commented so that 
making any additional edits shouldn't be (too) difficult. Note there may be some issues with dependencies when adding more than 3 theory parameters - this needs to be tested.  


## Running ReACT

### Python
An example jupyter notebook that demonstrates the usage of ReACT can be found in `notebooks/pyreact_demo.ipynb`.

### Stand-alone (ReACT and MG-Copter) 
One can also run ReACT and MG-Copter for specific cosmologies with specified transfer functions. A number of example output c++ scripts have been included in `reactions/examples` as well as a number of cosmologies in `reactions/examples/transfers`.

In particular, the bs.cpp, spt.cpp, halo_ps.cpp examples compute various quantities like the 1-loop bispectrum, halo-model quantities and redshift space power spectrum (TNS). Some example cosmologies are found in examples/transfers. One can produce their own cosmology by getting the transfer function from CAMB say, normalising it to 1 at small k, and also constructing a cosmology .ini file as in the examples. 

We can compile these examples with a command similar to : 

> gcc -I/Users/bbose/Desktop/ReACT-master/reactions/include -L/Users/bbose/Desktop/ReACT-master/reactions/lib -lcopter -lgsl -lstdc++ bs.cpp -o test

Then just run 

> ./test 


## Citation

When using ReACT in a publication, please acknowledge the code by citing the following papers:

arXiv:1812.05594 : "On the road to percent accuracy: non-linear reaction of the matter power spectrum to dark energy and modified gravity"

arXiv:2005.12184 : "On the road to per-cent accuracy IV:  ReACT -- computing the non-linear power spectrum beyond LCDM"

Respective bibtex entries:

```
 @article{Cataneo:2018cic,
   author = "Cataneo, Matteo and Lombriser, Lucas and Heymans, Catherine and Mead, Alexander and Barreira, Alexandre and Bose, Sownak and Li, Baojiu",
    title = "{On the road to percent accuracy: non-linear reaction of the matter power spectrum to dark energy and modified gravity}",
     eprint = "1812.05594",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stz1836",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "488",
    number = "2",
    pages = "2121--2142",
    year = "2019"
}
```

```
@article{Bose:2020wch,
    author = "Bose, Benjamin and Cataneo, Matteo and Tröster, Tilman and Xia, Qianli and Heymans, Catherine and Lombriser, Lucas",
    title = "{On the road to per-cent accuracy IV: ReACT -- computing the non-linear power spectrum beyond $\Lambda$CDM}",
    eprint = "2005.12184",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "5",
    year = "2020"
}
```

## Notes on parameter ranges (28/05/20)
* To optimise root finding within the spherical collapse portion of the code, the maximum redshift that one can solve the reaction for currently is z=2.5. 
* There are some current issues in the wCDM part of the code. Namely for very particular values of w0 and wa in the CPL evolving dark energy case, the spherical collapse library cannot solve the virial theorem. We advise sticking to the ranges 
-1.3<w0<-0.7 and -1.5<wa<0.6 to avoid these issues. 
* Spherical collapse will not solve for values of sigma_8(z=0)< 0.55 or 1.4<sigma_8(z=0).  

## Miscellaneous notes (28/05/20)
* One may need to add the sundials library directory to LDFLAGS in pyreact/Makefile for installation to complete correctly:
> LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L/home/bose/sundials/install_dir/lib64
* One may also need to add the sundials include directory as a CPPFLAG in pyreact/Makefile for installation to complete correctly:
> CPPFLAGS += -I/home/bose/sundials/install_dir/include
* If errors in spherical collapse are experienced for non-f(R) theories, try setting yenvf=0 in the scol_init function in reactions/src/HALO.cpp.
* Note if using the stand-alone version of ReACT, the reaction may have deviations away from unity of the order of ~0.1-0.3% for k<1e-3. Pyreact automatically sets it to unity at such large scales. 

## Additional Libraries from old versions of MG-Copter 

There are additional libraries which you can (try to!) add in in the reactions/src/extra_libraries folder:

Clustering dark energy (example: 1611.00542):
CDE.cpp

Gaussian streaming model for modified gravity (1705.09181):
CorrelationFunction.cpp

Heiarchichal moments in modified gravity (1703.03395):
HigherStat.cpp

CMB lensing statistics (1812.10635):
CMBc.cpp

Regularised perturbation theory (1408.4232):
RegPT.cpp

Effective field theory of LSS with resummation (1802.01566):
EFT.cpp

Works in progress (not tested):
Interpolator in modified gravity theory param space:
PMG_Interpolation.cpp 

Note that the SpecialFunctions.cpp and SPT.cpp libraries in the extra_libraries folder
contain some additional functions which are necessary for some of these extra libraries. Please contact benjamin.bose@unige.ch if you are interested in implementing these and encounter trouble in doing so. 
