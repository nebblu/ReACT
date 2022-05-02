
# ReACT with massive neutrinos
## This branch is a first implementation of massive neutrinos (1909.02561) into the basic ReACT code. 

## Table of contents: 

* [Introduction](https://github.com/nebblu/ReACT#introduction)
* [Requirements](https://github.com/nebblu/ReACT#requirements)
* [Installation](https://github.com/nebblu/ReACT#installation)
* [Models of gravity and dark energy](https://github.com/nebblu/ReACT#models-of-gravity-and-dark-energy)
* [Running ReACT](https://github.com/nebblu/ReACT#running-react)
* [Massive neutrinos](https://github.com/nebblu/ReACT#extended-react-massive-neutrinos-with-modified-gravity-andor-evolving-dark-energy)
* [Citation](https://github.com/nebblu/ReACT#citation)
* [Notes on parameter ranges](https://github.com/nebblu/ReACT#notes-on-parameter-ranges-updated-220321)
* [Miscellaneous](https://github.com/nebblu/ReACT#miscellaneous-notes-280520)
* [Additional Libraries](https://github.com/nebblu/ReACT#additional-libraries-from-old-versions-of-mg-copter)
* [Linux specific installtion](https://github.com/nebblu/ReACT#linux-specific-installation)

## Introduction

ReACT is an extension of the software package Copter (0905.0479) and MG-Copter (1606.02520) which allows for the calculation of redshift and real space large scale structure observables for a wide class of gravity and dark energy models. 

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

For Spherical Collapse module (`reactions/src/SCOL.cpp`) you will also need the [SUNDIALS](https://computing.llnl.gov/projects/sundials) package version 4.1.0 (also tested to work with version 5.0) **Note** that later versions of sundials may produce errors as they update their functions so please work with version 4.1.0 

Get these packages using your package manager of choice (e.g., `homebrew` on mac OS).

One should also make sure that sundials is on your library path, for example 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/install_dir/lib64:${LD_LIBRARY_PATH} 

### Python
A recent version of python should be installed. I have worked with Python 3.6.8 without issue. 


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

**See [here](https://github.com/nebblu/ReACT#linux-specific-installation) for a linux specific installation guide** 

## Models of gravity and dark energy

*Update (5/04/2021)* The parameter 'model' selects the model of gravity or dark energy to be assumed. This is passed to the 'initialise' function in the case of halo model reaction calculations (or to any other relevant function). The value of this parameter dictates which model: 

In Pyreact we specify the model as

gr : general relativity | f(r) : Hu-Sawicki f(R) | dgp : normal branch of DGP gravity | quintessence : Quintessence | cpl : CPL evolving dark energy

Model parameters are none, fR0, Omega_rc, w , {w,wa} respectively. 

In C++ code this is specified as an integer: 

1: GR | 2: Hu-Sawicki f(R) | 3: nDGP | 4: Quintessence | 5: CPL | 

Note that 1,2,3 all assume a LCDM background expansion. 

*Deprecated (5/04/2021)* :To choose a model of gravity or dark energy within the framework described in `arXiv:1606.02520` for example, open 
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

Note that the transfer function should have two columns: k[h/Mpc] and T(k). The supplied transfer function should also be a LCDM transfer function at z=0. MG-Copter then renormalises the P(k) using internally computed LCDM and MG/DE growth factors. 

We can compile these examples with a command similar to : 

> gcc -I/Users/bbose/Desktop/ReACT-master/reactions/include -L/Users/bbose/Desktop/ReACT-master/reactions/lib -lcopter -lgsl -lstdc++ bs.cpp -o test

Then just run 

> ./test 


## Extended ReACT: massive neutrinos with modified gravity and/or evolving dark energy. 

We have made extensions to the ReACT framework to include the effects of massive neutrinos, combining [1909.02561](https://arxiv.org/abs/1909.02561) and [1812.05594](https://arxiv.org/abs/1812.05594).  To accommodate massive neutrino effects, ReACT now has the option to take the 'real' cosmology's transfer function at the target redshift as produced by MGCAMB (as opposed to the LCDM transfer at z=0). The 1-loop perturbation theory part of the reaction is then approximated by neglecting beyond linear order massive neutrino effects as described in [1902.10692](https://arxiv.org/abs/1902.10692).

**Note on transfer function**: For the example file reactions/examples/reaction_mnu.cpp, the transfer function given to ReACT should have column 6,7,8 correspond to massive neutrinos, total matter and cdm+baryons respectively. This file should be headerless. 

We have introduced two new flags when computing the reaction. These are then fed to a global initialiser '$initialise$' that calculates linear growth rates, performs modified and unmodified spherical collapse, as well as calculation of $k_\star$ and $\mathcal{E}$.   

These flags are evident in the new example file reactions/examples/reactions_mnu.cpp which computes the wCDM + massive neutrino reaction and halofit pseudo spectrum for a target cosmology and transfer function. Specifically, the flags perform the following roles

**modg**: Tells ReACT to manually set $k_\star$ and $\mathcal{E}$ to LCDM values (1e-6 and 1. resp.). This is needed because of sensitivity of $\mathcal{E}$ to the ratio of 1-halo terms which may not be exactly equal at large scales for different cosmologies even when modified gravity isn't present. 

**modcamb**:  Tells ReACT whether or not to treat the input transfer function file as the full real transfer function at the target redshift (as produced by MGCAMB for example) - true value - or as a z=0, LCDM transfer (two columns k[h/Mpc] and T(k) normalised to 1 at small k) - false value. 

Note that the MGCAMB produced transfer should include columns for total matter (col 7), cdm + baryons (col 8) and massive neutrinos (col6). It does not need to be normalised to 1. at small k. There should be no header line in the transfer function file so this may need to be removed manually ( a segmentation fault will be thrown if it is not removed). 

If these flags are not specified, ReACT will assume a LCDM, z=0 transfer is being fed  (mgcamb=false) with modified gravity active (modg = true), as in original version of the code, which is compatible with pyreact and the cosmoSIS module.  

**Note** the cosmoSIS module has not yet been extended to include massive neutrinos. 

**Note** Pyreact currently only accepts three model parameters. If you wish to upgrade to more input parameters, you will need to edit pyreact/react_wrapper.cpp and pyreact/react.py. Note that to keep the react cosmosis module functional, you will also need to make the relevant adjustments in the cosmosis folder. 

## Citation

When using ReACT in a publication, please acknowledge the code by citing the following papers:

arXiv:1812.05594 : "On the road to percent accuracy: non-linear reaction of the matter power spectrum to dark energy and modified gravity"

arXiv:2005.12184 : "On the road to per-cent accuracy IV:  ReACT -- computing the non-linear power spectrum beyond LCDM"

arXiv:2105.12114 : "On the road to percent accuracy V: the non-linear power spectrum beyond LCDM with massive neutrinos and baryonic feedback"

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

```
@article{Bose:2021mkz,
    author = "Bose, Benjamin and Wright, Bill S. and Cataneo, Matteo and Pourtsidou, Alkistis and Giocoli, Carlo and Lombriser, Lucas and McCarthy, Ian G. and Baldi, Marco and Pfeifer, Simon and Xia, Qianli",
    title = "{On the road to percent accuracy V: the non-linear power spectrum beyond $\Lambda$CDM with massive neutrinos and baryonic feedback}",
    eprint = "2105.12114",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "5",
    year = "2021"
}
```


## Notes on parameter ranges (Updated: 22/03/21)
* To optimise root finding within the spherical collapse portion of the code, the maximum redshift that one can solve the reaction for currently is z=2.5. 
* There are some current issues in the wCDM part of the code. Namely for very particular values of w0 and wa in the CPL evolving dark energy case, the spherical collapse library cannot solve the virial theorem. We advise sticking to the ranges 
-1.3<w0<-0.7 and -1.5<wa<0.6 to avoid these issues. 
* Spherical collapse will not solve for values of sigma_8(z=0)< 0.55 or 1.4<sigma_8(z=0).  
* We have tested the massive neutrino code for m_nu=0.48eV (Omega_nu = 0.0105) without issue in the absence of modified gravity or evolving dark energy. 
* Spherical collapse may encounter issues for high fr0 values in combination with high Omega_nu. We have tested the code for Log[fr0]=-4 with m_nu=0.3eV with success but do not recommend higher values than these in combination.  

## Miscellaneous notes (28/05/20)
* One may need to add the sundials library directory to LDFLAGS in pyreact/Makefile for installation to complete correctly:
> LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L/home/bose/sundials/install_dir/lib64
* One may also need to add the sundials include directory as a CPPFLAG in pyreact/Makefile for installation to complete correctly:
> CPPFLAGS += -I/home/bose/sundials/install_dir/include
* If errors in spherical collapse are experienced for non-f(R) theories, try setting yenvf=0 in the scol_init function in reactions/src/HALO.cpp.
* Note if using the stand-alone version of ReACT, the reaction may have deviations away from unity of the order of ~0.1-0.3% for k<1e-3. Pyreact automatically sets it to unity at such large scales. 

## Additional Libraries (from old versions of MG-Copter) 

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



## Linux specific installation 

**In all the following change the paths accordingly**

1) Make sure the following are installed: python, sundials, g++, gsl. The versions I've tested with are:

Python 3.6.8

g++ (GCC) 8.5.0

sundials-4.1.0

2) Clone react_with_neutrinos branch:

```
$ git clone -b react_with_neutrinos git@github.com:nebblu/ReACT.git
```

3) Add in sundials directory to pyreact/Makefile :

**Example:** 


```
LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L/home/bose/sundials/instdir/lib64
```
```
CPPFLAGS += -I/home/bose/sundials/instdir/include
```
```
cd $(COPTER_DIR) && CXXFLAGS="$(CXXFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)" ./build_copter.sh
```

4) Export sundials library to LD_LIBRARY_PATH:

```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/instdir/lib64:${LD_LIBRARY_PATH}
```

5) install react

```
$ python3 setup.py develop --user
```

6) Check that everything is working. Go to the reactions/examples/ directory and try to run one of the example files. For example:

```
$ g++ -I/home/bose/react_tutorial/ReACT/reactions/include -L/home/bose/react_tutorial/ReACT/reactions/lib -lcopter -lgsl -lstdc++ halo_ps.cpp -o test
```

```
$ ./test
```

Example output:

```
$ 0 1.000000e-03 4.039658e+02 4.213977e+02 9.959108e-01 ...
```

**Note** you may also need to export the copter library path to run the example files: 

```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/react_tutorial/ReACT/reactions/lib:${LD_LIBRARY_PATH}
```
