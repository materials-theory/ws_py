ws_py
[![Build Status](https://travis-ci.com/materials-theory/ws_py.svg?branch=master)](https://travis-ci.com/materials-theory/ws_py)
[![Code Health](https://landscape.io/github/materials-theory/ws_py/master/landscape.svg?style=flat)](https://landscape.io/github/materials-theory/ws_py/master)
========
A set of python libraries to handle the VASP input / output files.


Requirements
------------
Based on [python 3](https://python.org), and have external dependencies on the following packages:

* [numpy](http://www.numpy.org)
* [scipy](https://scipy.org)
* [ase](https://wiki.fysik.dtu.dk/ase/)
* [spglib](https://atztogo.github.io/spglib/)
* [argparse](https://docs.python.org/3/library/argparse.html) (for command-line-tools)
* [pandas](https://pandas.pydata.org) (only for "poly" data export)

Due to its home-development environment (in OS X), all shebangs and $PATH might need manual modifications in order to use it properly. Especially, one might need to change **installpath** in `cmdline/vw.py`

All plotting functions are based on [Igor Pro](https://www.wavemetrics.com/products/igorpro/igorpro.htm), so one might need to modify plotting functions to use other plotting programs such as [gnuplot](http://www.gnuplot.info) or [matplotlib](https://matplotlib.org).

License and disclaimers
------------
Python codes inside ws_py package are licensed under the GNU General Public License (GPL) v3.

All responsibilities coming from using this packages are on the users, so please kindly double-check the results and use it with an extra care.

Examples / Tutorials
--------------------
Simple examples and tutorial sets are provided [HERE](https://github.com/materials-theory/ws_py/tree/master/examples).

Functions
------------
If paths are set properly, one might activate the command-line-tool with calling **vw.py**.

Can call `vw.py [function] --help` in terminal to see the usages in line.

----

* ### poly

Tool for handling the polyhedron-related structures and their analysis.

**Usage example :** `vw.py poly -c W -v O [additional options]`

#### Required arguments

  -c [center] : Setting the center atom of the polyhedron. 

  -v [vertex] : Setting the vertex atom of the polyhedron.

#### Options

  -i [input]  : Name of the input file. (Default : POSCAR)

  -o [output] : Name of the output file. Without this, it will not provide additional output file.

  -s [silent] : Silent option. Automatically turns on the -o option.

  -f [face]   : Setting the number of faces in polyhedron. (Default : 8)

  -l [length] : Setting the bond length limit in building a polyhedron. (Default : 3.0)

  -t          : Toggle on the octahedral tilt angle analysis.

#### Sub-options for -t

  -ideal [angle] : Setting the ideal angle in calculating tilt angle variance. (Default : 180.0)

  -centroid      : Use centroid in calculating CVC-O angle. Else, the center atom is used.

  -vo            : Exclude the center atom when calculating the centroid of octahedron.

  -edge          : Include edge-sharing octahedra in the tilt angle analysis.

  -tav           : Change the format of angle analysis in the TAV (tilt angle variance), rather than the average.

---

* ### plot

Tool for plotting (in Igor Pro compatible file formats).

**Usage example :** `vw.py plot [band|dos|pband|temp|diel] [additional options]`

#### Required arguments (But mutually exclusive)

  band  : Plot band structure.

  pband : Plot projected band structure.

  dos   : Plot density-of-states.

  diel  : Plot dielectric responses.

#### Sub-options for "band"

  -o [output]       : Setting the name of the output file. (Default : band.itx)

  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)

  -s                : Toggle on the automatic shifting of VBM to 0

  -g                : Draw the guiding line for high-symmetry points.

  -fake             : Toggle fake-weight method.

#### Sub-options for "dos"

  -o [output]       : Setting the name of the output file. (Default : dos.itx)

  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)

  -c                : Merges every DOS data into a one single .itx file.

#### Sub-options for "pband"

  -o [output]       : Setting the name of the output file. (Default : pband.itx)

  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)

  -s                : Toggle on the automatic shifting of VBM to 0

  -g                : Draw the guiding line for high-symmetry points.

  -fake             : Toggle fake-weight method.

  -atom [ints]      : Atom-projected bandstructure. 

  -orb              : Orbital-projected bandstructure.

  -spin             : Toggle spin(mz)-projected bandstructure mode in spin-polarized calculation.

#### Sub-options for "temp"

  -o [output]       : Setting the name of the output file. (Default : tband.itx)

  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)

  -t [float]        : Temperature in Kelvin scale. (Default: 300K)

  -emin [float]     : Minimum value of eigenvalues. (Default: -20eV)

  -emax [float]     : Maximum value of eigenvalues. (Default: 15eV)

  -egrid [int]      : Number of eigenvalue grids (y-axis). (Default: 2000)

  -etrunc [float]   : Upper and lower truncation energy. (Default: 1eV)

  -kgrid [int]      : Number of k-grids (x-axis). (Default: 800)

  -smin [float]     : Lower limit of smearing window. (Default: -2eV)

  -smax [float]     : Upper limit of smearing window. (Default: 2eV)

  -sgrid [int]      : Number of smearing grids (local x-axis). (Default: 2000)

  -m [string]       : Type of smearing (G, FD, BE, MB available). (Default: Gaussian)

  -fake             : Toggle fake-weight method.
   
#### Sub-options for "diel"

  -o [output]       : Setting the name of the output file. (Default : dielec.itx)

  -d                : Toggle on the transversal and longitudinal direction separation.

  -P                : Selecting the variables to plot in Igor itx file. (eI, eR, n, k, alpha, ELS, R, T, A)

  -D                : Toggle on the manual calculation / inclusion of the Drude peaks.

  -P                : Plasma frequency squared value required in calculating the Drude peaks.

  -t                : Tau (relaxation time) term required in calculating the Drude peaks. (Default : 0.1)

---

* ### conv

Tool for converting file formats. File format converter object is taken from [ASE](https://wiki.fysik.dtu.dk/ase/).

**Usage example :** `vw.py conv [cif2vasp|polyxyz|vasp2qe] -i [input] -o [output]`

#### Required arguments (But mutually exclusive)

  cif2vasp  : Converts cif file to VASP structure format.

  polyxyz   : Extract the octahedron data in VASP and creates xyz files. 

  vasp2qe   : Converts VASP structure file to QE(pwscf) input file format.

#### Common options

  -i [input]  : Name of the input file.

  -o [output] : Name of the output file.
  
#### Sub-options for "polyxyz"

  -i [input]  : Name of the input file.

  -o [output] : Name of the output file.

  -c [center] : Setting the center atom of the polyhedron. 

  -v [vertex] : Setting the vertex atom of the polyhedron.

  -f [face]   : Setting the number of faces in polyhedron. (Default : 8)

  -l [length] : Setting the bond length limit in building a polyhedron. (Default : 3.0)

  -C          : Replaces the center atom position to the centroid of octahedron

  -S          : Toggles the output format in [SHAPE](http://www.ee.ub.edu) program compatible 

---

* ### elec

Tool for handling the electronic structure result, mostly based on the bandstructure calculation.

**Usage example :** `vw.py elec [additional options]`

#### Options

  -fermi [energy] : Manually setting the fermi level.

  -m              : Toggle on the effective mass calculation tool.

  -e              : Toggle on the band edge position tool.

#### Sub-options for "-m" 

  -k              : Setting the k-point number to calculate carrier effective mass.

  -b              : Setting the band number to calculate carrier effective mass.

  -p              : Setting the number of eigenvalue points to interpolate.

---

* ### boltz

Tool for pre / postprocessing files for [BoltzTraP](https://www.sciencedirect.com/science/article/pii/S0010465506001305).

**Usage example :** `vw.py boltz [prep|post] [additional options]`

#### Required arguments (But mutually exclusive)

  prep  : Prepare files required to use BoltzTraP.

  post  : Postprocess BoltzTraP outputs.
  
#### Sub-options for "prep"

  -i [string]       : Name of file to extract structure data (Default : POSCAR)

  -o [string]       : Name of the output file (Default : boltztrap)

  -f [float]        : Fermi energy from scf run (Default : None, automatically reads from OUTCAR)

  -w [float]        : Energy window along Fermi E (Default : 5 eV)

  -g [float]        : Energy grid distance within given window (Default : 0.007 eV)

  -t [integer]      : Maximum temperature to consider (Default : 800 K)

  -step [integer]   : Temperature gradient (Default : 50 K)

  -mu [float]       : Chemical potential window (Default : 2 eV)

  -lpfac [integer]  : Number of lattice points per k-point (Default : 5)

#### Sub-options for "post"

  -o [string]       : Name of the output file (Default : boltzplot.itx)

  -s                : Toggles on the shifting Fermi level to zero

---
