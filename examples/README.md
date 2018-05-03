Example set
========
Simplest set of examples/tutorials for new users.

Pre-generated VASP input/output files are included in each folders.

After successful installation of **vw.py**, call `vw.py FUNCTION OPTIONS` following the tutorial provided below.


poly
----
Tool for handling the polyhedron-related structures and their analysis.

* Measuring intra-octahedral distortions in MoO3

`vw.py poly -c Mo -v O -i MoO3.vasp`

Octahedral center atom (Mo) and vertex atom (O) must be set by -c and -v option.

* Measuring inter-octahedral distortions (tilt angle analysis) in MoO3

`vw.py poly -c Mo -v O -i MoO3.vasp -t -s`

Toggled on tilt angle analysis by -t tag and turned on -s silent option to store results at polyresult.csv


elec
----
Tool for handling the electronic structure result, mostly based on the bandstructure calculation.

* Measuring band gap and calculated Fermi energy

`vw.py elec`

* Determining the band edge positions

`vw.py elec -e`

* Calculating effective mass from the bandstructure calculation

`vw.py elec -m -k 81 -b 4 -p 5`

Fermi energy is only the value provided from VASP calculation. One must use it with an extra care when interpreting any scientific data.

Provided band gap value is the result of sampling through the k-grid used in VASP calculation only. 

Carrier effective mass is based on the parabolic E-k relationship, and the convergence test with respect to the grid density near the band edges, and the number of point stencil used, is necessary.


plot
----
Tool for plotting (in Igor Pro compatible *.itx file formats).

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

Tool for converting file formats.

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

  -S          : Toggles the output format in SHAPE program compatible (see http://www.ee.ub.edu)

---

* ### boltz

Tool for pre / postprocessing files for BoltzTraP.

(see https://www.sciencedirect.com/science/article/pii/S0010465506001305)

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
