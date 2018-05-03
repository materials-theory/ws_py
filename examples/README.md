Example set
========
Simple set of examples/tutorials for new users.

Pre-generated VASP input/output files are included in each folders.

After successful installation of **vw.py**, call `vw.py FUNCTION OPTIONS` following the tutorial provided below.


poly
----
Tool for handling the polyhedron-related structures and their analysis.

* Measuring intra-octahedral distortions of MoO<sub>3</sub> (/examples/poly)

    `vw.py poly -c Mo -v O -i MoO3.vasp`

    Octahedral center atom (Mo) and vertex atom (O) must be set by -c and -v option.

* Measuring inter-octahedral distortions (tilt angle analysis) of MoO<sub>3</sub> (/examples/poly)

    `vw.py poly -c Mo -v O -i MoO3.vasp -t -s`

    Toggled on tilt angle analysis by -t tag and turned on -s silent option to store results at polyresult.csv

plot
----
Tool for plotting (in Igor Pro compatible *.itx file formats).

**Every plotting presets (e.g. style, layout, color scale, ...) can be manually adjusted by modifying `Plotter/plotter.py` and modules therein.**

* Plotting normal bandstructure of Si bulk (/examples/plot/0_band_normal)

    `vw.py plot band -g -s`

    High-symmetry point guide line is toggled on by -g tag, and VBM is shifted to 0 eV by -s tag.

* Plotting fake-weight bandstructure (e.g. HSE06 calculations) of Si bulk (/examples/plot/1_band_fake_weight)

    `vw.py plot band -g -s -fake`

* Plotting density-of-states of Si bulk (/examples/plot/2_dos)

    `vw.py plot dos -c`

    In order to generate only one .itx file, -c tag is used. Otherwise, it will be seperated out into a number of *.itx files.

* Plotting dielectric function of Si bulk (/examples/plot/3_dielectric)

    `vw.py plot diel -p eI eR`

    Plotting imaginary (eI) and real (eR) part of dielectric function. Other available plotting options are: n, k, alpha, ELS, R, T, and A.

* Plotting direction-dependent dielectric function of MoS<sub>2</sub> (/examples/plot/4_dielectric_direction)

    `vw.py plot diel -p eI eR -d`

    If one want to separate out the direction-dependent values (i.e. transversal and longitudinal), -d tag can be used.

* Plotting dielectric function of Al bulk with Drude term correction (/examples/plot/5_dielectric_drude)

    `vw.py plot diel -D -P 65.341`

    Plasma frequency squared value can be obtained from the OUTCAR file. 

* Plotting (orbital-, atom-)projected bandstructure of Sb<sub>2</sub>Te<sub>3</sub> (/examples/plot/6_pband)

    `vw.py plot pband -g -s -orb`

    Provides the orbital-projected bandstructure (s, p, d, and total).

    `vw.py plot pband -g -s -atom 0-7` 

    Provides the atom-projected bandstructure of Sb. Please be aware of python-style numbering -- the number starts from 0, not 1.

    Multiple atoms can be plotted at once, by using -atom tag like `-atom 0 9`, and `-atom 0-7 8-19`.

    Also note that -orb and -atom option can be used simultaneously. 

* Plotting spin-projected bandstructure of graphene (/examples/plot/7_pband_spin)

    `vw.py plot pband -g -s -spin`

    Red color is set for + occupation, while blue color is set for - occupation in default.

* Plotting temperature-broadened bandstructure with statistical models. (/examples/plot/8_temp)

    `vw.py plot temp -t 600 -m FD`

    This will generate the broadened bandstructure at 600K with Fermi-Dirac distribution applied.

elec
----
Tool for handling the electronic structure result, mostly based on the bandstructure calculation.

* Measuring band gap and calculated Fermi energy of bulk Si (/examples/elec)

    `vw.py elec`

* Determining the band edge positions of bulk Si (/examples/elec)

    `vw.py elec -e`

* Calculating hole effective mass of bulk Si from the bandstructure calculation (/examples/elec)

    `vw.py elec -m -k 81 -b 4 -p 5`

    Fermi energy is only the value provided from VASP calculation. One must use it with an extra care when interpreting any scientific data.

    Provided band gap value is the result of sampling through the k-grid used in VASP calculation only. 

    Carrier effective mass is based on the parabolic E-k relationship approximation, and the convergence test with respect to the grid density near the band edges, and the number of point stencil used, is necessary.

conv
----
Tool for converting file formats.

* Converting cif file to VASP compatible structure file (/examples/conv)

    `vw.py conv cif2vasp -i MoO3.cif -o MoO3_alpha.vasp`

* Converting VASP structure file to Quantum Espresso (pwscf) input file (/examples/conv)

    `vw.py conv vasp2qe -i MoO3.vasp -o pw.in`

* Converting VASP structure file to .xyz format or SHAPE program compatible format (/examples/conv)

    `vw.py conv polyxyz -i MoO3.vasp -c Mo -v O`

    One can toggle on the [SHAPE](http://www.ee.ub.edu) program formatter using -S tag.

boltz
-----

* Postprocessing [BoltzTraP](https://www.sciencedirect.com/science/article/pii/S0010465506001305) generated output

    `vw.py boltz post -s`
