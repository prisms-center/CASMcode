## CASM: A Clusters Approach to Statistical Mechanics


CASM [(https://github.com/prisms-center/CASMcode)](https://github.com/prisms-center/CASMcode) is an open source software package designed to perform first-principles statistical mechanical studies of multi-component crystalline solids. CASM interfaces with first-principles electronic structure codes, automates the construction and parameterization of effective Hamiltonians and subsequently builds highly optimized (kinetic) Monte Carlo codes to predict finite-temperature thermodynamic and kinetic properties. CASM uses group theoretic techniques that take full advantage of crystal symmetry in order to rigorously construct effective Hamiltonians for almost arbitrary degrees of freedom in crystalline solids. This includes cluster expansions for configurational disorder in multi-component solids and lattice-dynamical effective Hamiltonians for vibrational degrees of freedom involved in structural phase transitions.

This version of CASM supports:

- Constructing, fitting, and evaluating cluster expansion effective Hamiltonians with:
  - Occupational degrees of freedom. 
- High-throughput calculations using:
  - VASP: [https://www.vasp.at](https://www.vasp.at)  
- Semi-Grand canonical Monte Carlo calculations

CASM is updated frequently with support for new effective Hamiltonians, new interfaces for first-principles electronic structure codes, and new Monte Carlo methods. Collaboration is welcome and new features can be incorporated by forking the repository on GitHub, creating a new feature, and submitting pull requests. If you are interested in developing features that involve a significant time investment we encourage you to first contact the CASM development team at <casm-developers@lists.engr.ucsb.edu>.

CASM is currently beta software with very active development. Our goal is that the ``casm`` program interface, including file input and output formats, is mostly stable and backwards compatiblity will be taken into account as new features are added (though some breaking changes may occur as we learn from experience). The CASM library ``libcasm`` is much less stable and we anticipate significant changes will be incorporated in the near future. 

#### Citing CASM

CASM can be cited using the following four references:

- [ref1]        CASM, v0.2.1 (2017). Available from https://github.com/prisms-center/CASMcode. doi:[include doi here]
  - DOIs are generated after a release is archived so they cannot be included in the README for the current release. The appropriate DOI for a particular release can be obtained from the wiki page: <https://github.com/prisms-center/CASMcode/wiki/DOIs>.
- [ref2] 	J. C. Thomas, A. Van der Ven, *Finite-temperature properties of strongly anharmonic and mechanically unstable crystal phases from first principles*, *Physical Review B*, **88**, 214111 (2013).
- [ref3] 	B. Puchala, A. Van der Ven, *Thermodynamics of the Zr-O system from first-principles calculations*, *Physical Review B*, **88**, 094108 (2013).
- [ref4] 	A. Van der Ven, J. C. Thomas, Q. Xu, J. Bhattacharya “Linking the electronic structure of solids to their thermodynamic and kinetic properties”, *Mathematics and Computers in Simulation*, **80**(7) 1393-1410 (2010).

As an example, CASM can be acknowledged in a publication with:

“We used the CASM code [ref1], which automates the construction and parameterization of effective Hamiltonians and implements these Hamiltonians in Monte Carlo simulations [ref2,ref3,ref4].”

#### Citing Algorithms

CASM utilizes a wide variety of algorithms, many of which were developed by the CASM development team, and some of which have yet to be published. Please cite CASM [ref1] if you implement a particular algorithm from CASM in other software. 

CASM also relies on algorithms and methods that have been published in the literature. The cluster expansion for configurational degrees of freedom was rigorously formalized by Sanchez *et al.* [ref5, ref6]. The anharmonic potential cluster expansion as implemented in CASM was developed by Thomas *et al.* [ref2]. The local cluster expansion for diffusion barriers was introduced by Van der Ven *et al.* [ref7]. 

The algorithms in CASM that enumerate symmetrically distinct configurations rely on algebraic properties of principal ideal domains, which were brought to bear on the problem by Hart and Forcade [ref8]. The fitting of the interaction coefficients of a cluster expansion to first-principles data relies on a minimization of the cross-validation (CV) score, an approach introduced to cluster expansions by van de Walle *et al.* [ref9]. The approach of using a genetic algorithm to pick interaction coefficients that minimize the CV score was introduced by Hart *et al.* [ref10] while the depth first search approach is due to Puchala *et al.* [ref3]. The use of compressive sensing methods to parameterize a cluster expansion was introduced by Nelson *et al.* [ref11]. Convergence criteria for Monte Carlo sampling are due to van de Walle *et al.* [ref12]. Convex hulls are found using Qhull [ref13].


- [ref5]: 	J. Sanchez, F. Ducastelle and D. Gratias, *Phys. A* **128**, 334–350 (1984).
- [ref6]: 	D. deFontaine, in *Solid State. Phys.*, ed. H. Ehrenreich and D. Turnbull, Academic Press, vol. 47, pp. 33–176 (1994).
- [ref7]: 	A. Van der Ven, G. Ceder, M. Asta, and P. D. Tepesch, *Phys. Rev. B* **64**, 184307 (2001).
- [ref8]:       G.L.W. Hart and R.W. Forcade, *Phys. Rev. B* **77**, 224115 (2008).
- [ref9]: 	A. van de Walle and G. Ceder, *J. Phase Equilib.* **23**, 348 (2002).
- [ref10]:      G. L. W. Hart, V. Blum, M. J. Walorski, and A. Zunger, *Nat. Mater.* **4**, 391 (2005).
- [ref11]:      L.J. Nelson, G.L.W. Hart, F. Zhou, and V. Ozoliņš, *Phys. Rev. B* **87**, 035125 (2013).
- [ref12]: 	A. van de Walle, M. Asta, *Modell. Simul. Mater. Sci. Eng.* **10**, 521 (2002).
- [ref13]:  C.B. Barber, D.P. Dobkin, and H.T. Huhdanpaa, "The Quickhull algorithm for convex hulls," ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996, http://www.qhull.org.



#### Developers and Contributors:

CASM is developed by the Van der Ven group, originally at the University of Michigan and currently at the University of California Santa Barbara.

**Lead developers**:  John C. Thomas and Brian Puchala

**Developers**:  John Goiri and Anirudh Natarajan.

**Other contributors**: Min-Hua Chen, Jonathon Bechtel, Max Radin, Elizabeth Decolvenaere, Anna Belak, Liang Tian, and Naga Sri Harsha Gunda

#### Acknowledgements ####

The development of CASM was made possible with support from:

- The U.S. Department of Energy, Office of Basic Energy Sciences, Division of Materials Sciences and Engineering under Award #DE-SC0008637 that funds the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan.

- The National Science Foundation under Awards DMR-1410242, DMR-1105672 and DMR-1436154.


#### Contact:

Contact the developers at <casm-developers@lists.engr.ucsb.edu>.

The CASM development team will periodically send email notifications regarding new releases,
features, and bug fixes to the CASM users notification list. To join the list send an email to <CASM-Users-join@lists.engr.ucsb.edu> or visit <https://lists.engr.ucsb.edu/mailman/listinfo/casm-users> to sign up.

#### License

GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.


## Installation

See INSTALL.md


## Getting Started

The ``casm`` executable includes extensive help documentation describing the various commands and options. Simply executing ``casm`` will display a list of possible commands, and executing ``casm <cmd> -h`` will display help documentation particular to the chosen command.

For a beginner, the best place to start is to follow the suggestions printed when calling ``casm status -n``.  This provides step-by-step instructions for creating a CASM project, generating symmetry information, setting composition axes, enumerating configurations, calculating energies with VASP, setting reference states, and fitting an effective Hamiltonian using the program ``casm-learn``. ``casm-learn`` provides The subcommand ``casm format`` provides information on the directory structure of the CASM project and the format of all the CASM files.

All that is needed to start a new project is a ``prim.json`` file describing the crystal structure of the material being studied. See ``casm format --prim`` for a description and examples. Typically one will create a new project directory containing the ``prim.json`` file and then initialize the casm project. For example:

    # create a new CASM project directory
    mkdir /path/to/my_casm_project
    
    # cd into the CASM project directory
    cd /path/to/my_casm_project
    
    # create the prim.json file
    # ... create /path/to/my_casm_project/prim.json ...
    
    # initialize the new CASM project
    casm init


After initializing a casm project: 

- ``casm`` generates code that is compiled and linked at runtime in order to evaluate effective Hamiltonians in a highly optimized manner. If you installed the CASM header files and libraries in a location that is not in your default search path you must specify where to find them. Often the default compilation options work well, but there are some cases when the c++ compiler, compiler flags, or shared object construction flags might need to be customized. You can inspect the current settings via ``casm settings -l`` and options to change them via ``casm settings --desc``.





