## CASM: A Clusters Approach to Statistical Mechanics


CASM [(https://github.com/prisms-center/CASMcode)](https://github.com/prisms-center/CASMcode) is an open source software package designed to perform first-principles statistical mechanical studies of multi-component crystalline solids. CASM interfaces with first-principles electronic structure codes, automates the construction and parameterization of effective Hamiltonians and subsequently builds highly optimized (kinetic) Monte Carlo codes to predict finite-temperature thermodynamic and kinetic properties. CASM uses group theoretic techniques that take full advantage of crystal symmetry in order to rigorously construct effective Hamiltonians for almost arbitrary degrees of freedom in crystalline solids. This includes cluster expansions for configurational disorder in multi-component solids and lattice-dynamical effective Hamiltonians for vibrational degrees of freedom involved in structural phase transitions.

This version of CASM supports:

- Constructing, fitting, and evaluating cluster expansion effective Hamiltonians with:
  - Occupational degrees of freedom.
- High-throughput calculations using:
  - [VASP](https://www.vasp.at)  
  - [Quantum Espresso](https://www.quantum-espresso.org/)
  - [SeqQuest](https://dft.sandia.gov/Quest/SeqQ_Home.html)
- Monte Carlo calculations using:
  - Semi-grand canonical ensemble
  - Canonical ensemble

CASM is updated frequently with support for new effective Hamiltonians, new interfaces for first-principles electronic structure codes, and new Monte Carlo methods. Collaboration is welcome and new features can be incorporated by forking the repository on GitHub, creating a new feature, and submitting pull requests. If you are interested in developing features that involve a significant time investment we encourage you to first contact the CASM development team at <casm-developers@lists.engr.ucsb.edu>.

CASM is currently beta software with very active development. Our goal is that the ``casm`` program interface, including file input and output formats, is mostly stable and backwards compatiblity will be taken into account as new features are added (though some breaking changes may occur as we learn from experience). The CASM library ``libcasm`` is much less stable and we anticipate significant changes will be incorporated in the near future.


#### Getting Started

- [Installation and online documentation](https://prisms-center.github.io/CASMcode_docs/)

- [Workshop tutorial slides and demo projects](https://github.com/prisms-center/CASMcode_demo)


#### Developers and Contributors:

CASM is developed by the Van der Ven group, originally at the University of Michigan and currently at the University of California Santa Barbara.

**Lead developers**:  John C. Thomas and Brian Puchala

**Developers**:  John Goiri and Anirudh Natarajan

**Other contributors**: Min-Hua Chen, Jonathon Bechtel, Max Radin, Elizabeth Decolvenaere, Anna Belak, Liang Tian, Naga Sri Harsha Gunda, Julija Vinckeviciute, Sanjeev Kolli

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


#### For Developers

See INSTALL.md





