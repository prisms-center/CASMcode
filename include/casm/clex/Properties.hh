#ifndef PROPERTIES_HH
#define PROPERTIES_HH

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  class DeltaProperties;

  //***************************************************************
  /*                      PROPERTIES CLASS

    The Properties class works with the DeltaProperties class using
    'operator+' and 'operator-' as follows:

      Properties calc;
      Properties ref;
      DeltaProperties delta;

      delta = calc - ref;
      calc = ref + delta;


    The Properties class inherits from jsonParser. You can
    read/write/put/get from it just as you would from jsonParser.

    DATA MEMBERS:

    Currently, Properties does not strictly enforce which data
    members it contains.  Please use the following because it is
    expected elsewhere.

    double "energy":
      The energy in eV/(structure) from DFT calculations without
      allowing the atom coordinates to relax

    double "relaxed_energy":
      The energy in eV/(structure) after allowing DFT to relax
      the atom coordinates

    // The following are not implemented yet:

    Array< Vector3<double> > "forces":
      The forces on the atoms in structure, ordered in the same way
      that structure and config_index_to_bijk are arranged. The
      forces are in eV/Angstrom.

    BasicStructure<Coordinate> "structure":
      Contains the structure, ordered in the exact same way that
      config_index_to_bijk is ordered. The Coordinates are in
      Angstroms.

    \ingroup Configuration
   */
  //***************************************************************

  class Properties : public jsonParser {

  public:

    Properties() : jsonParser() {};

    Properties(const jsonParser &json) : jsonParser(json) {};

    Properties(std::istream &stream) : jsonParser(stream) {};

    Properties(const std::string &file_name) : jsonParser(file_name) {};

    Properties(const fs::path &mypath) : jsonParser(mypath) {};

    Properties &operator+=(const DeltaProperties &delta);

    Properties operator+(const DeltaProperties &delta) const;

    DeltaProperties operator-(const Properties &ref) const;

    Properties &operator=(const jsonParser &json) {
      jsonParser::operator=(json);
      return *this;
    }
  };


  //***************************************************************
  /*                      DeltaProperties Class

    The DeltaProperties class works with the Properties class using
    'operator+' and 'operator-' as follows:

      Properties calc;
      Properties ref;
      DeltaProperties delta;

      delta = calc - ref;
      calc = ref + delta;


    The DeltaProperties class inherits from jsonParser. You can
    read/write/put/get from it just as you would from jsonParser.

    DATA MEMBERS:

    Currently, DeltaProperties does not strictly enforce which data
    members it contains.  Please use the following because it is
    expected elsewhere.

    double "energy":
      The energy difference in eV/(structure)

    double "relaxed_energy":
      The relaxed energy difference in eV/(structure)

    // The following are not implemented yet:

    Array< Vector3<double> > "forces":
      The difference in forces on the atoms in structure, ordered
      in the same way that structure and config_index_to_bijk are
      arranged. The forces are in eV/Angstrom.

    Array< Vector3< double > > "displacement":
     The displacements of the atoms in fractional coordinates
     arranged in the same order as relaxed_structure,
     config_index_to_bijk and forces. The displacements are
     calculated such that they are applied after first applying the
     latticeStrain to the reference structure. The fractional
     coordinates are thus in the lattice of the strained structure.

    LatticeStrain "strain":
      The homogenous strain that needs to be applied to the
      reference structure.
   */
  //***************************************************************

  class DeltaProperties : public jsonParser {

  public:

    DeltaProperties() : jsonParser() {};

    DeltaProperties(const jsonParser &json) : jsonParser(json) {};

    DeltaProperties(std::istream &stream) : jsonParser(stream) {};

    DeltaProperties(const std::string &file_name) : jsonParser(file_name) {};

    DeltaProperties(const fs::path &mypath) : jsonParser(mypath) {};

    DeltaProperties(const Properties &calc, const Properties &ref);

    DeltaProperties &operator=(const jsonParser &json) {
      jsonParser::operator=(json);
      return *this;
    }
  };

}

#endif
