#include "casm/clex/Properties.hh"

namespace CASM {

  // Reference Properties + DeltaProperties -> Calculated/Clex Properties
  Properties &Properties::operator+=(const DeltaProperties &delta) {

    // add energy
    if(contains("energy") && delta.contains("energy"))
      (*this)["energy"] = (*this)["energy"].get<double>() + delta["energy"].get<double>();

    // add relaxed_energy
    if(contains("relaxed_energy") && delta.contains("relaxed_energy"))
      (*this)["relaxed_energy"] = (*this)["relaxed_energy"].get<double>() + delta["relaxed_energy"].get<double>();

    // add forces
    /*
    if(contains("forces") && delta.contains("forces")) {
      Array< Vector3<double> > f, f_delta;
      from_json(f, (*this)["forces"]);
      from_json(f_delta, delta["forces"]);
      for(int i = 0; i < f.size(); i++)
        f[i] += f_delta[i];
      (*this)["forces"] = f;
    }
    */

    // add deformation -- a guess at what this could be ...
    /*
    if(contains("structure") && delta.contains("deformation")) {
      BasicStructure<Coordinate> struc;
      Deformation deformation;
      from_json(struc, (*this)["structure"]);
      from_json(deformation, delta["deformation"]);
      (*this)["structure"] = struc * deformation.get_strain() + deformation.get_displacement();
    }
    */

    return *this;
  }

  // Reference Properties + DeltaProperties -> Calculated/Clex Properties
  Properties Properties::operator+(const DeltaProperties &delta) const {
    Properties calc(*this);
    return calc += delta;
  }

  // Calculated Properties - Reference Properties -> DeltaProperties
  DeltaProperties Properties::operator-(const Properties &ref) const {
    return DeltaProperties(*this, ref);
  }

  // Calculated Properties - Reference Properties -> DeltaProperties
  DeltaProperties::DeltaProperties(const Properties &calc, const Properties &ref) : jsonParser() {

    //std::cout << "begin DeltaProperties::DeltaProperties()" << std::endl;
    //std::cout << "calc:\n" << calc << std::endl;
    //std::cout << "ref:\n" << ref << std::endl;

    // energy
    if(calc.contains("energy") && ref.contains("energy"))
      (*this)["energy"] = calc["energy"].get<double>() - ref["energy"].get<double>();

    // relaxed_energy
    if(calc.contains("relaxed_energy") && ref.contains("relaxed_energy"))
      (*this)["relaxed_energy"] = calc["relaxed_energy"].get<double>() - ref["relaxed_energy"].get<double>();


    // forces
    /*
    if(calc.contains("forces") && ref.contains("forces")) {
      Array< Vector3<double> > f, f_ref;
      from_json(f, calc["forces"]);
      from_json(f_ref, ref["forces"]);
      for(int i = 0; i < f.size(); i++)
        f[i] -= f_ref[i];
      (*this)["forces"] = f;
    }
    */

    // deformation -- a guess at what this could be ...
    /*
    if(calc.contains("structure") && ref.contains("structure")) {
      (*this)["deformation"] = Deformation(from_json< BasicStructure<Coordinate> >(calc["structure"]),
                                           from_json< BasicStructure<Coordinate> >(ref["structure"]));
    }
    */
  }

}

