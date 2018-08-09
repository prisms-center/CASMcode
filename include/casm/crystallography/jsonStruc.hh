#ifndef JSONSTRUC_HH
#define JSONSTRUC_HH

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"


namespace CASM {

  /** \ingroup Structure
   *  @{
   */

  template< bool IsConst >
  class SimpleJSonSiteStructure {
  public:
    typedef CASM_TMP::ConstSwitch<IsConst, BasicStructure<Site> > StrucType;
    SimpleJSonSiteStructure(StrucType &_struc, const std::string &_prefix = std::string()) :
      m_struc_ptr(&_struc), m_prefix(_prefix) {};

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json) const {
      _from_json(*m_struc_ptr, json);
    }
  private:
    void _from_json(const BasicStructure<Site> &struc, const jsonParser &json) const {
      std::cerr << "WARNING: Attempting to populate a const Structure/BasicStructure from a JSON object.  This is not allowed.\n";
    }

    void _from_json(BasicStructure<Site> &struc, const jsonParser &json) const;

    StrucType *m_struc_ptr;
    std::string m_prefix;
  };


  template< bool IsConst >
  jsonParser &SimpleJSonSiteStructure<IsConst>::to_json(jsonParser &json) const {
    if(COORD_MODE::IS_FRAC())
      json["coord_mode"] = "direct";
    else
      json["coord_mode"] = "cartesian";

    StrucType &struc(*m_struc_ptr);
    std::map<std::string, std::vector<Site> > site_map;
    for(Index i = 0; i < struc.basis.size(); i++)
      site_map[struc.basis[i].occ_name()].push_back(struc.basis[i]);

    json[m_prefix + "basis"].put_array();
    json["atoms_per_type"].put_array();
    json["atom_type"].put_array();
    json[m_prefix + "lattice"] = struc.lattice();
    auto it = site_map.cbegin(), end_it = site_map.cend();
    for(; it != end_it; ++it) {
      json["atoms_per_type"].push_back(it->second.size());
      json["atom_type"].push_back(it->first);
      auto it2 = it->second.cbegin(), end_it2 = it->second.cend();
      for(; it2 != end_it2; ++it2)
        json[m_prefix + "basis"].push_back(*it2);
    }
    return json;
  }

  template< bool IsConst >
  void SimpleJSonSiteStructure<IsConst>::_from_json(BasicStructure<Site> &struc, const jsonParser &json) const {
    struc.basis.clear();
    struc.reset();
    try {
      std::string tstr;
      CASM::from_json(tstr, json["coord_mode"]);

      COORD_TYPE mode = CART;
      if(tstr == "direct" || tstr == "Direct")
        mode = FRAC;

      Lattice tlat;
      CASM::from_json(tlat, json[m_prefix + "lattice"]);
      struc.set_lattice(tlat, FRAC);
      Index l = 0;

      Eigen::Vector3d tvec;
      for(Index i = 0; i < json["atoms_per_type"].size(); i++) {
        CASM::from_json(tstr, json["atom_type"][i]);
        for(Index j = 0; j < json["atoms_per_type"][i].get<Index>(); j++) {
          CASM::from_json(tvec, json[m_prefix + "basis"][l++]);
          struc.basis.push_back(Site(Coordinate(tvec, struc.lattice(), mode), tstr));
        }
      }
    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Unable to parse Structure/BasicStructure from JSON object.  One or more tags were improperly specified:\n") + ex.what());
    }
    struc.update();

  }

  template< bool IsConst >
  void from_json(const SimpleJSonSiteStructure<IsConst> &jstruc, const jsonParser &json) {
    jstruc.from_json(json);
  }

  template< bool IsConst >
  jsonParser &to_json(const SimpleJSonSiteStructure<IsConst> &jstruc, jsonParser &json) {
    return jstruc.to_json(json);
  }

  inline
  SimpleJSonSiteStructure<true> simple_json(const BasicStructure<Site> &struc, const std::string &prefix) {
    return SimpleJSonSiteStructure<true>(struc, prefix);
  }

  inline
  SimpleJSonSiteStructure<false> simple_json(BasicStructure<Site> &struc, const std::string &prefix) {
    return SimpleJSonSiteStructure<false>(struc, prefix);
  }

  /** @} */
}

#endif
