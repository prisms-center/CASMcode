#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ConfigIO.hh"

/// What is being used to test it:

#include "casm/app/ProjectBuilder.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "Common.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigIOTest)

BOOST_AUTO_TEST_CASE(DatumFormatters) {

  using namespace ConfigIO;

  DataFormatterDictionary<Configuration> dict;

  auto check = [&](const BaseDatumFormatter<Configuration> &formatter) {
    auto prev_size = dict.size();
    dict.insert(formatter);
    BOOST_CHECK_MESSAGE(prev_size + 1 == dict.size(), formatter.name());
  };

  // String
  check(configname());
  check(scelname());

  // Boolean
  check(OnHull());
  check(OnClexHull());
  check(is_calculated());
  check(is_primitive());
  check(is_canonical());
  check(selected_in());

  // Integer
  check(scel_size());
  check(multiplicity());

  // Scalar
  //check(Clex());
  check(HullDist());
  check(ClexHullDist());
  check(Novelty());
  check(relaxed_energy());
  check(relaxed_energy_per_species());
  check(reference_energy());
  check(reference_energy_per_species());
  check(formation_energy());
  check(formation_energy_per_species());
  check(rms_force());
  check(basis_deformation());
  check(lattice_deformation());
  check(volume_relaxation());

  // VectorXd
  check(Corr());
  check(CompN());
  check(Comp());
  check(AtomFrac());
  check(RelaxationStrain());
  check(SiteFrac());
  check(StrucScore());

}

BOOST_AUTO_TEST_CASE(Make) {

  std::stringstream ss;
  typedef notstd::cloneable_ptr< BaseDatumFormatter<Configuration> > cloner;
  {
    DataFormatterDictionary<Configuration> dict0;

    using namespace ConfigIO;

    dict0.insert(
      Corr(),
      CompN(),
      format_operator_add<Configuration>()
    );

    DataFormatterDictionary<Configuration> dict;

    dict.insert(
      dict0,
      AtomFrac(),
      SiteFrac()
    );


    ss << "\n---\n";

    ss << "get iterator" << std::endl;
    DataFormatterDictionary<Configuration>::const_iterator it = dict.begin();
    ss << "get ref from transform iterator" << std::endl;
    BaseDatumFormatter<Configuration> &transform_ref = *it;
    (void) transform_ref;
    ss << "get const ref from transform iterator" << std::endl;
    const BaseDatumFormatter<Configuration> &const_transform_ref = *it;
    (void) const_transform_ref;

    ss << "get ref from iterator" << std::endl;
    const BaseDatumFormatter<Configuration> &ref = *(it.base()->second);
    (void) ref;
    ss << "get pair ref from base iterator" << std::endl;
    auto &pair_ref = *(it.base());
    (void) pair_ref;
    ss << "get const pair ref from base iterator" << std::endl;
    const auto &const_pair_ref = *(it.base());
    (void) const_pair_ref;
    ss << "use ref from iterator, name: " << ref.name() << std::endl;

    ss << "\n---\n";
    for(DataFormatterDictionary<Configuration>::const_iterator it = dict.begin(); it != dict.end(); ++it) {
      ss << "name: " << it.base()->second->name() << std::endl;
    }
  }

  ss << "\n-----------------------\n";
  {
    std::vector<cloner> dict;
    dict.emplace_back(ConfigIO::Corr());
    dict.emplace_back(ConfigIO::CompN());
    dict.emplace_back(format_operator_add<Configuration>());
    dict.emplace_back(format_operator_sub<Configuration>());
    dict.emplace_back(format_operator_eq<Configuration>());

    ss << "get iterator" << std::endl;
    auto it = dict.begin();
    ss << "get ref from iterator" << std::endl;
    const BaseDatumFormatter<Configuration> &ref = **it;
    ss << "use ref from iterator, name: " << ref.name() << std::endl;

    for(auto it = dict.begin(); it != dict.end(); ++it) {
      ss << "name: " << (*it)->name() << std::endl;
    }
  }

  ss << "\n-----------------------\n";
  {
    std::map<std::string, cloner > dict;
    dict["corr"] = cloner(ConfigIO::Corr());
    dict["comp_n"] =  cloner(ConfigIO::CompN());
    dict["add"] = cloner(format_operator_add<Configuration>());
    dict["sub"] = cloner(format_operator_sub<Configuration>());
    dict["eq"] = cloner(format_operator_eq<Configuration>());

    ss << "get iterator" << std::endl;
    auto it = dict.begin();
    ss << "get ref from iterator" << std::endl;
    const BaseDatumFormatter<Configuration> &ref = *(it->second);
    ss << "use ref from iterator, name: " << ref.name() << std::endl;

    for(auto it = dict.begin(); it != dict.end(); ++it) {
      ss << "name: " << it->second->name() << std::endl;
    }
  }

  ss << "\n-----------------------\n";
  auto adder = format_operator_add<Configuration>();
  ss << "adder name: " << adder.name() << std::endl;

  BaseDatumFormatter<Configuration> *ptr = &adder;
  ss << "ptr adder name: " << ptr->name() << std::endl;

  notstd::cloneable_ptr< BaseDatumFormatter<Configuration> > cl_ptr(adder);
  ss << "cl_ptr adder name: " << cl_ptr->name() << std::endl;

  auto other = cl_ptr;
  ss << "other adder name: " << other->name() << std::endl;

  notstd::cloneable_ptr< BaseDatumFormatter<Configuration> > another(cl_ptr);
  ss << "another adder name: " << another->name() << std::endl;

  ss << "\n-----------------------\n";
  {

    DataFormatterDictionary<Configuration> dict = make_dictionary<Configuration>();
    for(auto it = dict.begin(); it != dict.end(); ++it) {
      ss << "name: " << it->name() << std::endl;
    }
  }

  //std::cout << ss.str() << std::endl;
}

BOOST_AUTO_TEST_CASE(AllTest) {

  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();
  proj.check_enum();

  Log log = null_log();
  PrimClex primclex(proj.dir, log);

  log << "---- Comp -------------" << std::endl;
  ConfigIO::Comp comp;
  for(auto it = primclex.config_begin(); it != primclex.config_end(); ++it) {
    log << "name: " << it->name() << "  comp: " << comp(*it).transpose() << "  print: ";
    comp.print(*it, log);
    log << std::endl;
  }

  log << "---- BaseValueFormatter Ptr -------------" << std::endl;
  BaseValueFormatter<Eigen::VectorXd, Configuration> *value_ptr = &comp;
  for(auto it = primclex.config_begin(); it != primclex.config_end(); ++it) {
    log << "name: " << it->name() << "  value: " << (*value_ptr)(*it).transpose() << "  print: ";
    value_ptr->print(*it, log);
    log << std::endl;
  }

  log << "---- BaseDatumFormatter Ptr -------------" << std::endl;
  BaseDatumFormatter<Configuration> *datum_ptr = &comp;
  for(auto it = primclex.config_begin(); it != primclex.config_end(); ++it) {
    log << "name: " << it->name() << "  print: ";
    datum_ptr->print(*it, log);
    log << std::endl;
  }

  //std::cout << ss.str() << std::endl;

}

BOOST_AUTO_TEST_SUITE_END()
