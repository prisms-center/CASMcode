#include "casm/basis_set/OccupationDoFTraits.hh"

#include <memory>

#include "casm/basis_set/Adapter.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/json_io.hh"  // TODO: remove
#include "casm/casm_io/json/InputParser.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/string_io.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

namespace parse_OccupationDoFSpecs_impl {

/// Reads std::vector<DoF_impl::SublatValues> from JSON
void from_json(std::vector<DoF_impl::SublatValues> &sublat_values,
               const jsonParser &json, std::string key);

void to_json(std::vector<DoF_impl::SublatValues> const &sublat_values,
             jsonParser &json);

/// Reads the "sites_basis_functions" parameter from JSON
void _require_site_basis_functions(
    InputParser<DoF_impl::OccupationDoFSpecs> &parser, const Structure &prim);

}  // namespace parse_OccupationDoFSpecs_impl

namespace DoF_impl {

namespace {
double pretty(double value, double tol) {
  double pretty_value = value;
  if (std::abs(value) < tol) {
    return 0.;
  }
  if (std::abs(std::round(value) - value) < tol) {
    pretty_value = std::round(value);
  }
  return pretty_value;
}
}  // namespace

std::string OccupationDoFSpecs::_name() const {
  return OccupationDoFTraits().name();
}

std::vector<double> chebychev_sublat_prob_vec(Index occupant_dof_size) {
  return std::vector<double>(occupant_dof_size, 1. / double(occupant_dof_size));
}

std::vector<double> occupation_sublat_prob_vec(Index occupant_dof_size) {
  std::vector<double> tprob(occupant_dof_size, 0.);
  tprob[0] = 1.;
  return tprob;
}

/// \brief Returns sublat_values for the given sublattice, ordered to match the
/// order of allowed_occupants
std::vector<double> sublat_values_vec(
    std::vector<SublatValues> const &sublat_values, Index sublat_index,
    std::string key, const std::vector<std::string> &allowed_occupants,
    bool normalize) {
  for (const auto &sublat_value : sublat_values) {
    if (sublat_value.indices.find(sublat_index) != sublat_value.indices.end()) {
      std::vector<double> tprob(allowed_occupants.size(), 0.);
      double tsum{0};
      for (Index ns = 0; ns < allowed_occupants.size(); ns++) {
        auto it = sublat_value.values.at(key).find(allowed_occupants[ns]);
        if (it == sublat_value.values.at(key).end()) {
          throw std::runtime_error("BasisFunctionSpecs error: basis site " +
                                   std::to_string(sublat_index) +
                                   " must have a value specified for species " +
                                   allowed_occupants[ns] + "\n");
        }
        tprob[ns] = it->second;
        tsum += tprob[ns];
      }
      if (normalize) {
        for (Index j = 0; j < tprob.size(); j++) {
          tprob[j] /= tsum;
        }
      }
      return tprob;
    }
  }
  throw std::runtime_error(
      "BasisFunctionSpecs error: values are not specified for basis "
      "site " +
      std::to_string(sublat_index) + "\n");
}

namespace OccupationDoFSpecs_impl {

/// struct for writing validation routines, only used temporarily to construct
/// Validator
struct OccupationDoFSpecsValidator : public Validator {
  OccupationDoFSpecsValidator(OccupationDoFSpecs const &_occ_specs,
                              Structure const &prim);

  OccupationDoFSpecs const &occ_specs;
  std::vector<std::set<std::string>> allowed_molecule_names;

  void check_sublat_molecule_names(
      Index b, DoF_impl::SublatValues const &sublat_value,
      const std::set<std::string> &allowed_molecule_names);
  void check_molecule_names();
  void check_sublat_indices();
};

OccupationDoFSpecsValidator::OccupationDoFSpecsValidator(
    OccupationDoFSpecs const &_occ_specs, Structure const &prim)
    : occ_specs(_occ_specs) {
  // check_dof_keys();
  if (occ_specs.site_basis_function_type ==
      SITE_BASIS_FUNCTION_TYPE::COMPOSITION) {
    /// get allowed_molecule_names and convert to vector of set
    auto tmp{CASM::xtal::allowed_molecule_names(prim)};
    allowed_molecule_names.clear();
    for (const auto &vec : tmp) {
      allowed_molecule_names.push_back(
          std::set<std::string>{vec.begin(), vec.end()});
    }

    check_molecule_names();
    check_sublat_indices();
  }
}

void OccupationDoFSpecsValidator::check_sublat_molecule_names(
    Index b, DoF_impl::SublatValues const &sublat_value,
    std::set<std::string> const &allowed_molecule_names) {
  for (const auto &pair : sublat_value.values) {
    std::string key = pair.first;
    for (const auto &name : allowed_molecule_names) {
      if (sublat_value.values.at(key).find(name) ==
          sublat_value.values.at(key).end()) {
        error.insert("No '" + key + "' is given for '" + name +
                     " on sublattice " + std::to_string(b) +
                     ". Allowed molecules are: " +
                     to_string(allowed_molecule_names) + ".");
      }
    }
    auto begin = allowed_molecule_names.begin();
    auto end = allowed_molecule_names.end();
    for (const auto &pair : sublat_value.values.at(key)) {
      if (std::find(begin, end, pair.first) == end) {
        error.insert(to_string(pair.first) + " is not allowed on sublattice " +
                     std::to_string(b) + ". Allowed molecules are: " +
                     to_string(allowed_molecule_names) + ".");
      }
    }
  }
}

void OccupationDoFSpecsValidator::check_molecule_names() {
  for (const auto &sublat_value : occ_specs.sublat_values) {
    for (const auto &b : sublat_value.indices) {
      if (b >= allowed_molecule_names.size()) {
        error.insert("Sublattice index " + std::to_string(b) +
                     " is out of range.");
      }
      check_sublat_molecule_names(b, sublat_value, allowed_molecule_names[b]);
    }
  }
}

void OccupationDoFSpecsValidator::check_sublat_indices() {
  std::set<Index> found;
  for (const auto &sublat_value : occ_specs.sublat_values) {
    for (const auto &b : sublat_value.indices) {
      if (b >= allowed_molecule_names.size()) {
        error.insert("Sublattice index " + std::to_string(b) +
                     " is out of range.");
      }
      if (found.count(b)) {
        error.insert("Sublattice index " + std::to_string(b) +
                     " is duplicated.");
      }
      found.insert(b);
    }
  }
  for (Index b = 0; b < allowed_molecule_names.size(); ++b) {
    if (!found.count(b) && allowed_molecule_names[b].size() > 1) {
      error.insert("Missing sublattice index " + std::to_string(b) + ".");
    }
  }
}

}  // namespace OccupationDoFSpecs_impl

Validator validate(OccupationDoFSpecs const &occ_specs, const Structure &prim) {
  return OccupationDoFSpecs_impl::OccupationDoFSpecsValidator(occ_specs, prim);
}

// helper for _construct_orthonormal_discrete_functions
std::vector<double> _make_occ_probs(
    Index b_ind, std::vector<std::string> const &site_unique_names,
    OccupationDoFSpecs const &occ_specs) {
  std::vector<double> tprob;

  // options: CHEBYCHEV, OCCUPATION, COMPOSITION
  if (occ_specs.site_basis_function_type ==
      SITE_BASIS_FUNCTION_TYPE::CHEBYCHEV) {
    return chebychev_sublat_prob_vec(site_unique_names.size());
  } else if (occ_specs.site_basis_function_type ==
             SITE_BASIS_FUNCTION_TYPE::OCCUPATION) {
    return occupation_sublat_prob_vec(site_unique_names.size());
  } else if (occ_specs.site_basis_function_type ==
             SITE_BASIS_FUNCTION_TYPE::COMPOSITION) {
    return sublat_values_vec(occ_specs.sublat_values, b_ind, "composition",
                             site_unique_names, true);
  } else {
    throw std::runtime_error(
        "BasisFunctionSpecs error: unknown SITE_BASIS_FUNCTION_TYPE");
  }
}

// helper for OccupationDoFTraits::construct_site_bases
void _construct_orthonormal_discrete_functions(
    BasisSet &site_basis,
    Orbit<PrimPeriodicSymCompare<IntegralCluster>> const &orbit,
    Structure const &prim, DoF_impl::OccupationDoFSpecs const &occ_specs) {
  UnitCellCoord const &uccoord = orbit.prototype()[0];
  Site const &site = uccoord.sublattice_site(prim);
  Index sublattice_index = uccoord.sublattice();
  std::vector<std::vector<std::string>> unique_names =
      prim.structure().unique_names();
  if (unique_names.empty()) {
    unique_names = xtal::allowed_molecule_unique_names(prim);
  }
  // convert vector<xtal::Molecule> to OccupantDoF<xtal::Molecule>
  adapter::Adapter<OccupantDoF<xtal::Molecule>, std::vector<xtal::Molecule>>
      adapter_f;
  OccupantDoF<xtal::Molecule> allowed_occs =
      adapter_f(site.occupant_dof(),
                prim.occupant_symrep_IDs()[sublattice_index], sublattice_index);

  // construct site invariant group
  auto eq_map_range = orbit.equivalence_map(0);
  SymGroup invariant_group{eq_map_range.first, eq_map_range.second};

  if (occ_specs.site_basis_function_type == SITE_BASIS_FUNCTION_TYPE::DIRECT) {
    Index b = sublattice_index;
    Index N = unique_names[b].size();
    // expect N-1 functions (phi0, phi1, etc.) of size N
    std::vector<std::vector<double>> site_basis_functions;

    for (Index i = 0; i < N - 1; ++i) {
      std::string key = "phi" + std::to_string(i);
      site_basis_functions.push_back(sublat_values_vec(
          occ_specs.sublat_values, b, key, unique_names[b], false));
    }

    jsonParser tmp;
    to_json(site_basis_functions, tmp["site_basis_functions"]);
    parse_OccupationDoFSpecs_impl::to_json(occ_specs.sublat_values,
                                           tmp["sublat_values"]);

    // construct site bases for the prototype site
    site_basis.construct_discrete_functions(allowed_occs, site_basis_functions,
                                            sublattice_index, invariant_group);
    return;
  }
  // SITE_BASIS_FUNCTION_TYPE::CHEBYCHEV, ::OCCUPATION, or ::COMPOSITION
  else {
    // construct reference occupation probabilities
    std::vector<double> tprob = _make_occ_probs(
        sublattice_index, unique_names[sublattice_index], occ_specs);

    if (tprob.size() != unique_names[sublattice_index].size()) {
      std::stringstream ss;
      ss << "Error in OccupationDoFTraits::construct_site_bases: occ_probs "
            "error.";
      ss << "  sublattice_index=" << sublattice_index << "  ";
      ss << "  occ_probs=[";
      for (Index i = 0; i < tprob.size(); ++i) {
        if (i != 0) ss << ",";
        ss << tprob[i];
      }
      ss << "]";
      throw std::runtime_error(ss.str());
    }

    // construct site bases for the prototype site
    site_basis.construct_orthonormal_discrete_functions(
        allowed_occs, tprob, sublattice_index, invariant_group);
    return;
  }
}

/// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its
/// site
std::vector<BasisSet> OccupationDoFTraits::construct_site_bases(
    Structure const &_prim,
    std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> &_asym_unit,
    BasisFunctionSpecs const &_basis_function_specs) const {
  // site bases, to be populated
  std::vector<BasisSet> site_bases{_prim.basis().size()};

  // specifications, for how to construct site basis functions
  auto const &occ_specs =
      get<OccupationDoFSpecs>(this->name(), _basis_function_specs);

  // loop over orbits in the asymmetric unit
  for (auto const &orbit : _asym_unit) {
    UnitCellCoord const &prototype_uccoord = orbit.prototype()[0];
    Site const &prototype_site = prototype_uccoord.sublattice_site(_prim);
    Index prototype_sublattice_index = prototype_uccoord.sublattice();

    if (prototype_site.occupant_dof().size() < 2) {
      continue;
    }

    // construct site basis for prototype site in orbit
    auto &prototype_site_basis = site_bases[prototype_sublattice_index];
    _construct_orthonormal_discrete_functions(prototype_site_basis, orbit,
                                              _prim, occ_specs);

    // construct site bases for all equivalent sites by:
    // - applying symmetry to the site basis already constructed for the
    // prototype site
    // - setting necessary indices & ids
    for (Index ne = 1; ne < orbit.size(); ne++) {
      Index other_sublattice_index = orbit[ne][0].sublattice();
      auto &other_site_basis = site_bases[other_sublattice_index];
      other_site_basis = prototype_site_basis;
      other_site_basis.apply_sym(orbit.equivalence_map()[ne][0]);
      other_site_basis.accept(OccFuncBasisIndexer{(int)other_sublattice_index});
      other_site_basis.set_dof_IDs({other_sublattice_index});
    }
  }
  return site_bases;
}

//************************************************************
std::string OccupationDoFTraits::clexulator_point_prepare_string(
    Structure const &_prim,
    std::map<UnitCellCoord, std::set<UnitCellCoord>> const &_nhood,
    PrimNeighborList &_nlist, std::vector<BasisSet> const &site_bases,
    std::string const &indent) const {
  std::stringstream ss;
  ss << indent << "switch(nlist_ind) {\n";
  for (auto const &nbor : _nhood) {
    ss << indent << "case " << _nlist.neighbor_index(nbor.first) << ":\n";
    std::stringstream ssfunc;
    // Put neighborhood in a sensible order:
    std::map<Index, std::set<Index>> sublat_nhood;
    for (auto const &ucc : nbor.second) {
      sublat_nhood[ucc.sublattice()].insert(_nlist.neighbor_index(ucc));
    }

    for (auto const &sublat : sublat_nhood) {
      Index b = sublat.first;
      for (Index n : sublat.second) {
        for (Index f = 0; f < site_bases[b].size(); f++) {
          ssfunc << indent << "    ParamPack::Val<Scalar>::set(m_params, m_"
                 << site_basis_name() << "_param_key, " << f << ", " << n
                 << ", eval_occ_func_" << b << "_" << f << "(" << n << "));\n";
        }
      }
    }
    if (ssfunc.str().size()) {
      ss << indent << "  if(m_params.eval_mode(m_" << site_basis_name()
         << "_param_key) != ParamPack::READ) {\n"
         << ssfunc.str() << indent << "  }\n";
    }
    ss << indent << "  break;\n";
  }
  ss << indent << "}\n";
  return ss.str();
}

//************************************************************

std::string OccupationDoFTraits::clexulator_global_prepare_string(
    Structure const &_prim,
    std::map<UnitCellCoord, std::set<UnitCellCoord>> const &_nhood,
    PrimNeighborList &_nlist, std::vector<BasisSet> const &site_bases,
    std::string const &indent) const {
  std::stringstream ss, ssfunc;
  std::map<Index, std::set<Index>> tot_nhood;
  for (auto const &nbor : _nhood)
    for (auto const &ucc : nbor.second)
      tot_nhood[ucc.sublattice()].insert(_nlist.neighbor_index(ucc));

  for (auto const &nbor : tot_nhood) {
    Index b = nbor.first;
    for (Index n : nbor.second) {
      for (Index f = 0; f < site_bases[b].size(); f++) {
        ssfunc << indent << "  ParamPack::Val<Scalar>::set(m_params, m_"
               << site_basis_name() << "_param_key, " << f << ", " << n
               << ", eval_occ_func_" << b << "_" << f << "(" << n << "));\n";
      }
    }
  }
  if (ssfunc.str().size()) {
    ss << indent << "if(m_params.eval_mode(m_" << site_basis_name()
       << "_param_key) != ParamPack::READ) {\n"
       << ssfunc.str() << indent << "}\n";
  }

  return ss.str();
}
//************************************************************

std::string OccupationDoFTraits::clexulator_member_declarations_string(
    Structure const &_prim, std::vector<BasisSet> const &_site_bases,
    std::string const &indent) const {
  std::stringstream stream;
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> asym_unit;
  std::ostream nullstream(0);

  // TODO: This is a temporary fix to avoid changing the interface of
  // OccupationDoFTraits. What can happen is making these functions take a
  // shared_ptr, but what *really* needs to happen is fixing all the classes and
  // functions (such as those in SymCompare) that are requesting Structures when
  // all they really need is SymRepIDs. I can't fix all that right now though,
  // so check this out: it's a shared pointer to an existing Structure, that has
  // a custom destructor that does nothing.
  auto _prim_ptr =
      std::shared_ptr<const Structure>(&_prim, [](const Structure *) {});
  make_prim_periodic_asymmetric_unit(_prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  for (Index no = 0; no < asym_unit.size(); no++) {
    Index nb = asym_unit[no][0][0].sublattice();
    if (_site_bases[nb].size() == 0) continue;
    stream
        << indent
        << "// Occupation Function tables for basis sites in asymmetric unit "
        << no << ":\n";
    for (Index ne = 0; ne < asym_unit[no].size(); ne++) {
      nb = asym_unit[no][ne][0].sublattice();
      stream << indent << "//   - basis site " << nb << ":\n";
      for (Index f = 0; f < _site_bases[nb].size(); f++) {
        stream << indent << "double "
               << "m_occ_func_" << nb << '_' << f << '['
               << _prim.basis()[nb].occupant_dof().size() << "];\n";
      }
      stream << '\n';
    }
  }
  return stream.str();
}

//************************************************************

std::string OccupationDoFTraits::clexulator_private_method_declarations_string(
    Structure const &_prim, std::vector<BasisSet> const &_site_bases,
    const std::string &indent) const {
  std::stringstream stream;
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> asym_unit;
  std::ostream nullstream(0);

  // TODO: This is a temporary fix to avoid changing the interface of
  // OccupationDoFTraits. What can happen is making these functions take a
  // shared_ptr, but what *really* needs to happen is fixing all the classes and
  // functions (such as those in SymCompare) that are requesting Structures when
  // all they really need is SymRepIDs. I can't fix all that right now though,
  // so check this out: it's a shared pointer to an existing Structure, that has
  // a custom destructor that does nothing.
  auto _prim_ptr =
      std::shared_ptr<const Structure>(&_prim, [](const Structure *) {});
  make_prim_periodic_asymmetric_unit(_prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  for (Index no = 0; no < asym_unit.size(); no++) {
    Index nb = asym_unit[no][0][0].sublattice();
    if (_site_bases[nb].size() == 0) continue;

    for (Index ne = 0; ne < asym_unit[no].size(); ne++) {
      nb = asym_unit[no][ne][0].sublattice();

      stream
          << indent
          << "// Occupation Function evaluators and accessors for basis site "
          << nb << ":\n";
      for (Index f = 0; f < _site_bases[nb].size(); f++) {
        stream << indent << "double const &eval_occ_func_" << nb << '_' << f
               << "(const int &nlist_ind) const {\n"
               << indent << "  return "
               << "m_occ_func_" << nb << '_' << f << "[_occ(nlist_ind)];\n"
               << indent << "}\n\n"
               <<

            indent << "double const &occ_func_" << nb << '_' << f
               << "(const int &nlist_ind) const {\n"
               << indent << "  return "
               << "m_params.read(m_" << site_basis_name() << "_param_key, " << f
               << ", nlist_ind);\n"
               << indent << "}\n";
      }
      stream << '\n';
    }
  }
  return stream.str();
}

//************************************************************

std::vector<DoFType::ParamAllocation>
OccupationDoFTraits::param_pack_allocation(
    Structure const &_prim, std::vector<BasisSet> const &_bases) const {
  std::vector<DoFType::ParamAllocation> result;
  Index NB = 0;
  for (BasisSet const &basis : _bases) {
    NB = max(basis.size(), NB);
  }
  if (NB)
    result.push_back(
        DoFType::ParamAllocation(site_basis_name(), NB, Index(-1), true));

  return result;
}

//************************************************************

std::string OccupationDoFTraits::clexulator_constructor_string(
    Structure const &_prim, std::vector<BasisSet> const &_site_bases,
    const std::string &indent) const {
  std::stringstream stream;
  stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
  stream.precision(10);

  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> asym_unit;
  std::ostream nullstream(0);

  // TODO: This is a temporary fix to avoid changing the interface of
  // OccupationDoFTraits. What can happen is making these functions take a
  // shared_ptr, but what *really* needs to happen is fixing all the classes and
  // functions (such as those in SymCompare) that are requesting Structures when
  // all they really need is SymRepIDs. I can't fix all that right now though,
  // so check this out: it's a shared pointer to an existing Structure, that has
  // a custom destructor that does nothing.
  auto _prim_ptr =
      std::shared_ptr<const Structure>(&_prim, [](const Structure *) {});
  make_prim_periodic_asymmetric_unit(_prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  for (const auto &asym : asym_unit) {
    for (const auto &equiv : asym) {
      Index nb = equiv[0].sublattice();
      for (Index f = 0; f < _site_bases[nb].size(); f++) {
        for (Index s = 0; s < _prim.basis()[nb].occupant_dof().size(); s++) {
          OccFuncEvaluator t_eval(s);
          _site_bases[nb][f]->accept(t_eval);

          if (s == 0) stream << indent;
          stream << "m_occ_func_" << nb << '_' << f << '[' << s
                 << "] = " << t_eval.value();
          if (s + 1 == _prim.basis()[nb].occupant_dof().size())
            stream << ";\n\n";
          else
            stream << ", ";
        }
      }
    }
  }
  return stream.str();
}

std::string OccupationDoFTraits::site_basis_description(BasisSet site_bset,
                                                        Site site,
                                                        Index site_ix) const {
  std::stringstream ss;
  if (site_bset.size() == 0) ss << "        [No site basis functions]\n\n";
  /* std::vector<DoF::RemoteHandle> handles(1, site.occupant_dof().handle()); */
  // TODO: What is "s"? Do I need to know about this "s" elsewhere to use the
  // remote handles
  //!!TODO!! is site_ix the right thing to pass?
  std::vector<DoF::RemoteHandle> handles{
      DoF::RemoteHandle(DoFType::occupation().name(), "s", site_ix)};
  int s;
  handles[0] = s;
  site_bset.register_remotes(handles);
  for (Index f = 0; f < site_bset.size(); f++) {
    for (s = 0; s < site.occupant_dof().size(); s++) {
      if (s == 0) ss << "    ";
      ss << "    \\phi_{" << site_ix << "," << f << "}["
         << site.occupant_dof()[s].name()
         << "] = " << pretty(site_bset[f]->remote_eval(), 1e-10);
      if (s + 1 == site.occupant_dof().size())
        ss << "\n";
      else
        ss << ",   ";
    }
  }
  return ss.str();
}

std::vector<std::unique_ptr<FunctionVisitor>>
OccupationDoFTraits::site_function_visitors(
    std::string const &nlist_specifier) const {
  std::vector<std::unique_ptr<FunctionVisitor>> result;
  result.push_back(std::unique_ptr<FunctionVisitor>(
      new OccFuncLabeler("occ_func_%b_%f(" + nlist_specifier + ")")));
  return result;
}

std::vector<std::unique_ptr<FunctionVisitor>>
OccupationDoFTraits::clust_function_visitors() const {
  std::vector<std::unique_ptr<FunctionVisitor>> result;
  result.push_back(std::unique_ptr<FunctionVisitor>(
      new OccFuncLabeler("occ_func_%b_%f(%n)")));
  return result;
}

void OccupationDoFTraits::parse_dof_specs(
    InputParser<BasisFunctionSpecs> &parser, Structure const &prim) const {
  fs::path dof_specs_path{"dof_specs"};
  fs::path subpath = dof_specs_path / this->name();
  auto subparser = parser.subparse<OccupationDoFSpecs>(subpath, prim);
  if (subparser->value) {
    parser.value->dof_specs.push_back(std::move(subparser->value));
  }
}

void OccupationDoFTraits::dof_specs_to_json(
    const BasisFunctionSpecs &basis_function_specs, jsonParser &json,
    Structure const &prim) const {
  OccupationDoFSpecs const &occupation_dof_specs =
      get<OccupationDoFSpecs>(this->name(), basis_function_specs);
  CASM::to_json(occupation_dof_specs, json["dof_specs"][this->name()]);
}

DoFType::Traits *OccupationDoFTraits::_clone() const {
  return new OccupationDoFTraits(*this);
}

}  // namespace DoF_impl

namespace DoFType {
DoF_impl::OccupationDoFTraits occupation() {
  return DoF_impl::OccupationDoFTraits();
}

std::unique_ptr<DoFSpecs> chebychev_bfuncs() {
  return notstd::make_unique<DoF_impl::OccupationDoFSpecs>(
      DoF_impl::SITE_BASIS_FUNCTION_TYPE::CHEBYCHEV);
}

std::unique_ptr<DoFSpecs> occupation_bfuncs() {
  return notstd::make_unique<DoF_impl::OccupationDoFSpecs>(
      DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION);
}

std::unique_ptr<DoFSpecs> composition_bfuncs(
    std::vector<DoF_impl::SublatValues> sublat_values) {
  return notstd::make_unique<DoF_impl::OccupationDoFSpecs>(
      DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION, sublat_values);
}
}  // namespace DoFType

namespace parse_OccupationDoFSpecs_impl {

/// Reads std::vector<DoF_impl::SublatValues> from JSON
///
/// Example form: (using sublat compositions)
/// \code
/// [
///   {
///     "sublat_indices": [0, 1],
///     "composition": {"A": 0.25, "B": 0.75}
///   },
///   {
///     "sublat_indices": [2, 3],
///     "composition": {"A": 0.75, "B": 0.25}
///   }
/// ]
/// \endcode
///
/// Example form: (using site basis functions, phi0, phi1,)
/// \code
/// [
///   {
///     "sublat_indices": [0, 1],
///     "phi0": {"A": 0.0, "B": 1.0}
///   },
///   {
///     "sublat_indices": [2, 3],
///     "phi0": {"A": 0.0, "B": 1.0, "C": 0.0}
///     "phi1": {"A": 0.0, "B": 0.0, "C": 1.0}
///   }
/// ]
/// \endcode
void from_json(std::vector<DoF_impl::SublatValues> &sublat_values,
               const jsonParser &json, std::string key) {
  sublat_values.clear();
  for (auto it = json.begin(); it != json.end(); ++it) {
    std::map<std::string, std::map<std::string, double>> all_values;
    std::set<Index> sublat_indices;
    it->find("sublat_indices")->get(sublat_indices);
    for (auto inner_it = it->begin(); inner_it != it->end(); ++inner_it) {
      std::string key = inner_it.name();
      if (key == "sublat_indices") {
        continue;
      }
      std::map<std::string, double> single_values;
      inner_it->get(single_values);
      all_values[key] = single_values;
    }
    sublat_values.emplace_back(sublat_indices, all_values);
  }
}

void to_json(std::vector<DoF_impl::SublatValues> const &sublat_values,
             jsonParser &json) {
  json.put_array();
  for (auto const &sublat_value : sublat_values) {
    jsonParser tjson;
    tjson["sublat_indices"] = sublat_value.indices;
    for (auto const &pair : sublat_value.values) {
      tjson[pair.first] = pair.second;
    }
    json.push_back(tjson);
  }
}

/// Reads the "sites_basis_functions" parameter from JSON
///
/// - only validates JSON format to the extent needed to assign values to the
///   OccupationDoFSpecs object, does not validate the values against the prim
void _require_site_basis_functions(
    InputParser<DoF_impl::OccupationDoFSpecs> &parser, const Structure &prim) {
  auto it = parser.self.find("site_basis_functions");
  if (it == parser.self.end()) {
    parser.error.insert("Error: Required 'site_basis_functions' not found.");
  } else {
    if (it->is_string()) {
      DoF_impl::SITE_BASIS_FUNCTION_TYPE site_basis_function_type;
      it->get<DoF_impl::SITE_BASIS_FUNCTION_TYPE>(site_basis_function_type);

      if (site_basis_function_type ==
              DoF_impl::SITE_BASIS_FUNCTION_TYPE::CHEBYCHEV ||
          site_basis_function_type ==
              DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION) {
        parser.value = notstd::make_unique<DoF_impl::OccupationDoFSpecs>(
            site_basis_function_type);
      } else {
        // For SITE_BASIS_FUNCTION_TYPE::COMPOSITION and ::DIRECT
        auto sublat_values_it = parser.self.find("sublat_values");
        if (sublat_values_it == parser.self.end()) {
          parser.error.insert("Error: missing required 'sublat_values'");
        } else {
          std::vector<DoF_impl::SublatValues> sublat_values;
          parse_OccupationDoFSpecs_impl::from_json(sublat_values,
                                                   *sublat_values_it, "values");
          parser.value = notstd::make_unique<DoF_impl::OccupationDoFSpecs>(
              site_basis_function_type, sublat_values);
        }
      }
    } else if (it->is_array()) {  // deprecated
      try {
        std::vector<DoF_impl::SublatValues> sublat_values;
        parse_OccupationDoFSpecs_impl::from_json(sublat_values, *it,
                                                 "composition");
        parser.value = notstd::make_unique<DoF_impl::OccupationDoFSpecs>(
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION, sublat_values);
      } catch (std::exception &e) {
        parser.error.insert(
            "Error: 'site_basis_functions' is an array, but failed to read "
            "sublattice compositions.");
      }
    } else {
      parser.error.insert(
          "Error: 'site_basis_functions' is not a recognized value.");
    }
  }
}

}  // namespace parse_OccupationDoFSpecs_impl

// -- OccupationDoFSpecs IO --

const std::string traits<DoF_impl::SITE_BASIS_FUNCTION_TYPE>::name =
    "site_basis_function_type";

const std::multimap<DoF_impl::SITE_BASIS_FUNCTION_TYPE,
                    std::vector<std::string>>
    traits<DoF_impl::SITE_BASIS_FUNCTION_TYPE>::strval = {
        {DoF_impl::SITE_BASIS_FUNCTION_TYPE::CHEBYCHEV,
         {"CHEBYCHEV", "Chebychev", "chebychev"}},
        {DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION,
         {"OCCUPATION", "Occupation", "occupation"}},
        {DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION,
         {"COMPOSITION", "Composition", "composition"}},
        {DoF_impl::SITE_BASIS_FUNCTION_TYPE::DIRECT,
         {"DIRECT", "Direct", "direct"}}};

ENUM_IO_DEF(DoF_impl::SITE_BASIS_FUNCTION_TYPE)
ENUM_JSON_IO_DEF(DoF_impl::SITE_BASIS_FUNCTION_TYPE)

void parse(InputParser<DoF_impl::OccupationDoFSpecs> &parser,
           const Structure &prim) {
  // read site_basis_functions from JSON (may be string or sublat compositions)
  parse_OccupationDoFSpecs_impl::_require_site_basis_functions(parser, prim);

  // validate the sublat compositions w.r.t. the prim
  parser.insert(DoF_impl::validate(*parser.value, prim));
}

void to_json(DoF_impl::OccupationDoFSpecs const &occupation_dof_specs,
             jsonParser &json) {
  json.put_obj();
  json["site_basis_functions"] = occupation_dof_specs.site_basis_function_type;
  if (occupation_dof_specs.site_basis_function_type ==
          DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION ||
      occupation_dof_specs.site_basis_function_type ==
          DoF_impl::SITE_BASIS_FUNCTION_TYPE::DIRECT) {
    parse_OccupationDoFSpecs_impl::to_json(occupation_dof_specs.sublat_values,
                                           json["sublat_values"]);
  }
}

}  // namespace CASM
