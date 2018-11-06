#include "casm/app/AppIO_impl.hh"
#include "casm/app/HamiltonianModules.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace CASM {


  // --------- PrimIO Definitions --------------------------------------------------

  BasicStructure<Site> read_prim(fs::path filename, HamiltonianModules const &_modules, double xtal_tol) {

    try {
      jsonParser json(filename);
      return read_prim(json, _modules, xtal_tol);
    }
    catch(...) {
      std::cerr << "Error reading prim from " << filename << std::endl;
      /// re-throw exceptions
      throw;
    }
  }

  /// \brief Read prim.json
  BasicStructure<Site> read_prim(const jsonParser &json, HamiltonianModules const &_modules, double xtal_tol) {

    try {

      // read lattice
      Eigen::Matrix3d latvec_transpose;

      from_json(latvec_transpose, json["lattice_vectors"]);

      Lattice lat(latvec_transpose.transpose(), xtal_tol);

      // create prim using lat
      BasicStructure<Site> prim(lat);

      // read title
      prim.set_title(json["title"].get<std::string>());

      Eigen::Vector3d vec;

      // Global DoFs
      {
        std::map<std::string, DoFSet> _dof_map;
        if(json.contains("dof")) {
          auto it = json["dof"].begin(), end_it = json["dof"].end();
          for(; it != end_it; ++it) {
            if(_dof_map.count(it.name()))
              throw std::runtime_error("Error parsing global field \"dof\" from JSON. DoF type " + it.name() + " cannot be repeated.");

            try {
              _dof_map.emplace(std::make_pair(it.name(), it->get<DoFSet>(DoF::traits(it.name()))));
            }
            catch(std::exception &e) {
              throw std::runtime_error("Error parsing global field \"dof\" from JSON. Failure for DoF type " + it.name() + ": " + e.what());
            }

          }
          prim.set_global_dofs(_dof_map);
        }
      }

      // read basis coordinate mode
      COORD_TYPE mode;
      from_json(mode, json["coordinate_mode"]);

      // read basis sites
      for(jsonParser const &bjson : json["basis"]) {
        Site tsite = bjson.get<Site>(prim.lattice(), mode);

        std::vector<Molecule> t_occ;
        std::string occ_key;
        if(bjson.contains("occupants"))
          occ_key = "occupants";
        else if(bjson.contains("occupant_dof"))
          occ_key = "occupant_dof";

        if(!occ_key.empty()) {
          for(std::string const &occ : bjson[occ_key].get<std::vector<std::string> >()) {
            //std::cout << "CREATING OCCUPANT " << occ << "\n";
            //REPLACE THIS LINE FOR GENERALIZED OCCUPANTS
            t_occ.push_back(Molecule::make_atom(occ));
          }
        }
        //std::cout << "t_occ.size() = " << t_occ.size() << "\n";
        if(t_occ.empty())
          t_occ = {Molecule::make_unknown()};
        tsite.set_allowed_occupants(t_occ);
        prim.push_back(tsite);
      }

      return prim;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }

  }

  /// \brief Write prim.json to file
  void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode) {

    SafeOfstream outfile;
    outfile.open(filename);

    jsonParser json;
    write_prim(prim, json, mode);
    json.print(outfile.ofstream());

    outfile.close();
  }

  /// \brief Write prim.json as JSON
  void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode) {

    json = jsonParser::object();

    json["title"] = prim.title();

    json["lattice_vectors"] = prim.lattice().lat_column_mat().transpose();

    if(mode == COORD_DEFAULT) {
      mode = COORD_MODE::CHECK();
    }

    if(mode == FRAC) {
      json["coordinate_mode"] = "Fractional";
    }
    else if(mode == CART) {
      json["coordinate_mode"] = "Cartesian";
    }

    json["basis"] = jsonParser::array(prim.basis().size());
    for(int i = 0; i < prim.basis().size(); i++) {
      json["basis"][i] = jsonParser::object();
      json["basis"][i]["coordinate"].put_array();
      if(mode == FRAC) {
        json["basis"][i]["coordinate"].push_back(prim.basis()[i].frac(0));
        json["basis"][i]["coordinate"].push_back(prim.basis()[i].frac(1));
        json["basis"][i]["coordinate"].push_back(prim.basis()[i].frac(2));
      }
      else if(mode == CART) {
        json["basis"][i]["coordinate"].push_back(prim.basis()[i].cart(0));
        json["basis"][i]["coordinate"].push_back(prim.basis()[i].cart(1));
        json["basis"][i]["coordinate"].push_back(prim.basis()[i].cart(2));
      }

      json["basis"][i]["occupant_dof"] = jsonParser::array(prim.basis()[i].site_occupant().size());

      for(int j = 0; j < prim.basis()[i].site_occupant().size(); j++) {
        json["basis"][i]["occupant_dof"][j] = prim.basis()[i].site_occupant()[j].name();
      }

    }
  }


  // --------- SymmetryIO Declarations --------------------------------------------------

  void write_symop(const SymGroup &grp, Index i, jsonParser &j) {
    j = jsonParser::object();

    const SymOp &op = grp[i];

    to_json(op.matrix(), j["matrix"]["CART"]);
    to_json(grp.lattice().inv_lat_column_mat()*op.matrix()*grp.lattice().lat_column_mat(), j["matrix"]["FRAC"]);

    to_json_array(op.tau(), j["tau"]["CART"]);
    to_json_array(grp.lattice().inv_lat_column_mat()*op.tau(), j["tau"]["FRAC"]);

    to_json(grp.class_of_op(i), j["conjugacy_class"]);
    to_json(grp.ind_inverse(i), j["inverse"]);

    add_sym_info(grp.info(i), j);
  }

  void write_symgroup(const SymGroup &grp, jsonParser &json) {
    json = jsonParser::object();

    json["symop"] = jsonParser::array(grp.size());
    for(int i = 0; i < grp.size(); i++) {
      write_symop(grp, i, json["symop"][i]);
    }

    json["name"] = grp.get_name();
    json["latex_name"] = grp.get_latex_name();
    json["periodicity"] = grp.periodicity();
    if(grp.periodicity() == PERIODIC) {
      json["possible_space_groups"] = grp.possible_space_groups();
    }
    json["conjugacy_class"] = grp.get_conjugacy_classes();
    json["inverse"] = jsonParser::array(grp.size());
    for(int i = 0; i < grp.size(); i++) {
      json["inverse"][i] = grp.ind_inverse(i);
    }
    json["multiplication_table"] = grp.get_multi_table();
    json["character_table"] = grp.character_table();
  }


  // --------- ChemicalReference IO Definitions --------------------------------------------------

  /// \brief Read chemical reference states from JSON file
  ///
  /// See documentation in related function for expected form of the JSON
  ChemicalReference read_chemical_reference(fs::path filename,
                                            const Structure &prim,
                                            double tol) {
    try {
      jsonParser json(filename);
      return read_chemical_reference(json, prim, tol);
    }
    catch(...) {
      std::cerr << "Error reading chemical reference states from " << filename << std::endl;
      /// re-throw exceptions
      throw;
    }
  }

  /// \brief Read chemical reference states from JSON
  ///
  /// Example expected form:
  /// \code
  /// {
  ///   "chemical_reference" : {
  ///     "global" : ...,
  ///     "supercell": {
  ///       "SCELX": ...,
  ///       "SCELY": ...
  ///     },
  ///     "config": {
  ///       "SCELX/I": ...,
  ///       "SCELY/J": ...
  ///     }
  ///   }
  /// }
  /// \endcode
  ///
  /// See one_chemical_reference_from_json for documentation of the \code {...} \endcode expected form.
  ChemicalReference read_chemical_reference(const jsonParser &json,
                                            const Structure &prim,
                                            double tol) {

    if(json.find("chemical_reference") == json.end()) {
      throw std::runtime_error("Error reading chemical reference states: Expected \"chemical_reference\" entry");
    }

    return jsonConstructor<ChemicalReference>::from_json(json["chemical_reference"], prim, tol);

  }

  void write_chemical_reference(const ChemicalReference &chem_ref, fs::path filename) {
    SafeOfstream outfile;
    outfile.open(filename);

    jsonParser json;
    write_chemical_reference(chem_ref, json);
    json.print(outfile.ofstream());

    outfile.close();
  }

  void write_chemical_reference(const ChemicalReference &chem_ref, jsonParser &json) {
    json.put_obj();
    to_json(chem_ref, json["chemical_reference"]);
  }


  // --------- CompositionAxes Definitions --------------------------------------------------

  /// \brief Read CompositionAxes from file
  ///
  /// - Stores err_code and err_message for non-fatal errors
  /// - Catches and throws any exceptions
  CompositionAxes::CompositionAxes(fs::path _filename) {
    read(_filename);
  }

  /// \brief Read CompositionAxes from file
  ///
  /// - Stores err_code and err_message for non-fatal errors
  /// - Catches and throws any exceptions
  CompositionAxes::CompositionAxes(const jsonParser &json) {
    read(json);
  }

  void CompositionAxes::read(fs::path _filename) {

    try {

      read(jsonParser(_filename), _filename);

    }
    catch(...) {

      std::cerr << "Error reading composition axes from " << _filename << std::endl;
      throw;
    }
  }

  void CompositionAxes::read(const jsonParser &json, fs::path _filename) {

    try {

      *this = CompositionAxes();

      filename = _filename;

      filename = _filename;

      if(json.contains("standard_axes")) {
        read_composition_axes(std::inserter(standard, standard.begin()), json["standard_axes"]);
      }

      if(json.contains("custom_axes")) {
        read_composition_axes(std::inserter(custom, custom.begin()), json["custom_axes"]);
      }

      std::string key;
      has_current_axes = json.get_if(key, "current_axes");

      if(has_current_axes) {
        select(key);
      }

    }
    catch(...) {
      throw;
    }
  }

  /// \brief Write CompositionAxes to file
  void CompositionAxes::write(fs::path _filename) {

    filename = _filename;
    SafeOfstream outfile;
    outfile.open(filename);
    jsonParser json;
    write(json);
    json.print(outfile.ofstream());
    outfile.close();
  }

  /// \brief Write CompositionAxes to JSON
  void CompositionAxes::write(jsonParser &json) const {

    json = jsonParser::object();

    json["standard_axes"][standard.cbegin()->first] = true;

    for(auto it = standard.cbegin(); it != standard.cend(); ++it) {
      json["standard_axes"][it->first] = it->second;
    }

    for(auto it = custom.cbegin(); it != custom.cend(); ++it) {
      json["custom_axes"][it->first] = it->second;
    }

    if(has_current_axes) {
      json["current_axes"] = curr_key;
    }

  }

  /// \brief Set this->curr using key
  void CompositionAxes::select(std::string key) {
    if(standard.find(key) != standard.cend() &&
       custom.find(key) != custom.cend()) {

      std::stringstream ss;
      ss << "Error: The composition axes " << key <<
         " can found in both the standard and custom compostion axes.\n\n" <<
         "Please edit the custom composition axes to remove this ambiguity.";

      err_message = ss.str();

      err_code = 1;

    }
    else if(standard.find(key) == standard.cend() &&
            custom.find(key) == custom.cend()) {

      std::stringstream ss;
      ss << "Warning: The composition axes " << key <<
         " can not be found in the standard or custom compostion axes.\n\n" <<
         "Please use 'casm composition --select' to re-select your composition axes,\n" <<
         "Please use 'casm composition --calc' to re-calc your standard axes,\n" <<
         "or edit the custom composition axes.";

      err_message = ss.str();

      err_code = 2;
    }
    else if(standard.find(key) != standard.cend()) {
      curr = standard[key];
    }
    else if(custom.find(key) != custom.cend()) {
      curr = custom[key];
    }
    curr_key = key;
    has_current_axes = true;

  }

  // ---------- prim_nlist.json IO -------------------------------------------------------------
  /*
    void write_prim_nlist(const Array<UnitCellCoord> &prim_nlist, const fs::path &nlistpath) {

      SafeOfstream outfile;
      outfile.open(nlistpath);
      jsonParser nlist_json;
      nlist_json = prim_nlist;
      nlist_json.print(outfile.ofstream());
      outfile.close();

      return;
    }

    Array<UnitCellCoord> read_prim_nlist(const fs::path &nlistpath, const BasicStructure<Site>& prim) {

      try {
        jsonParser json(nlistpath);
        Array<UnitCellCoord> prim_nlist(json.size(), UnitCellCoord(prim));
        auto nlist_it = prim_nlist.begin();
        for(auto it=json.begin(); it!=json.end(); ++it) {
          from_json(*nlist_it++, *it);
        }
        return prim_nlist;
      }
      catch(...) {
        std::cerr << "Error reading: " << nlistpath;
        throw;
      }
    }
  */

  // ---------- Orbit<IntegralCluster> & ClexBasis IO ------------------------------------------------------------------

  const std::string traits<ORBIT_PRINT_MODE>::name = "orbit_print_mode";

  const std::multimap<ORBIT_PRINT_MODE, std::vector<std::string> > traits<ORBIT_PRINT_MODE>::strval = {
    {ORBIT_PRINT_MODE::PROTO, {"PROTO", "Proto", "proto"} },
    {ORBIT_PRINT_MODE::FULL, {"FULL", "Full", "full"} }
  };

  ENUM_IO_DEF(ORBIT_PRINT_MODE)


  PrinterBase::PrinterBase(int _indent_space, char _delim, COORD_TYPE _mode) :
    indent_space(_indent_space),
    indent_level(0),
    delim(_delim),
    mode(_mode) {}

  std::string PrinterBase::indent() const {
    return std::string(indent_space * indent_level, ' ');
  }

  void PrinterBase::coord_mode(Log &out) const {
    out << out.indent_str() << "COORD_MODE = " << mode << std::endl << std::endl;
  }


  const std::string Printer<IntegralCluster>::element_name = "Clusters";

  void Printer<IntegralCluster>::print(const IntegralCluster &clust, Log &out) const {
    if(!out.print()) {
      return;
    }
    COORD_TYPE _mode = mode;
    if(_mode == COORD_DEFAULT) {
      _mode = COORD_MODE::CHECK();
    }
    COORD_MODE printer_mode(_mode);
    if(_mode != INTEGRAL) {
      // calculate nice widths
      int prec = 7;
      int width = prec;
      Eigen::Vector3d vec;
      out.ostream().precision(prec);
      out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      for(const auto &coord : clust) {
        if(_mode == CART) vec = coord.coordinate().cart();
        else if(_mode == FRAC) vec = coord.coordinate().frac();
        width = print_matrix_width(out, vec.transpose(), width);
      }

      // calculate nice widths
      Eigen::IOFormat format(prec, width + 1);
      for(const auto &coord : clust) {
        out << out.indent_str() << indent();
        coord.site().print(out, format);
        if(delim) out << delim;
        out << std::flush;
      }
    }
    else {
      // calculate nice widths
      int prec = 1;
      int width = prec;
      out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      for(const auto &coord : clust) {
        width = print_matrix_width(out, coord.unitcell().transpose(), width);
      }

      // print
      Eigen::IOFormat format(prec, width);
      for(const auto &coord : clust) {
        out << out.indent_str() << indent();
        out << coord << " ";
        coord.site().site_occupant().print(out);
        out << std::flush;
        if(delim) out << delim;
        out << std::flush;
      }
    }
  }

  ProtoFuncsPrinter::ProtoFuncsPrinter(ClexBasis const &_clex_basis, int _indent_space, char _delim, COORD_TYPE _mode) :
    SitesPrinter(_indent_space, _delim, _mode),
    clex_basis(_clex_basis) {
    for(auto const &dofset : clex_basis.site_bases()) {
      for(BasisSet const &bset : dofset.second) {
        if(dofset.first != "occ"
           && bset.size()
           && bset[0]->type_name() != "Variable") {
          labelers.push_back(SubExpressionLabeler(bset.name(), "\\phi^{(" + dofset.first + ")}_%n_%l"));
        }
      }
    }
  }

  void print_site_basis_funcs(Structure const &prim,
                              ClexBasis const &clex_basis,
                              Log &out,
                              Index indent_space,
                              COORD_TYPE mode) {
    std::string indent(indent_space, ' ');

    std::ostream nullstream(0);
    COORD_MODE printer_mode(mode);
    std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
    make_prim_periodic_asymmetric_unit(prim,
                                       CASM_TMP::ConstantFunctor<bool>(true),
                                       TOL,
                                       std::back_inserter(asym_unit),
                                       nullstream);


    for(auto const &dofset : clex_basis.site_bases()) {
      out << indent << indent << "Site basis functions for DoF \"" << dofset.first << "\":\n";
      for(Index no = 0; no < asym_unit.size(); no++) {
        out << indent << indent << "Asymmetric unit " << no + 1 << ":\n";
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].sublat();
          out << indent << indent << "  Basis site " << b << ":\n"
              << "  ";
          if(printer_mode.check() == INTEGRAL) {
            out << indent << indent << asym_unit[no][ne][0] << ' ';
            asym_unit[no][ne][0].site().site_occupant().print(out);
            out << std::flush;
          }
          else
            asym_unit[no][ne][0].site().print(out);

          out << "\n";
          if(dofset.second[b].size() == 0)
            out << "        [No site basis functions]\n\n";
          if(dofset.first == "occ") {
            for(Index f = 0; f < dofset.second[b].size(); f++) {

              BasisSet tbasis(dofset.second[b]);

              int s;
              std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ", "s", prim.basis()[b].site_occupant().ID()));
              remote[0] = s;
              tbasis.register_remotes(remote);

              for(s = 0; s < prim.basis()[b].site_occupant().size(); s++) {
                if(s == 0)
                  out << "    ";
                out << "    \\phi_" << b << '_' << f << '[' << prim.basis()[b].site_occupant()[s].name() << "] = "
                    << tbasis[f]->remote_eval();
                if(s + 1 == prim.basis()[b].site_occupant().size())
                  out << "\n";
                else
                  out << ",   ";
              }
            }
          }
          else {
            //std::string formula;
            //bool relabel=false;
            for(Index f = 0; f < dofset.second[b].size(); f++) {
              if(dofset.second[b][f]) {
                out << "        ";
                if(dofset.second[b][f]->type_name() != "Variable") {
                  out << "\\phi^{" << (dofset.first) << "}_" << b << '_' << f << " = ";
                }
                out << dofset.second[b][f]->tex_formula() << "\n";
              }
            }
          }
        }
      }
    }
    out << "\n\n";


  }

  void write_site_basis_funcs(Structure const &prim,
                              ClexBasis const &clex_basis,
                              jsonParser &json) {


    //   "site_functions":[
    //     {
    //       "asym_unit": X,
    //       "sublat": 2,
    //       "basis": {
    //         "phi_b_0": {"Va":0.0, "O":1.0},
    //         "phi_b_1": {"Va":0.0, "O":1.0}
    //       }
    //     },
    //     ...
    //   ],

    jsonParser &sitef = json["site_functions"];
    sitef = jsonParser::array(prim.basis().size(), jsonParser::object());

    std::ostream nullstream(0);
    std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
    make_prim_periodic_asymmetric_unit(prim,
                                       CASM_TMP::ConstantFunctor<bool>(true),
                                       TOL,
                                       std::back_inserter(asym_unit),
                                       nullstream);


    for(auto const &dofset : clex_basis.site_bases()) {
      for(Index no = 0; no < asym_unit.size(); no++) {

        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].sublat();
          sitef[b]["sublat"] = b;
          sitef[b]["asym_unit"] = no;

          if(dofset.second[b].size() == 0) {
            sitef[b][dofset.first]["basis"].put_null();
            continue;
          }


          for(Index f = 0; f < dofset.second[b].size(); f++) {
            std::stringstream fname;

            if(dofset.first == "occ") {
              fname << "\\phi_" << b << '_' << f;
              BasisSet tbasis(dofset.second[b]);
              int s;
              std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ", "s", prim.basis()[b].site_occupant().ID()));
              remote[0] = s;
              tbasis.register_remotes(remote);

              for(s = 0; s < prim.basis()[b].site_occupant().size(); s++) {
                sitef[b][dofset.first]["basis"][fname.str()][prim.basis()[b].site_occupant()[s].name()] = tbasis[f]->remote_eval();
              }
            }
            else {
              if(dofset.second[b][f]) {
                if(dofset.second[b][f]->type_name() != "Variable") {
                  fname << "\\phi^{" << (dofset.first) << "}_" << b << '_' << f;
                }
                else
                  fname << "var_" << b << '_' << f;
              }
              sitef[b][dofset.first]["basis"][fname.str()] = dofset.second[b][f]->tex_formula();
            }
          }
        }
      }
    }
  }

  // explicit template instantiations

#define PRINT_CLUST_INST(ITERATOR,INSERTER,PRINTER) \
  template void print_clust<ITERATOR,PRINTER>(ITERATOR begin, ITERATOR end, Log &out, PRINTER printer); \
  template jsonParser &write_clust<ITERATOR>(ITERATOR begin, ITERATOR end, jsonParser &json, PRINTER printer); \
  template jsonParser &write_clust<ITERATOR>(ITERATOR begin, ITERATOR end, jsonParser &json, PRINTER printer, const jsonParser &bspecs);

#define ORBIT_CONTAINER_INST(ITERATOR,INSERTER,ORBIT) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,ProtoSitesPrinter) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,FullSitesPrinter) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,ProtoFuncsPrinter) \
  template INSERTER read_clust<INSERTER, typename ORBIT::SymCompareType>(\
    INSERTER result,\
    const jsonParser &json,\
    const Structure &prim,\
    const SymGroup& generating_grp,\
    const typename ORBIT::SymCompareType &sym_compare,\
    double xtal_tol);

#define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
#define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >

#define _SET_IT(ORBIT) std::set<ORBIT>::iterator
#define _SET_INSERTER(ORBIT) std::insert_iterator<std::set<ORBIT> >

#define ORBIT_VECTOR_INST(ORBIT) ORBIT_CONTAINER_INST(_VECTOR_IT(ORBIT),_VECTOR_INSERTER(ORBIT), ORBIT)
#define ORBIT_SET_INST(ORBIT) ORBIT_CONTAINER_INST(_SET_IT(ORBIT),_SET_INSERTER(ORBIT), ORBIT)

  ORBIT_VECTOR_INST(LocalIntegralClusterOrbit)
  ORBIT_VECTOR_INST(PrimPeriodicIntegralClusterOrbit)
  ORBIT_VECTOR_INST(ScelPeriodicIntegralClusterOrbit)
  ORBIT_VECTOR_INST(WithinScelIntegralClusterOrbit)

  ORBIT_SET_INST(LocalIntegralClusterOrbit)
  ORBIT_SET_INST(PrimPeriodicIntegralClusterOrbit)
  ORBIT_SET_INST(ScelPeriodicIntegralClusterOrbit)
  ORBIT_SET_INST(WithinScelIntegralClusterOrbit)


#define DIFFTRANS_VECTOR_INST(ORBIT) \
  PRINT_CLUST_INST(_VECTOR_IT(ORBIT), _VECTOR_INSERTER(ORBIT), PrototypePrinter<Kinetics::DiffusionTransformation>)

  DIFFTRANS_VECTOR_INST(PrimPeriodicDiffTransOrbit)
}

