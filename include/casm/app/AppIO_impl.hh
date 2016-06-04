#ifndef CASM_AppIO_impl
#define CASM_AppIO_impl

#include "casm/app/AppIO.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
//#include "casm/symmetry/SymInfo.hh"
//#include "casm/symmetry/SymGroup.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ClexBasis.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/clex.hh"
#include "casm/casm_io/SafeOfstream.hh"


namespace CASM {

  // --------- CompositionAxes Definitions --------------------------------------------------

  /// \brief Read standard axes from JSON, and output to std::map<std::string, CompositionConverter>
  ///
  /// - json Format: {"standard_axes": {"key0" : {CompositionConverter}, ... more axes ...},
  ///                 "custom_axes":   {"key0" : {CompositionConverter}, ... more axes ...},
  ///                 "current_axes": "key"}
  ///   - "axes_type" is either "standard_axes" or "custom_axes"
  ///   - "key" is a string indicating which composition axes in the "standard_axes" or "custom_axes" JSON object
  ///
  /// - This read json["standard_axes"] or json["custom_axes"]
  ///
  template<typename OutputIterator>
  OutputIterator read_composition_axes(OutputIterator result, const jsonParser &json) {

    CompositionConverter conv;
    for(auto it = json.cbegin(); it != json.cend(); ++it) {
      from_json(conv, *it);
      *result++ = std::make_pair(it.name(), conv);
    }
    return result;
  }


  // ---------- casm bspecs --orbits -----------------------------------------------------------



  /// \brief Print IntegralCluster orbits
  ///
  /// \param begin,end Range of Orbit<IntegralCluster, SymCompareType>
  /// \param out output stream
  /// \param mode Coordinate output mode
  /// \param printer A functor to control printing for the orbit
  ///
  /// Printer is expected to have:
  /// - \code std::string Printer::indent(); \endcode
  /// - \code void Printer::coord_mode(std::ostream& out); \endcode
  /// - \code void Printer::operator()(const Orbit<IntegralCluster, SymCompareType>& orbit, std::ostream& out, Index orbit_index, Index Norbits); \endcode
  ///
  template<typename ClusterOrbitIterator, typename OrbitPrinter>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    std::ostream &out,
    OrbitPrinter printer) {

    printer.coord_mode(out);

    out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.precision(5);

    Index branch = -1;
    Index orbit_index = 0;
    Index Norbits = std::distance(begin, end);
    std::string indent = printer.indent();

    for(auto it = begin; it != end; ++it) {
      if(it->prototype().size() != branch) {
        branch = it->prototype().size();
        out << "** Branch " << branch << " ** " << std::endl;
      }
      out << printer.indent() << "** " << orbit_index << " of " << Norbits << " Orbits **"
          << "  Points: " << it->prototype().size()
          << "  Mult: " << it->size()
          << "  MinLength: " << it->prototype().min_length()
          << "  MaxLength: " << it->prototype().max_length() << std::endl;
      printer(*it, out, orbit_index, Norbits);
      out << std::endl;
      ++orbit_index;
    }

  }

  /// \brief Print site basis functions, as for 'casm bset --functions'
  ///
  /// \param begin,end Range of Orbit<IntegralCluster, SymCompareType>
  /// \param clex_basis A ClexBasis object generated from the provided orbits
  /// \param out output stream
  /// \param mode Coordinate output mode
  ///
  template<typename ClusterOrbitIterator>
  void print_site_basis_funcs(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    const ClexBasis &clex_basis,
    std::ostream &out,
    COORD_TYPE mode) {

    out << "COORD_MODE = " << mode << std::endl << std::endl;

    out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.precision(5);
    Index asym_unit_count = 1;

    auto bset_orb_it = clex_basis.begin();
    for(auto it = begin; it != end; ++it, ++bset_orb_it) {
      if(it->prototype().size() != 1) {
        continue;
      }

      out << "Asymmetric unit " << asym_unit_count++ << ":\n";


      auto bset_it = bset_orb_it->begin();
      for(const auto &equiv : *it) {

        auto b = equiv[0].sublat();
        const auto &site_occ = equiv[0].site().site_occupant();
        const auto &bset = *bset_it++;

        out << "  Basis site " << b << ":\n"
            << "  ";
        equiv[0].site().print(out);
        out << "\n";
        if(bset.size() == 0) {
          out << "        [No site basis functions]\n\n";
        }
        for(Index f = 0; f < bset.size(); ++f) {
          for(Index s = 0; s < site_occ.size(); ++s) {
            if(s == 0)
              out << "    ";
            out << "    \\phi_" << b << '_' << f << '[' << site_occ[s].name << "] = "
                << bset[f]->eval(Array<Index>(1, site_occ.ID()), Array<Index>(1, s));
            if(s + 1 == site_occ.size())
              out << "\n";
            else
              out << ",   ";
          }
        }
      }
    }
  }


  // ---------- clust.json IO ------------------------------------------------------------------

  /// \brief Read JSON containing Orbit<IntegralCluster, SymCompareType> prototypes
  ///
  /// - Uses 'prim', 'generating_grp', and 'sym_compare' to generate orbits from
  ///   prototypes read from the JSON
  /// - Ignores "prim" and "bspecs" info in the JSON
  ///
  template<typename ClusterOutputIterator, typename SymCompareType>
  ClusterOutputIterator read_clust(
    ClusterOutputIterator result,
    jsonParser &json,
    const BasicStructure<Site> &prim,
    const SymGroup &generating_grp,
    const SymCompareType &sym_compare) {

    for(auto it = json["orbits"].begin(); it != json["orbits"].end(); ++it) {

      // read prototype
      auto &j = (*it)["prototype"]["sites"];
      std::vector<UnitCellCoord> coord(j.size(), UnitCellCoord(prim));
      auto coord_it = coord.begin();
      for(auto j_it = j.begin(); it != j.end(); ++j_it) {
        from_json(*coord_it++, *j_it);
      }

      IntegralCluster proto(prim, coord.begin(), coord.end());
      *result++ = Orbit<IntegralCluster, SymCompareType>(proto, generating_grp, sym_compare);
    }
    return result;
  }

  /// \brief Write JSON containing Orbit<IntegralCluster, SymCompareType> prototypes
  ///
  /// Format:
  /// \code
  /// {
  ///   "orbits": [
  ///     {
  ///       "prototype" : {
  ///         "min_length" : number,
  ///         "max_length" : number,
  ///         "sites" : (JSON array of UnitCellCoord)
  ///       }
  ///     },
  ///     ... for each orbit ...
  ///   ],
  ///   "bspecs" : (JSON object, the contents of bspecs.json)
  ///   "prim" : (JSON object, the contents of prim.json)
  /// }
  /// \endcode
  ///
  template<typename ClusterOrbitIterator>
  jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, jsonParser &bspecs, jsonParser &json) {
    json = jsonParser::object();
    json["orbits"] = jsonParser::array(std::distance(begin, end));
    auto j_it = json["orbits"].begin();
    for(auto it = begin; it != end; ++it) {
      const auto &proto = it->prototype();
      auto &j = (*j_it++)["prototype"];
      j["min_length"] = proto.min_length();
      j["max_length"] = proto.max_length();
      j["sites"].put_array(proto.begin(), proto.end());
    }
    json["bspecs"] = bspecs;
    write_prim(begin->prototype().prim(), json["prim"], FRAC);

    return json;
  }


  // ---------- basis.json IO ------------------------------------------------------------------

  /// \brief Write summary of cluster expansion basis
  ///
  /// Format:
  /// \code
  /// {
  ///   "site_functions":[
  ///     {
  ///       "asym_unit": X,
  ///       "sublat_indices: [2, 3],
  ///       "phi_b_0": {"Va":0.0, "O":1.0},
  ///       "phi_b_1": {"Va":0.0, "O":1.0},
  ///        ...
  ///     },
  ///     ...
  ///   ],
  ///   "cluster_functions":[
  ///     {
  ///       "eci": X.XXXXX,
  ///       "prototype_function": "\phi_b_i(s_j)...",
  ///       "orbit": [branch_index, orbit_index],
  ///       "linear_orbit_index": I,
  ///       "mult": X,
  ///       "prototype": [
  ///         [b, i, j, k],
  ///         ...
  ///       ]
  ///     },
  ///     ...
  ///   ]
  /// }
  /// \endcode
  ///
  /// \brief Write summary of cluster expansion basis
  template<typename ClusterOrbitIterator>
  void write_basis(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    const ClexBasis &clex_basis,
    jsonParser &json,
    double tol) {

    throw std::runtime_error("write_basis needs to be reimplemented");
    /*
      json = jsonParser::object();

      //   "site_functions":[
      //     {
      //       "asym_unit": X,
      //       "sublat": [2, 3],
      //       "basis": {
      //         "phi_b_0": {"Va":0.0, "O":1.0},
      //         "phi_b_1": {"Va":0.0, "O":1.0}
      //       }
      //     },
      //     ...
      //   ],

      jsonParser &sitef = json["site_functions"];
      sitef = jsonParser::array(prim.basis.size(), jsonParser::object());
      for(Index no = 0; no < tree.asym_unit().size(); no++) {
        for(Index ne = 0; ne < tree.asym_unit()[no].size(); ne++) {

          const SiteCluster &equiv = tree.asym_unit()[no][ne];
          const Site &site = equiv[0];

          Index b = site.basis_ind();
          sitef[b]["sublat"] = b;
          sitef[b]["asym_unit"] = no;

          if(equiv.clust_basis.size() == 0) {
            sitef[b]["basis"].put_null();
          }
          else {
            for(Index f = 0; f < equiv.clust_basis.size(); f++) {
              std::stringstream fname;
              fname << "\\phi_" << b << '_' << f;
              for(Index s = 0; s < site.site_occupant().size(); s++) {

                // "\phi_b_f": {"Zr":0.0, ...}
                sitef[b]["basis"][fname.str()][site.site_occupant()[s].name] =
                  equiv.clust_basis[f]->eval(
                    Array<Index>(1, site.site_occupant().ID()),
                    Array<Index>(1, s)
                  );
              }
            }
          }
        }
      }

      //   "cluster_functions":[
      //     {
      //       ("eci": X.XXXXX,) <-- is included after fitting
      //       "prototype_function": "\phi_b_i(s_j)...",
      //       "orbit": [branch_index, cluster_orbit_index, bfunc_index],
      //       "linear_orbit_index": I,
      //       "mult": X,
      //       "prototype": {
      //         "max_length": X.X,
      //         "min_length": X.X,
      //         "sites": [
      //           [b, i, j, k],
      //           ...
      //         ]
      //       }
      //     },
      //     ...
      //   ]
      // }

      jsonParser &orbitf = json["cluster_functions"];
      orbitf = jsonParser::array();
      for(Index i = 0; i < tree.size(); i++) {
        jsonParser tjson;
        for(Index j = 0; j < tree.size(i); j++) { //Loops over all i sized Orbits of clusters

          // cluster specific info
          tjson["orbit"] = std::vector<Index>({i, j, 0});
          tjson["mult"] = tree.orbit(i, j).size();
          to_json(jsonHelper(tree.orbit(i, j)[0], prim), tjson["prototype"]);

          // basis function info
          BasisSet tbasis(tree.orbit(i, j)[0].clust_basis);
          tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));

          for(Index nf = 0; nf < tbasis.size(); ++nf) {
            tjson["orbit"][2] = nf;
            tjson["prototype_function"] = tbasis[nf]->tex_formula();
            orbitf.push_back(tjson);
          }
        }
      }

      for(Index i = 0; i < orbitf.size(); ++i) {
        orbitf[i]["linear_function_index"] = i;
      }
      */
  }
}

#endif
