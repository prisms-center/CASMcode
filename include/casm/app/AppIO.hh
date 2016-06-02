#ifndef CASM_AppIO
#define CASM_AppIO

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/clex.hh"
#include "casm/casm_io/SafeOfstream.hh"

namespace CASM {

  // --- These functions are for casm I/O -----------

  class ClexBasis;

  // --------- PrimIO Declarations --------------------------------------------------

  BasicStructure<Site> read_prim(fs::path filename);

  BasicStructure<Site> read_prim(const jsonParser &json);

  /// \brief Write prim.json to file
  void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode);

  /// \brief Write prim.json as JSON
  void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode);



  // --------- SymmetryIO Declarations --------------------------------------------------

  void write_symop(const SymGroup &grp, Index i, jsonParser &j);

  void write_symgroup(const SymGroup &grp, jsonParser &json);


  // --------- ChemicalReference IO Declarations --------------------------------------------------

  ChemicalReference read_chemical_reference(fs::path filename, const Structure &prim, double tol);

  ChemicalReference read_chemical_reference(const jsonParser &json, const Structure &prim, double tol);

  void write_chemical_reference(const ChemicalReference &chem_ref, fs::path filename);

  void write_chemical_reference(const ChemicalReference &chem_ref, jsonParser &json);

  void write_chemical_reference(const ChemicalReference &chem_ref, jsonParser &json);


  // --------- CompositionAxes Declarations --------------------------------------------------

  struct CompositionAxes {

    CompositionAxes() {}

    /// \brief Read CompositionAxes from file
    CompositionAxes(fs::path _filename);

    /// \brief Read CompositionAxes from JSON
    CompositionAxes(const jsonParser &json);

    /// \brief Read CompositionAxes from file
    void read(fs::path _filename);

    /// \brief Read CompositionAxes from JSON
    void read(const jsonParser &json, fs::path filename = fs::path());

    /// \brief Write CompositionAxes to file
    void write(fs::path _filename);

    /// \brief Write CompositionAxes to JSON
    void write(jsonParser &json) const;

    /// \brief Set this->curr using key
    void select(std::string key);


    std::map<std::string, CompositionConverter> standard;
    std::map<std::string, CompositionConverter> custom;
    bool has_current_axes = false;
    std::string curr_key;
    CompositionConverter curr;
    fs::path filename;

    int err_code = 0;
    std::string err_message;
  };


  /// \brief Read standard axes from JSON, and output to std::map<std::string, CompositionConverter>
  template<typename OutputIterator>
  OutputIterator read_composition_axes(OutputIterator result, const jsonParser &json);


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

  struct SitesPrinter {

    int indent_space;
    char delim;
    COORD_TYPE mode;


    SitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      indent_space(_indent_space),
      delim(_delim),
      mode(_mode) {}

    std::string indent() const {
      return std::string(indent_space, ' ');
    }

    void coord_mode(std::ostream &out) {
      out << "COORD_MODE = " << mode << std::endl << std::endl;
    }

    void print_sites(const IntegralCluster &clust, std::ostream &out) {
      for(const auto &coord : clust) {
        out << indent() << indent() << indent();
        out.setf(std::ios::showpoint, std::ios_base::fixed);
        out.precision(5);
        out.width(9);
        coord.site().print(out);
        if(delim)
          out << delim;
        out << std::flush;
      }
    }
  };

  struct ProtoSitesPrinter : public SitesPrinter {

    ProtoSitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      SitesPrinter(_indent_space, _delim, _mode) {}


    void operator()(const Orbit<IntegralCluster> &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      out << indent() << indent() << "Prototype" << " of " << orbit.size()
          << " Equivalent Clusters in Orbit " << orbit_index << std::endl;
      print_sites(orbit.prototype(), out);
    }
  };

  struct FullSitesPrinter : public SitesPrinter {

    FullSitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      SitesPrinter(_indent_space, _delim, _mode) {}


    void operator()(const Orbit<IntegralCluster> &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      for(Index equiv_index = 0; equiv_index != orbit.size(); ++equiv_index) {
        out << indent() << indent() << equiv_index << " of " << orbit.size()
            << " Equivalent Clusters in Orbit " << orbit_index << std::endl;
        print_sites(orbit.prototype(), out);
      }
    }
  };

  struct ProtoFuncsPrinter : public SitesPrinter {

    const ClexBasis &clex_basis;

    ProtoFuncsPrinter(const ClexBasis &_clex_basis, int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      SitesPrinter(_indent_space, _delim, _mode),
      clex_basis(_clex_basis) {}


    void operator()(const Orbit<IntegralCluster> &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      out << indent() << indent() << "Prototype" << " of " << orbit.size()
          << " Equivalent Clusters in Orbit " << orbit_index << std::endl;
      print_sites(orbit.prototype(), out);

      throw std::runtime_error("Error printing basis functions: ProtoFuncsPrinter not implemented");
      //print_clust_basis(out, nf, 8, '\n');
      //nf += prototype(i, j).clust_basis.size();
      //out << "\n\n" << std::flush;
    }
  };

  /// \brief Print IntegralCluster orbits
  ///
  /// \param begin,end Range of Orbit<IntegralCluster>
  /// \param out output stream
  /// \param mode Coordinate output mode
  /// \param printer A functor to control printing for the orbit
  ///
  /// Printer is expected to have:
  /// - \code std::string Printer::indent(); \endcode
  /// - \code void Printer::coord_mode(std::ostream& out); \endcode
  /// - \code void Printer::operator()(const Orbit<IntegralCluster>& orbit, std::ostream& out, Index orbit_index, Index Norbits); \endcode
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

      Index b = equiv[0].sublat();
      const auto &site_occ = equiv[0].site().site_occupant();

      auto bset_it = bset_orb_it->begin();
      for(const auto &equiv : *it) {
        const auto &bset = *bset_it++;
        out << "  Basis site " << b << ":\n"
            << "  ";
        equiv[0].site().print(out);
        out << "\n";
        if(bset.size() == 0) {
          out << "        [No site basis functions]\n\n";
        }
        for(Index f = 0; f < bset.size(); ++f) {
          for(Index s = 0; f < site_occ.size(); ++s) {
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

  template<typename ClusterOrbitsIterator>
  jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, jsonParser &bspecs, jsonParser &json) {
    json = jsonParser::object();
    json["orbits"] = jsonParser::array(std::distance(begin, end));
    auto j_it = json["orbits"].begin();
    for(auto it = begin; it != end; ++it) {
      const auto &proto = it->prototype();
      auto &j = (*j_it++)["prototype"];
      j["min_length"] = proto.min_length();
      j["max_length"] = proto.max_length();
      j["sites"] = jsonParser::array();
      for(auto coord_it = proto.begin(); coord_it != proto.end(); ++coord_it) {
        j["sites"].push_back(*coord_it);
      }
    }
    json["bspecs"] = bspecs;
    write_prim(begin->prototype().prim(), json["prim"], FRAC);

    return json;
  }


  // ---------- basis.json IO ------------------------------------------------------------------


  /// \brief Write summary of cluster expansion basis
  template<typename ClusterOrbitsIterator>
  void write_basis(ClusterOrbitsIterator begin, ClusterOrbitsIterator end, const ClexBasis &clex_basis, jsonParser &json, double tol);
}

#endif
