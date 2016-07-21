#ifndef CLEXBASISWRITER_HH
#define CLEXBASISWRITER_HH

#include "casm/CASM_global_definitions.hh"

namespace CASM {
  class ClexBasis;
  class PrimNeighborList;


  class ClexBasisWriter {

  };


  class OrbitFunctionWriter {


  private:
    std::string m_name;
    std::vector<std::string> m_signature;
    std::vector<std::string> arg_names;
    std::string m_short_desc;
    std::string m_long_desc;

  };

  namespace ClexBasisWriter_impl {


    std::string clexulator_member_definitions(std::string const &class_name,
                                              ClexBasis const &clex,
                                              Index N_corr,
                                              Index N_sublat,
                                              std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_private_method_definitions(std::string const &class_name,
                                                      ClexBasis const &clex,
                                                      Index N_corr,
                                                      Index N_sublat,
                                                      std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_public_method_definitions(std::string const &class_name,
                                                     ClexBasis const &clex,
                                                     Index N_corr,
                                                     Index N_sublat,
                                                     std::string const &indent);
    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_orbit_function_strings(ClexBasis::BSetOrbit const &_bset_orbit,
                                                                           OrbitType const &_clust_orbit,
                                                                           std::function<std::string(Index, Index)> method_namer,
                                                                           PrimNeighborList &_nlist,
                                                                           std::vector<std::unique_ptr<FunctionVisitor> > const &labelers);
    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_flower_function_strings(ClexBasis::BSetOrbit const &_bset_orbit,
                                                                            OrbitType const &_clust_orbit,
                                                                            std::function<std::string(Index, Index)> method_namer,
                                                                            PrimNeighborList &_nlist,
                                                                            std::vector<std::unique_ptr<FunctionVisitor> > const &labelers);

    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string, std::string> clexulator_dflower_function_strings(ClexBasis::BSetOrbit const &_bset_orbit,
        OrbitType const &_clust_orbit,
        std::function<std::string(Index, Index)> method_namer,
        PrimNeighborList &_nlist,
        std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
        FunctionVisitor const &prefactor_labeler);
    //*******************************************************************************************

    std::string clexulator_interface_implementation(std::string const &class_name,
                                                    ClexBasis const &clex,
                                                    PrimNeighborList &_nlist,
                                                    std::string const &indent);

    //*******************************************************************************************

    std::string clexulator_constructor_implementation(std::string const &class_name,
                                                      ClexBasis const &clex,
                                                      PrimNeighborList &_nlist,
                                                      std::vector<std::string> const &orbit_method_names,
                                                      std::vector< std::vector<std::string> > const &flower_method_names,
                                                      std::vector< std::vector<std::string> > const &dflower_method_names,
                                                      std::string const &indent);

    //*******************************************************************************************
    /// \brief Print clexulator
    void print_clexulator(const Structure &prim,
                          std::vector<OrbitType > const &_tree
                          const PrimNeighborList &_nlist,
                          std::string class_name,
                          std::ostream &stream,
                          double xtal_tol);

    //*******************************************************************************************

    // Divide by multiplicity. Same result as evaluating correlations via orbitree.
    template<typename OrbitType>
    std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                        OrbitType const &_clust_orbit,
                                                        PrimNeighborList &_nlist,
                                                        std::vector<std::unique_ptr<FunctionVisitor> > const &labelers) {
      std::string prefix, suffix;
      std::vector<std::string> formulae(_bset_orbit[0].size(), std::string());
      if(_clust_orbit.size() > 1) {
        prefix = "(";
        suffix = ")/" + std::to_string(_clust_orbit.size()) + ".";
      }

      for(Index ne = 0; ne < _bset_orbit.size(); ne++) {
        std::vector<PrimNeighborList::Scalar> nbor_IDs =
          _nlist.neighbor_indices(_clust_orbit[ne].elements().begin(),
                                  _clust_orbit[ne].elements().end());
        _bset_orbit[ne].set_dof_IDs(std::vector<Index>(nbor_IDs.begin(), nbor_IDs.end()));
        for(Index nl = 0; nl < labelers.size(); nl++)
          _bset_orbit[ne].accept(*labelers[nl]);
      }

      for(Index nf = 0; nf < _bset_orbit[0].size(); nf++) {
        for(Index ne = 0; ne < _bset_orbit.size(); ne++) {
          if(!_bset_orbit[ne][nf] || (_bset_orbit[ne][nf]->is_zero()))
            continue;

          if(formulae[nf].empty())
            formulae[nf] += prefix;
          else if((_bset_orbit[ne][nf]->formula())[0] != '-' && (_bset_orbit[ne][nf]->formula())[0] != '+')
            formulae[nf] += " + ";

          formulae[nf] += _bset_orbit[ne][nf]->formula();
        }

        if(!formulae[nf].empty())
          formulae[nf] += suffix;
      }
      return formulae;
    }

    //*******************************************************************************************
    /// nlist_index is the index of the basis site in the neighbor list
    template<typename OrbitType>
    std::map< UnitCell, std::vector< std::string > > flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                                                 std::function<BasisSet(BasisSet const &)> _bset_transform,
                                                                                 OrbitType const &_clust_orbit,
                                                                                 PrimNeighborList &_nlist,
                                                                                 std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
                                                                                 Index sublat_index) {


      typedef std::vector<std::string> Formulae;
      typedef std::map<UnitCell, std::vector<std::string> > TransFormulae;
      TransFormulae result;

      std::string prefix, suffix;
      std::set<UnitCellCoord> trans_set;
      for(IntegralCluster const &equiv : _clust_orbit) {
        for(UnitCellCoord const &site : equiv.elements()) {
          if(site.sublat() == sublat_index)
            trans_set.insert(site);
        }
      }

      std::map<UnitCellCoord, std::set<UnitCellCoord> > unique_trans = ClexBasis_impl::unique_ucc(trans_set.begin(),
                                                                       trans_set.end(),
                                                                       _clust_orbit.sym_compare());

      //normalize by multiplicity (by convention)
      if(_clust_orbit.size() > 1) {
        prefix = "(";
        suffix = ")/" + std::to_string(_clust_orbit.size()) + ".";
      }


      for(std::pair<UnitCellCoord, std::set<UnitCellCoord> > const &trans_orbit : unique_trans) {

        Formulae &formulae(result.emplace(std::make_pair(trans_orbit.first, Formulae(_bset_orbit[0].size()))).first->second);

        // loop over equivalent clusters
        for(Index ne = 0; ne < _clust_orbit.size(); ne++) {

          // loop over cluster sites
          auto it(_clust_orbit[ne].cbegin()), end_it(_clust_orbit[ne].cbegin());
          for(; it != end_it; ++it) {

            // Continue if the cluster site doesn't belong to the target sublattice, or if we are working on a different translation
            if(sublat_index != it -> sublat() || trans_orbit.second.find(*it) == trans_orbit.second.end())
              continue;

            IntegralCluster trans_clust = _clust_orbit[ne] - it->unitcell();

            std::vector<PrimNeighborList::Scalar> nbor_IDs =
              _nlist.neighbor_indices(trans_clust.elements().begin(),
                                      trans_clust.elements().end());
            _bset_orbit[ne].set_dof_IDs(std::vector<Index>(nbor_IDs.begin(), nbor_IDs.end()));

            BasisSet transformed_bset(_bset_transform(_bset_orbit[ne]));
            for(Index nl = 0; nl < labelers.size(); nl++)
              transformed_bset.accept(*labelers[nl]);

            for(Index nf = 0; nf < transformed_bset.size(); nf++) {
              if(!transformed_bset[nf] || (transformed_bset[nf]->is_zero()))
                continue;

              if(formulae[nf].empty())
                formulae[nf] += prefix;
              else if((transformed_bset[nf]->formula())[0] != '-' && (transformed_bset[nf]->formula())[0] != '+')
                formulae[nf] += " + ";

              formulae[nf] += transformed_bset[nf]->formula();

            }
          }
        }
      }

      // append suffix to all formulae
      TransFormulae::iterator trans_it = result.begin(),
                              trans_end = result.end();
      for(; trans_it != trans_end; ++trans_it) {
        Formulae &formulae(trans_it->second);
        for(Index nf = 0; nf < formulae.size(); nf++) {
          if(!formulae[nf].empty())
            formulae[nf] += suffix;
        }
      }

      return result;

    }

    //*******************************************************************************************
    template<typename OrbitType>
    void print_proto_clust_funcs(ClexBasis const &_clex_basis,
                                 std::ostream &out,
                                 BasicStructure<Site> const &_prim,
                                 std::vector<OrbitType > const &_tree) {
      //Prints out all prototype clusters (CLUST file)

      out << "COORD_MODE = " << COORD_MODE::NAME() << std::endl << std::endl;

      auto const &site_bases(_clex_basis.site_bases());

      out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
      out.precision(5);
      auto it = site_bases.begin(), end_it = site_bases.end();
      for(; it != end_it; ++it) {
        out << "Basis site definitions for DoF " << it->first << ".\n";

        for(Index nb = 0; nb < (it->second).size(); ++nb) {
          out << "  Basis site " << nb + 1 << ":\n"
              << "  ";
          _prim.basis[nb].print(out);
          out << "\n";
          out << DoFType::traits(it->first).site_basis_description((it->second)[nb], _prim.basis[nb]);
        }
      }

      out << "\n\n";
      Index nf = 0;
      for(Index i = 0; i < _tree.size(); i++) {
        if(i == 0 || _tree[i].prototype().size() != _tree[i - 1].prototype().size())
          out << "** " << _tree[i].prototype().size() << "-site clusters ** \n" << std::flush;

        out << "      ** Orbit " << i + 1 << " of " << _tree.size() << " **"
            << "  Points: " << _tree[i].prototype().size()
            << "  Mult: " << _tree[i].size()
            << "  MinLength: " << _tree[i].prototype().min_length()
            << "  MaxLength: " << _tree[i].prototype().max_length()
            << '\n' << std::flush;

        out << "            " << "Prototype" << " of " << _tree[i].size() << " Equivalent Clusters in Orbit " << i << '\n' << std::flush;
        print_clust_basis(out,
                          _clex_basis.bset_orbit(i)[0],
                          _tree[i].prototype(),
                          nf, 8, '\n');

        nf += _clex_basis.bset_orbit(i)[0].size();
        out << "\n\n" << std::flush;

        if(_tree[i].prototype().size() != 0) out << '\n' << std::flush;
      }
    }

    //*******************************************************************************************


    template<typename UCCIterType, typename IntegralClusterSymCompareType>
    std::map<UnitCellCoord, std::set<UnitCellCoord> > unique_ucc(UCCIterType begin,
                                                                 UCCIterType end,
                                                                 IntegralClusterSymCompareType const &sym_compare) {
      std::map<UnitCellCoord, std::set<UnitCellCoord> > tresult;

      typedef IntegralCluster cluster_type;
      typedef Orbit<cluster_type, IntegralClusterSymCompareType> orbit_type;

      if(begin == end)
        return tresult;
      // store orbits as we find them
      std::set<orbit_type> orbits;

      SymGroup identity_group(begin->prim().factor_group.begin(), begin->prim().factor_group.begin()++);
      orbit_type empty_orbit(cluster_type(begin->prim()), identity_group, sym_compare);

      // by looping over each site in the grid,
      for(; begin != end; ++begin) {

        // create a test cluster from prototype
        cluster_type test(empty_orbit.prototype());

        // add the new site
        test.elements().push_back(*begin);


        // try to find test cluster in already found orbits
        auto it = find_orbit(orbits.begin(), orbits.end(), test);
        if(it != orbits.end()) {
          tresult[it->prototype().element[0]].insert(*begin);
          continue;
        }

        tresult[test.element(0)].insert(*begin);

        // if not yet found, use test to generate a new Orbit
        orbits.insert(orbit_type(test, identity_group, sym_compare));
      }

      std::map<UnitCellCoord, std::set<UnitCellCoord> > result;
      for(auto &_set_pair : tresult)
        result.emplace(std::make_pair(*_set_pair.second.begin(), std::move(_set_pair.second)));
      return result;

    }
  }

}
#endif
