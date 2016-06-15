#include "casm/basis_set/DoFTraits.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {
  template<typename OrbitIteratorType>
  void ClexBasis::generate(OrbitIteratorType _orbit_begin,
                           OrbitIteratorType _orbit_end,
                           jsonParser const &_bspecs,
                           Index _max_poly_order /*= -1*/) {
    std::vector<DoFKey> dof_keys;
    _bspecs.get_if(dof_keys, "dofs");
    std::vector<DoFKey> global_keys;
    std::vector<DoFKey> local_keys;

    //separate local_args from global_args
    for(DoFKey const &key : dof_keys) {
      if(m_global_bases.find(key) != m_global_bases.end()) {
        global_keys.push_back(key);
      }
      else if(m_site_bases.find(key) != m_site_bases.end()) {
        local_keys.push_back(key);
      }
      else {
        throw std::runtime_error(std::string("Attempting to build Clex basis set, but missing degree of freedom \"") + key + "\n");
      }
    }
    m_bset_tree.resize(std::distance(_orbit_begin, _orbit_end));

    auto bset_it = m_bset_tree.begin();
    for(; _orbit_begin != _orbit_end; ++_orbit_begin, ++bset_it) {
      bset_it->reserve(_orbit_begin->size());
      bset_it->push_back(_construct_prototype_basis(*_orbit_begin,
                                                    local_keys,
                                                    global_keys,
                                                    -1/* polynomial_order */));
      for(Index j = 1; j < _orbit_begin->size(); j++) {
        bset_it->push_back((*(_orbit_begin->equivalence_map(j).first)) * (*bset_it)[0]);
      }
    }
  }


  //*************************************************
  // @param local_args[i][j] is BasisSet for i'th DoFspace at j'th site of cluster
  //                                site  0            site 1            site 2
  // DoFKey 0: displacement       {x0,y0,z0}       {x1, y1, z1}       {x2,y2,z2}
  // DoFKey 1: configuration       {pA0,pB0}            {}            {pA2, pB2}
  //
  // Step 1:  Get the kroenecker product of cluster permutation with DoF symrep
  //
  //  permutation  |    kronecker prod  |    DoF Symrep (e.g., x--y displacement)
  //   [ 0  1 ]              v                     [cos -sin]
  //   [ 1  0 ]             XkX                    [sin  cos]
  //
  //       [x0]      [  0    0   cos -sin ]   [x0]
  //       [y0]      [  0    0   sin  cos ]   [y0]
  //   S * [x1]  =   [ cos -sin   0    0  ]   [x1]
  //       [x2]      [ sin  cos   0    0  ]   [x2]
  //
  // ----------------------------------------------------------------------
  //
  // Step 2: mix-in @param global_args to get all_argsets
  //
  // GLOBAL ARGS
  //
  // strain             { e1, e2, e3, e4, e5, e6}
  // composition        { comp_a, comp_b }
  //
  // arg_subsets =   [ {x0,y0,z0,x1,y1,z1,x2,y2,z2},
  //                   {pA0,pB0,pA2,pB2},
  //                   {e1,e2,e3,e4,e5,e6},
  //                   {comp_a,comp_b}]
  //
  template<typename OrbitType>
  BasisSet ClexBasis::_construct_prototype_basis(OrbitType const &_orbit,
                                                 std::vector<DoFKey> const &local_keys,
                                                 std::vector<DoFKey> const &global_keys,
                                                 Index max_poly_order) const {
    //std::cout<<"In IntegralCluster::generate_clust_basis, the size of this cluster is:"<<size()<<std::endl;
    //std::cout<<"valid_index evaluates to:"<<valid_index(max_poly_order)<<std::endl;

    // Default polynomial order is cluster size
    if(!valid_index(max_poly_order))
      max_poly_order = _orbit.prototype().size();

    //std::cout<<"Max_poly_order "<<max_poly_order<<std::endl;

    std::vector<BasisSet const *> arg_subsets;
    for(DoFKey const &key : global_keys) {
      auto find_it = m_global_bases.find(key);
      if(find_it != m_global_bases.end())
        arg_subsets.push_back(&(find_it->second));
      else
        throw std::runtime_error("Unable to construct basis sets. No known global DoF: " + key + "\n");
    }

    // copy local site bases to
    std::vector<BasisSet> all_local;
    all_local.reserve(local_keys.size());

    //Loop over dof's
    for(DoFKey const &key : local_keys) {
      // Make copies of local arguments to ensure that they are distinguishable by their DoF_IDs
      // i.e., make copies in 'tlocal' and reset the DoF_IDs to {0,1,2,etc...}
      auto find_it = m_site_bases.find(key);
      if(find_it == m_site_bases.end())
        throw std::runtime_error("Unable to construct basis sets. No known local DoF: " + key + "\n");

      std::vector<BasisSet> const &arg_vec(find_it->second);
      std::vector<BasisSet> tlocal;
      tlocal.reserve(_orbit.prototype().size());
      std::vector<BasisSet const *> site_args(_orbit.prototype().size(), nullptr);
      //Loop over sites
      for(Index i = 0; i < _orbit.prototype().size(); i++) {
        if(arg_vec[_orbit.prototype()[i].sublat()].size()) {
          tlocal.push_back(arg_vec[_orbit.prototype()[i].sublat()]);
          tlocal.back().set_dof_IDs(std::vector<Index>(1, i));
          site_args[i] = &tlocal.back();
        }
      }
      all_local.push_back(ClexBasis_impl::construct_clust_dof_basis(_orbit.prototype(), site_args));
      if(all_local.back().size())
        arg_subsets.push_back(&(all_local.back()));
    }

    return m_basis_builder->build(_orbit.prototype(), arg_subsets, max_poly_order, 1);
  }

  //********************************************************************

  // Divide by multiplicity. Same result as evaluating correlations via orbitree.
  template<typename OrbitType>
  std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                      OrbitType const &_clust_orbit,
                                                      PrimNeighborList &_nlist,
                                                      std::vector<FunctionVisitor *> const &labelers) {
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

  //********************************************************************
  /// nlist_index is the index of the basis site in the neighbor list
  template<typename OrbitType>
  std::map< UnitCell, std::vector< std::string > > flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                                               std::function<BasisSet(BasisSet const &)> _bset_transform,
                                                                               OrbitType const &_clust_orbit,
                                                                               PrimNeighborList &_nlist,
                                                                               std::vector<FunctionVisitor *> const &labelers,
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

  //***********************************************
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
    out << "Basis site definitions for DoF " << it->first << ".\n";
    Index nb;
    for(Index nb = 0; nb < (it->second).size(); ++nb) {
      out << "  Basis site " << nb + 1 << ":\n"
          << "  ";
      _prim.basis[nb].print(out);
      out << "\n";
      //HERE
      out << DoFType::traits(it->first).site_basis_description((it->second)[nb], _prim.basis[nb]);
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

  //***********************************************
  namespace ClexBasis_impl {
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
