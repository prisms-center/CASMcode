#ifndef CASM_ClexBasis_impl
#define CASM_ClexBasis_impl

#include "casm/clex/ClexBasis.hh"
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


  //*******************************************************************************************
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

}

#endif
