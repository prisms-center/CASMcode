#ifndef SITECLUSTER_HH
#define SITECLUSTER_HH

#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"
#include "casm/container/Tensor.hh"
#include "casm/clusterography/Cluster.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {

  class Supercell;
  class Site;
  class SymOp;
  class Permutation;
  class jsonParser;

  // For future consideration - maybe we shouldn't allow non-const access of sites
  class SiteCluster : public GenericCluster<Site> {
    Array<Array<Index> > m_trans_nlist_inds;
  public:
    BasisSet clust_basis;

    //~~~~~~~~~

    SiteCluster(const Lattice &init_home);

    SiteCluster &permute(const Array<Index> &iperm);
    SiteCluster &permute(const Permutation &perm);
    SiteCluster &apply_sym(const SymOp &op);
    SiteCluster &apply_sym_no_trans(const SymOp &op);

    void push_back(const Site &new_site);

    //void fill_discrete_basis_tensors();   //John G 011013
    void prepare_prototype();

    /// Use this method (and only this method) to set the nlist ind of each site
    /// it automatically reassigns DoF IDs for the associated basis functions.
    void set_nlist_inds(const Array<Index> &new_indices);

    /// Easily collect the current nlist_inds of the cluster's sites.
    ReturnArray<Index> nlist_inds() const;

    /// Access and assign trans_nlists
    const Array<Array<Index> > &trans_nlists() const;
    const Array<Index> &trans_nlist(Index i) const;
    void add_trans_nlist(const Array<Index> &new_nlist);

    void generate_clust_basis(Array<BasisSet const *> local_args, Array<BasisSet const *> global_args, Index max_poly_order = -1);

    inline void decorate(const Array<int> decor);
    ReturnArray<Array<int> > get_decor_map() const;
    ReturnArray<Array<int> > get_full_decor_map() const;
    ReturnArray<SiteCluster> get_decorations(const Array<Array<int> > &dmap) const;

    ///Extracts bits in bitstring corresponding to the cluster and returns them as an array
    ReturnArray<int> get_occ_array(const Array<int> &bitstring) const;

    void print_clust_basis(std::ostream &stream, int space, char delim = 0, COORD_TYPE mode = COORD_DEFAULT) const;

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  };

  SiteCluster operator*(const SymOp &LHS, const SiteCluster &RHS);

  //*************************************************
  /*
   * set site_occupant for each site according to decor
   */
  //*************************************************

  inline void SiteCluster::decorate(const Array<int> decor) {
    for(Index i = 0; i < decor.size(); i++) {
      at(i).set_occ_value(decor[i]);
    }
  };

  jsonParser &to_json(const SiteCluster &clust, jsonParser &json);
  void from_json(SiteCluster &clust, const jsonParser &json);
};

#endif
