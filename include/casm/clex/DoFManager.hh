#ifndef DOFMANAGER_HH
#define DOFMANAGER_HH

#include "casm/container/Array.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  class DoFEnvironment;

  ///DoFManager holds multiple DoFEnvironments, and provides a simple interface for adding and managing DoFs.
  class DoFManager {
    Array<DoFEnvironment *> m_environs;
  public:
    DoFManager() {};
    DoFManager(const DoFManager &RHS);
    ~DoFManager();
    bool contains(const std::string &dof_name) const;
    void add_dof(const std::string &dof_name);
    void set_local_dof_state(const Configuration &config, Index l);
    void set_global_dof_state(const Configuration &config);

    void resize_neighborhood(Index nlist_size);

    template<typename ClustType>
    void register_dofs(GenericOrbitree<ClustType> &tree)const;

    //Clexulator printing routines
    ReturnArray<FunctionVisitor *> get_function_label_visitors() const;
    void print_clexulator_member_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;
    void print_clexulator_private_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;
    void print_clexulator_private_method_implementations(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;

    void print_clexulator_public_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;
    void print_clexulator_public_method_implementations(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;

    void print_to_clexulator_constructor(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;
  };

  /// DoFEnvironment is abstract class that handles interface between Configurations and BasisSets
  /// and allows for efficient evaluation of functions on a Configuration.
  /// Evaluation is structured this way to preserve polymorphism of DoF and allow easy addition of new types of DoFs.
  /// At the same time, we wish to prevent circular dependencies among DoF, Configuration, BasisSet and Function.
  class DoFEnvironment {
    std::string m_dof_name;

  public:
    DoFEnvironment(std::string dof_name) : m_dof_name(dof_name) {};
    virtual ~DoFEnvironment() {};

    const std::string &name() const {
      return m_dof_name;
    };

    ///set the state of GLOBAL parameters based on their value in config
    virtual void set_global_state(const Configuration &config) = 0;

    ///set the state of LOCAL parameters based on their value in config
    virtual void set_local_state(const Configuration &config, Index l) = 0;

    virtual void resize_neighborhood(Index nlist_size) {};

    /// register_dofs sizes internal datastructures and then loops over all Orbitrees to register
    /// the remote values in each DoF whose name matches m_dof_name
    virtual int register_dofs(BasisSet &basis)const = 0;

    // Help with Clexulator printing
    virtual FunctionVisitor *get_function_label_visitor() const {
      return nullptr;
    };

    virtual void print_clexulator_member_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const {};

    virtual void print_clexulator_private_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const {};
    virtual void print_clexulator_private_method_implementations(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const {};

    virtual void print_clexulator_public_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const {};
    virtual void print_clexulator_public_method_implementations(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const {};

    virtual void print_to_clexulator_constructor(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const {};
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class OccupationDoFEnvironment : public DoFEnvironment {
    //inherited from DoFEnvironment:
    //  std::string m_dof_name;

    Array<int> m_neighbor_occ;
  public:
    OccupationDoFEnvironment(std::string dof_name = "p") : DoFEnvironment(dof_name) {};

    void set_global_state(const Configuration &config);
    void set_local_state(const Configuration &config, Index l);

    void resize_neighborhood(Index nlist_size) {
      m_neighbor_occ.resize(nlist_size, 0);
    };

    int register_dofs(BasisSet &basis) const {
      Array<DoF::RemoteHandle> handles;
      for(Index i = 0; i < m_neighbor_occ.size(); i++)
        handles.push_back(DoF::RemoteHandle(m_neighbor_occ[i]));
      return basis.register_remotes(name(), handles);
    };

    FunctionVisitor *get_function_label_visitor() const {
      //std::cout << "*****************MAKING OCC_FUNC LABELER\n";
      return new OccFuncLabeler("occ_func_%b_%f(%n)");
    };

    void print_clexulator_member_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const;

    void print_clexulator_private_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;

    void print_clexulator_public_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;

    void print_to_clexulator_constructor(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DisplacementDoFEnvironment : public DoFEnvironment {
    //inherited from DoFEnvironment:
    //  std::string m_dof_name;

    // disp_index indicates which displacement component of Vector3<double> to select
    // e.g., to get the value of y_l (e.g., the y-displacement at site l), disp_index should be 1,
    // and the value is at config.get_disp(l)[disp_index]
    int m_disp_index;

    Array<double> m_neighbor_disp;
  public:
    DisplacementDoFEnvironment(std::string dof_name) : DoFEnvironment(dof_name) {
      if(dof_name == "disp_x")
        m_disp_index = 0;
      else if(dof_name == "disp_y")
        m_disp_index = 1;
      else if(dof_name == "disp_z")
        m_disp_index = 2;
      else
        m_disp_index = -1;
    };

    void set_global_state(const Configuration &config);
    void set_local_state(const Configuration &config, Index l);


    void resize_neighborhood(Index nlist_size) {
      m_neighbor_disp.resize(nlist_size, 0);
    };

    int register_dofs(BasisSet &basis)const {
      Array<DoF::RemoteHandle> handles;
      for(Index i = 0; i < m_neighbor_disp.size(); i++)
        handles.push_back(DoF::RemoteHandle(m_neighbor_disp[i]));
      return basis.register_remotes(name(), handles);

    };

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  class StrainDoFEnvironment : public DoFEnvironment {
    //inherited from DoFEnvironment:
    //  std::string m_dof_name;

    Array<double> m_strain_vals;
  public:
    StrainDoFEnvironment(std::string dof_name) : DoFEnvironment(dof_name) {};

    void set_global_state(const Configuration &config);
    void set_local_state(const Configuration &config, Index l);


    int register_dofs(BasisSet &basis)const {
      Array<DoF::RemoteHandle> handles;
      for(Index i = 0; i < m_strain_vals.size(); i++)
        handles.push_back(DoF::RemoteHandle(m_strain_vals[i]));
      return basis.register_remotes(name(), handles);
    };

  };


  //************************************************************

  template<typename ClustType>
  void DoFManager::register_dofs(GenericOrbitree<ClustType> &tree) const {
    int sum(0);
    Index nb, no, ne, i;
    for(nb = 0; nb < tree.size(); nb++) {
      for(no = 0; no < tree[nb].size(); no++) {
        for(ne = 0; ne < tree[nb][no].size(); ne++) {

          for(i = 0; i < m_environs.size(); i++)
            sum += m_environs[i]->register_dofs(tree[nb][no][ne].clust_basis);

        }
      }
    }
    //std::cout << "Finished registering " << sum << " remotes.\n";
  }

};

#endif
