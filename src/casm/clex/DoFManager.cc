#include "casm/clex/DoFManager.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {


  DoFManager::DoFManager(const DoFManager &RHS) {
    for(Index i = 0; i < m_environs.size(); i++) {
      add_dof(m_environs[i]->name());
    }
  }

  //************************************************************

  DoFManager::~DoFManager() {
    for(Index i = 0; i < m_environs.size(); i++) {
      delete m_environs[i];
    }
  }

  //************************************************************
  bool DoFManager::contains(const std::string &dof_name) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      if(dof_name == (m_environs[i]->name()))
        return true;
    }
    return false;
  }
  //************************************************************
  void DoFManager::add_dof(const std::string &dof_name) {
    if(contains(dof_name))
      return;

    if(dof_name == "occupation")
      m_environs.push_back(new OccupationDoFEnvironment(dof_name));
    else if(dof_name == "x-displacemnt")
      m_environs.push_back(new DisplacementDoFEnvironment(dof_name, 0));
    else if(dof_name == "y-displacemnt")
      m_environs.push_back(new DisplacementDoFEnvironment(dof_name, 1));
    else if(dof_name == "z-displacemnt")
      m_environs.push_back(new DisplacementDoFEnvironment(dof_name, 2));
    else if(dof_name == "strain")
      m_environs.push_back(new StrainDoFEnvironment(dof_name));
    else {
      std::cerr << "CRITICAL ERROR: Unknown option '" << dof_name << "' passed to DoFManager::add_dof()!\n"
                << "                Exiting...\n";
      exit(1);
    }
    return;
  }

  //************************************************************

  void DoFManager::set_local_dof_state(const Configuration &config, Index l) {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->set_local_state(config, l);
    }
  }

  //************************************************************

  void DoFManager::set_global_dof_state(const Configuration &config) {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->set_global_state(config);
    }
  }

  //************************************************************

  void DoFManager::resize_neighborhood(Index nlist_size) {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->resize_neighborhood(nlist_size);
    }
  }

  //************************************************************

  ReturnArray<FunctionVisitor *> DoFManager::get_function_label_visitors() const {
    Array<FunctionVisitor *> tlabels;
    for(Index i = 0; i < m_environs.size(); i++) {
      tlabels.push_back(m_environs[i]->get_function_label_visitor());
      if(!tlabels.back())
        tlabels.pop_back();
    }
    return tlabels;
  }

  //************************************************************

  void DoFManager::print_clexulator_member_definitions(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->print_clexulator_member_definitions(stream, primclex, indent);
    }
  }

  //************************************************************

  void DoFManager::print_clexulator_private_method_definitions(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->print_clexulator_private_method_definitions(stream, primclex, indent);
    }
  }

  //************************************************************

  void DoFManager::print_clexulator_private_method_implementations(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->print_clexulator_private_method_implementations(stream, primclex, indent);
    }
  }

  //************************************************************

  void DoFManager::print_clexulator_public_method_definitions(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->print_clexulator_public_method_definitions(stream, primclex, indent);
    }
  }

  //************************************************************

  void DoFManager::print_clexulator_public_method_implementations(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->print_clexulator_public_method_implementations(stream, primclex, indent);
    }
  }

  //************************************************************

  void DoFManager::print_to_clexulator_constructor(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index i = 0; i < m_environs.size(); i++) {
      m_environs[i]->print_to_clexulator_constructor(stream, primclex, indent);
    }
  }

  //************************************************************

  void OccupationDoFEnvironment::print_clexulator_member_definitions(std::ostream &stream, const PrimClex &primclex, const std::string &indent)const {
    for(Index b = 0; b < primclex.get_prim().basis.size(); b++) {
      if(!primclex.get_prim().basis[b].occupant_basis().size())
        continue;

      stream <<
             indent << "// Occupation Function table for basis site " << b << ":\n";
      for(Index f = 0; f < primclex.get_prim().basis[b].occupant_basis().size(); f++) {
        stream <<
               indent << "double " << "m_occ_func_" << b << '_' << f << '[' << primclex.get_prim().basis[b].site_occupant().size() << "];\n";
      }
      stream << '\n';
    }

  }

  //************************************************************
  void OccupationDoFEnvironment::print_clexulator_private_method_definitions(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    for(Index b = 0; b < primclex.get_prim().basis.size(); b++) {
      if(!primclex.get_prim().basis[b].occupant_basis().size())
        continue;

      stream <<
             indent << "// Occupation Function accessors for basis site " << b << ":\n";
      for(Index f = 0; f < primclex.get_prim().basis[b].occupant_basis().size(); f++) {
        stream <<
               indent << "const double &occ_func_" << b << '_' << f << "(const int &nlist_ind)const{return " << "m_occ_func_" << b << '_' << f << "[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}\n";
      }
      stream << '\n';
    }
  }

  //************************************************************
  void OccupationDoFEnvironment::print_clexulator_public_method_definitions(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {

  }

  //************************************************************

  void OccupationDoFEnvironment::print_to_clexulator_constructor(std::ostream &stream, const PrimClex &primclex, const std::string &indent) const {
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream.precision(10);

    for(Index b = 0; b < primclex.get_prim().basis.size(); b++) {
      for(Index f = 0; f < primclex.get_prim().basis[b].occupant_basis().size(); f++) {
        for(Index s = 0; s < primclex.get_prim().basis[b].site_occupant().size(); s++) {
          if(s == 0)
            stream << indent;
          stream << "m_occ_func_" << b << '_' << f << '[' << s << "] = "
                 << primclex.get_prim().basis[b].occupant_basis()[f]->eval(Array<Index>(1, primclex.get_prim().basis[b].site_occupant().ID()), Array<Index>(1, s));
          if(s + 1 == primclex.get_prim().basis[b].site_occupant().size())
            stream << ";\n\n";
          else
            stream << ", ";
        }
      }
    }
  }

  //************************************************************

  void OccupationDoFEnvironment::set_global_state(const Configuration &config) {
    //do nothing for global set_state
  }

  //************************************************************

  void OccupationDoFEnvironment::set_local_state(const Configuration &config, Index l) {
    const Supercell &scel(config.get_supercell());
    for(Index n = 0; n < m_neighbor_occ.size(); n++) {
      m_neighbor_occ[n] = config.occ(scel.get_nlist_l(l, n));
    }
    return;
  }

  //************************************************************

  void DisplacementDoFEnvironment::set_global_state(const Configuration &config) {
    //do nothing for global set_state
  }

  //************************************************************

  void DisplacementDoFEnvironment::set_local_state(const Configuration &config, Index l) {
    const Supercell &scel(config.get_supercell());
    for(Index n = 0; n < m_neighbor_disp.size(); n++) {
      m_neighbor_disp[n] = config.disp(scel.get_nlist_l(l, n))[m_disp_index];
    }
    return;
  }

  //************************************************************

  void StrainDoFEnvironment::set_global_state(const Configuration &config) {
    //something like
    /*
    for(Index i=0; i<m_strain_vals.size(); i++){
      m_strain_vals[i]=config.get_unrolled_strain_metric()[i];
    }
    return;
    */
  }
  //************************************************************

  void StrainDoFEnvironment::set_local_state(const Configuration &config, Index l) {
    //Do nothing for local state

  }

};

