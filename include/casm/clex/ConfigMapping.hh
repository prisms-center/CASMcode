#ifndef CONFIGMAPPING_HH
#define CONFIGMAPPING_HH
#include "casm/CASM_global_definitions.hh"
namespace CASM {
  class Supercell;
  class Lattice;
  class SymGroup;
  class PrimClex;

  class ConfigMapper{
  public:
    enum NullInitializer{null_initializer};
    ConfigMapper(NullInitializer) : m_pclex(NULL) {}
    ConfigMapper(PrimClex &_pclex,double _lattice_weight, double _max_vol_change=0.25, double _tol=TOL);


    PrimClex& primclex() const{
      return *m_pclex;
    }

    double lattice_weight() const{
      return m_lattice_weight;
    }
    bool import_structure_occupation(const fs::path &pos_path,
                                     std::string &imported_name,
                                     jsonParser &relaxation_properties,
                                     bool robust_flag,
                                     bool rotate_flag,
                                     bool strict_flag,
                                     double vol_tol = 0.25) const;

    bool import_structure_occupation(const BasicStructure<Site> &_struc,
                                     std::string &imported_name,
                                     jsonParser &relaxation_properties,
                                     bool robust_flag,
                                     bool rotate_flag,
                                     bool strict_flag,
                                     double vol_tol = 0.25) const;

    bool import_structure_occupation(const BasicStructure<Site> &_struc,
                                     const Configuration *hint_ptr,
                                     std::string &imported_name,
                                     jsonParser &relaxation_properties,
                                     bool robust_flag,
                                     bool rotate_flag,
                                     bool strict_flag,
                                     double vol_tol = 0.25) const;

    bool import_structure(const fs::path &pos_path,
                          std::string &imported_name,
                          jsonParser &relaxation_properties,
                          bool robust_flag,
                          bool rotate_flag,
                          bool strict_flag,
                          double vol_tol = 0.25) const;

    bool import_structure(const BasicStructure<Site> &_struc,
                          std::string &imported_name,
                          jsonParser &relaxation_properties,
                          bool robust_flag,
                          bool rotate_flag,
                          bool strict_flag,
                          double vol_tol = 0.25) const;

    bool struc_to_configdof(const BasicStructure<Site> &_struc,
                            ConfigDoF &mapped_configdof,
                            Lattice &mapped_lat,
                            bool robust_flag,
                            bool rotate_flag,
                            double vol_tol = 0.25) const;


    bool ideal_struc_to_configdof(BasicStructure<Site> struc,
                                  ConfigDoF &mapped_config_dof,
                                  Lattice &mapped_lat) const;


    bool deformed_struc_to_configdof(const BasicStructure<Site> &_struc,
                                     ConfigDoF &mapped_config_dof,
                                     Lattice &mapped_lat,
                                     bool rotate_flag,
                                     double vol_tol = 0.25) const;
  
  private:
    PrimClex * m_pclex;
    mutable std::map<Index, std::vector<Lattice> > m_superlat_map;
    double m_lattice_weight, m_max_volume_change;
    double m_tol;
    std::vector<std::pair<std::string,Index> > m_fixed_components;
    const std::vector<Lattice>& _lattices_of_vol(Index prim_vol) const;
  };
  
  Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                     const Lattice &relaxed_lat,
                                     const SymGroup &sym_group,
                                     Eigen::MatrixXd &trans_mat,
                                     Eigen::MatrixXd &deformation,
                                     double &best_cost,
                                     Index min_vol,
                                     Index max_vol,
                                     double _tol);

  Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                     const Lattice &relaxed_lat,
                                     const SymGroup &sym_group,
                                     Eigen::MatrixXd &trans_mat,
                                     Eigen::MatrixXd &deformation,
                                     double &best_cost,
                                     const std::vector<Lattice> &from_range,
                                     double _tol);
  namespace ConfigMap_impl{

    // Assignment Problem Routines
    // Find cost matrix for displacements between POS and relaxed structures.
    // Returns false if lattices are incompatible
    bool calc_cost_matrix(const Supercell &scel,
                          const BasicStructure<Site> &rstruc,
                          const Coordinate &trans,
                          Eigen::MatrixXd &cost_matrix);
    
    //\JSB
    
    
    
    // mapping routine. Return an ideal configuration corresponding to a relaxed structure.
    // Return false if 'rstruc' is incompatible with supercell (can happen frequently when vacancies are allowed)
    // Options:
    //   TRANSLATE = true -> rigid-translations are removed. (typically this option should be used, especially if you care about vacancies)
    //
    //   TRANSLATE = false -> rigid translations are not considered. (less robust but more efficient -- use only if you know rigid translations are small or zero)
    
    bool struc_to_configdof(const Supercell &scel,
                            BasicStructure<Site> rstruc,
                            ConfigDoF &config_dof,
                            const bool translate_flag,
                            const double _tol);
  }
  namespace ConfigMapping {
    double strain_cost(const Lattice &relaxed_lat, const ConfigDoF &_dof);

    double basis_cost(const ConfigDoF &_dof);

  }


}

#endif
