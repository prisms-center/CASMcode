#ifndef CONFIGMAPPING_HH
#define CONFIGMAPPING_HH
#include "casm/CASM_global_definitions.hh"
namespace CASM {
  class Supercell;
  class Lattice;
  class SymGroup;
  class PrimClex;

  Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                     const Lattice &relaxed_lat,
                                     const SymGroup &sym_group,
                                     Eigen::MatrixXd &trans_mat,
                                     Eigen::MatrixXd &deformation,
                                     double &best_cost,
                                     Index min_vol,
                                     Index max_vol,
                                     double _tol);

  bool import_structure_occupation(const fs::path &pos_path,
                                   PrimClex &pclex,
                                   std::string &imported_name,
                                   jsonParser &relaxation_properties,
                                   bool robust_flag,
                                   bool rotate_flag,
                                   bool strict_flag,
                                   double _tol,
                                   double lattice_weight = 0.5,
                                   double vol_tol = 0.25);

  bool import_structure_occupation(const BasicStructure<Site> &_struc,
                                   PrimClex &pclex,
                                   std::string &imported_name,
                                   jsonParser &relaxation_properties,
                                   bool robust_flag,
                                   bool rotate_flag,
                                   bool strict_flag,
                                   double _tol,
                                   double lattice_weight = 0.5,
                                   double vol_tol = 0.25);

  bool import_structure_occupation(const BasicStructure<Site> &_struc,
                                   const Configuration *hint_ptr,
                                   PrimClex &pclex,
                                   std::string &imported_name,
                                   jsonParser &relaxation_properties,
                                   bool robust_flag,
                                   bool rotate_flag,
                                   bool strict_flag,
                                   double _tol,
                                   double lattice_weight = 0.5,
                                   double vol_tol = 0.25);

  bool import_structure(const fs::path &pos_path,
                        PrimClex &pclex,
                        std::string &imported_name,
                        jsonParser &relaxation_properties,
                        bool robust_flag,
                        bool rotate_flag,
                        bool strict_flag,
                        double _tol,
                        double lattice_weight = 0.5,
                        double vol_tol = 0.25);

  bool import_structure(const BasicStructure<Site> &_struc,
                        PrimClex &pclex,
                        std::string &imported_name,
                        jsonParser &relaxation_properties,
                        bool robust_flag,
                        bool rotate_flag,
                        bool strict_flag,
                        double _tol,
                        double lattice_weight = 0.5,
                        double vol_tol = 0.25);

  bool struc_to_configdof(const BasicStructure<Site> &_struc,
                          PrimClex &pclex,
                          ConfigDoF &mapped_configdof,
                          Lattice &mapped_lat,
                          bool robust_flag,
                          bool rotate_flag,
                          double _tol,
                          double lattice_weight = 0.5,
                          double vol_tol = 0.25);


  bool ideal_struc_to_configdof(BasicStructure<Site> struc,
                                PrimClex &pclex,
                                ConfigDoF &mapped_config_dof,
                                Lattice &mapped_lat,
                                double _tol);


  bool deformed_struc_to_configdof(const BasicStructure<Site> &_struc,
                                   PrimClex &pclex,
                                   ConfigDoF &mapped_config_dof,
                                   Lattice &mapped_lat,
                                   bool rotate_flag,
                                   double _tol,
                                   double lattice_weight = 0.5,
                                   double vol_tol = 0.25);



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

  namespace ConfigMapping {
    double strain_cost(const Lattice &relaxed_lat, const ConfigDoF &_dof);

    double basis_cost(const ConfigDoF &_dof);

  }


}

#endif
