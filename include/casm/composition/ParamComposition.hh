#ifndef CASM_ParamComposition
#define CASM_ParamComposition

#include <vector>

#include "casm/casm_io/container/stream_io.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
// defines the enum type used in composition.
enum COMPOSITION_TYPE { PARAM_COMP = 0, NUMBER_ATOMS = 1 };

class ParamComposition {
 public:
  using AllowedOccupants = std::vector<std::vector<std::string> >;

  static std::vector<std::string> string_components(
      ParamComposition::AllowedOccupants const &_allowed_occs);
  //*************************************************************
  // CONSTRUCTORS
  ParamComposition() : m_prim_end_members(0, 0) {
    std::cerr << "WARNING in ParamComposition, you have not initialized a "
                 "PrimClex. Things could go horribly wrong if you dont. I hope "
                 "you know what you are doing."
              << std::endl;
    m_comp.resize(2);
    m_rank_of_space = -1;
    m_comp[0].resize(0, 0);
    m_comp[1].resize(0, 0);
    m_origin.resize(0);
    return;
  }

  ParamComposition(AllowedOccupants _allowed_occs);

  ParamComposition(AllowedOccupants _allowed_occs,
                   const Eigen::MatrixXd &transf_mat,
                   const Eigen::VectorXd &_origin, const int &_rank_of_space,
                   const int &COMP_TYPE);

  //*************************************************************
  // GENERATE Routines
  void generate_prim_end_members();
  void generate_composition_transf();
  void generate_composition_space(bool verbose = false);

  //*************************************************************
  // CALC Routines
  ParamComposition calc_composition_object(
      const Eigen::VectorXd &torigin,
      const std::vector<Eigen::VectorXd> tspanning);
  Eigen::VectorXd calc(const Eigen::VectorXd &tcomp, const int &MODE);
  //    ptree calc_composition_ptree() const; //returns a ptree object with all
  //    the data from this composition object
  void calc_spanning_end_members();
  void calc_transformation_matrices();

  Eigen::VectorXd calc_param_composition(
      const Eigen::VectorXd &num_atoms_per_prim) const;
  Eigen::VectorXd calc_num_atoms(
      const Eigen::VectorXd &param_composition) const;

  //*************************************************************
  // PRINT Routines
  void print_prim_end_members(std::ostream &stream) const {
    stream << "Number of end members: " << m_prim_end_members.rows()
           << std::endl;
    stream << m_components << std::endl << "------" << std::endl;
    stream << m_prim_end_members << std::endl;
  }
  void print_components(std::ostream &stream) const {
    stream << "Components: " << m_components << std::endl;
  }

  void print_composition_axes(std::ostream &stream) const;
  void print_curr_composition_axes(std::ostream &stream) const;
  void print_end_member_formula(const int &end_member_index,
                                std::ostream &stream,
                                const int &stream_width) const;
  void print_member_formula(const Eigen::VectorXd &member, std::ostream &stream,
                            const int &stream_width) const;
  void print_origin_formula(std::ostream &stream,
                            const int &stream_width) const;
  void print_composition_formula(std::ostream &stream,
                                 const int &stream_width) const;

  void print_composition_matrices(std::ostream &stream) const {
    stream << "components: ";
    print_components(stream);
    stream << "comp[PARAM_COMP] " << m_comp[PARAM_COMP] << std::endl;
    stream << "comp[NUMBER_ATOMS] " << m_comp[NUMBER_ATOMS] << std::endl;
    stream << "origin: " << m_origin << std::endl;
  }

  //*************************************************************
  // READ
  //    void read(const std::string &comp_filename);
  //    void read(std::istream &stream);
  //    void read(ptree comp_ptree);

  //*************************************************************
  // MISCELLANEOUS
  void select_composition_axes(const Index &choice);

  //*************************************************************
  // ACCESSORS

  const std::vector<Eigen::VectorXd> &spanning_end_members() const {
    return m_spanning_end_members;
  }

  /// \brief Return all possible end members as row matrix
  Eigen::MatrixXd prim_end_members() const { return m_prim_end_members; }

  const std::vector<Eigen::MatrixXd> &comp() const { return m_comp; }

  const Eigen::VectorXd &origin() const { return m_origin; }

  const int &rank_of_space() const { return m_rank_of_space; }

  const int &number_of_references() const { return m_rank_of_space; }

  /// \brief For each sublattice, a list of occupants allowed on that sublattice
  const std::vector<std::vector<std::string> > &allowed_occs() const {
    return m_allowed_occs;
  }

  /// \brief Components are in order of appearance precedence in allowed_occs()
  const std::vector<std::string> &components() const { return m_components; }

  std::string composition_formula() const;

  const std::vector<ParamComposition> &allowed_list() const {
    return m_allowed_list;
  }

  //*************************************************************
  // TEST FUNCTIONS

  // A ParamComposition is defined as 'set' if origin is non-empty, the
  // prim pointer does not point to nullptr, rank_of_space is greater
  // than 0, the comp matrices are non-empty and have their matrices
  // such that they are square
  bool is_set() const {
    if (m_rank_of_space <= 0 || m_components.size() == 0 ||
        m_origin.size() != m_components.size()) {
      return false;
    }
    if (m_comp.size() == 2) {
      return (m_comp[0].rows() == m_components.size() &&
              m_comp[0].cols() == m_components.size() &&
              m_comp[1].rows() == m_components.size() &&
              m_comp[1].cols() == m_components.size());
    } else {
      return false;
    }
    return true;
  }

 private:
  /// \brief m_allowed_occs[b] is list of occupants allowed at site [b]
  AllowedOccupants m_allowed_occs;

  /// The list of all allowed components, based on allowed_occs()
  std::vector<std::string> m_components;

  /// holds the transformation matrix to go from NUMBER_ATOMS to
  /// PARAM_COMP and vice versa
  std::vector<Eigen::MatrixXd> m_comp;

  /// The origin of the composition space
  Eigen::VectorXd m_origin;

  /// the number of variables you need to define to give the position
  /// in this space
  int m_rank_of_space;

  /// holds the list of end members as defined in this space by the
  /// comp and origin matrices
  std::vector<Eigen::VectorXd> m_spanning_end_members;

  /// holds the list of all allowed end_members in the PRIM
  ///   each row is an end_member
  Eigen::MatrixXd m_prim_end_members;

  /// the list of possible composition axes that have positive
  /// composition axes as computed by generate_composition_axes
  std::vector<ParamComposition> m_allowed_list;
};

}  // namespace CASM
#endif
