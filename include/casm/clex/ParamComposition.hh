#ifndef PARAMCOMPOSITION_HH
#define PARAMCOMPOSITION_HH

#include "casm/external/boost.hh"
using boost::property_tree::ptree;

#include "casm/container/Array.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  class Structure;

  //defines the enum type used in composition.
  enum COMPOSITION_TYPE {PARAM_COMP = 0, NUMBER_ATOMS = 1};

  class ParamComposition {
    // holds the transformation matrix to go from NUMBER_ATOMS to
    // PARAM_COMP and vice versa
    Array< Eigen::MatrixXd > m_comp;

    // hold the list of allowed components, based on the PRIM
    Array< std::string > m_components;

    // The origin of the composition space
    Eigen::VectorXd m_origin;

    //the number of variables you need to define to give the position
    //in this space
    int m_rank_of_space;

    // The prim that we are considering compositions of
    const Structure *m_prim_struc;

    //holds the list of end members as defined in this space by the
    //comp and origin matrices
    Array< Eigen::VectorXd > m_spanning_end_members;

    // holds the list of all allowed end_members in the PRIM
    //   each row is an end_member
    Eigen::MatrixXd m_prim_end_members;

    // gives the sublattices that components[i] can be allowed on  [component][sublat] (0/1)
    Eigen::MatrixXi m_sublattice_map;

    // the list of possible composition axes that have positive
    // composition axes as computed by generate_composition_axes
    Array< ParamComposition > m_allowed_list;
  public:
    //*************************************************************
    //CONSTRUCTORS
    ParamComposition() : m_prim_end_members(0, 0), m_sublattice_map(0, 0) {
      std::cerr << "WARNING in ParamComposition, you have not initialized a PrimClex. Things could go horribly wrong if you dont. I hope you know what you are doing." << std::endl;
      m_comp.resize(2);
      m_rank_of_space = -1;
      m_comp[0].resize(0, 0);
      m_comp[1].resize(0, 0);
      m_origin.resize(0);
      return;
    };

    ParamComposition(const Structure &_prim) : m_prim_end_members(0, 0), m_sublattice_map(0, 0) {
      m_comp.resize(2);
      m_comp[0].resize(0, 0);
      m_comp[1].resize(0, 0);
      m_origin.resize(0);
      m_prim_struc = &_prim;
      m_rank_of_space = -1;
    }

    ParamComposition(const Array< std::string > &_components, const Eigen::MatrixXd &transf_mat, const Eigen::VectorXd &_origin, const int &_rank_of_space, const Structure &_prim, const int &COMP_TYPE) {
      m_rank_of_space = _rank_of_space;
      m_components = _components;
      m_origin = _origin;
      m_prim_struc = &_prim;
      m_comp.resize(2);
      if(COMP_TYPE == PARAM_COMP) {
        m_comp[PARAM_COMP] = transf_mat;
        m_comp[NUMBER_ATOMS] = m_comp[PARAM_COMP].inverse();
      }
      else if(COMP_TYPE == NUMBER_ATOMS) {
        m_comp[NUMBER_ATOMS] = transf_mat;
        m_comp[PARAM_COMP] = m_comp[NUMBER_ATOMS].inverse();
      }
      calc_spanning_end_members();

    };

    /*
        // Add in the primclex pointer
        ParamComposition(const std::string &json_filename, const Structure &_prim) {
          m_prim_struc = &_prim;
          m_comp.resize(2);
          m_comp[0].resize(0, 0);
          m_comp[1].resize(0, 0);
          m_origin.resize(0);
          std::ifstream comp_json_file;
          comp_json_file.open(json_filename.c_str());
          if(!comp_json_file) {
            std::cout << "ERROR\n";
            return;
          }
          read(comp_json_file);
          return;
        };

        ParamComposition(ptree comp_ptree, const Structure &_prim) {
          m_prim_struc = &_prim;
          m_comp.resize(2);
          m_comp[0].resize(0, 0);
          m_comp[1].resize(0, 0);
          m_origin.resize(0);
          read(comp_ptree);
        }
    */
    //*************************************************************
    //GENERATE Routines
    void generate_components();
    void generate_prim_end_members();
    //void generate_composition_axes();
    void generate_composition_transf();
    void generate_sublattice_map();
    void generate_composition_space(bool verbose = false);

    //*************************************************************
    //CALC Routines
    ParamComposition calc_composition_object(const Eigen::VectorXd &torigin, const Array< Eigen::VectorXd> tspanning);
    Eigen::VectorXd calc(const Eigen::VectorXd &tcomp, const int &MODE);
    //    ptree calc_composition_ptree() const; //returns a ptree object with all the data from this composition object
    void calc_spanning_end_members();
    void calc_transformation_matrices();

    Eigen::VectorXd calc_param_composition(const Eigen::VectorXd &num_atoms_per_prim) const;
    Eigen::VectorXd calc_num_atoms(const Eigen::VectorXd &param_composition) const;

    // Lists components (species) of crystal whose compositions are fixed (i.e., are not involved in alloying)
    // each pair gives (species_name, #in_prim)
    std::vector<std::pair<std::string, Index> > fixed_components();
    //*************************************************************
    //PRINT Routines

    //    void print(std::ostream &stream, bool print_comp_axes_flag = false) const {
    //      if(print_comp_axes_flag) {
    //        stream << "/*" << std::endl;
    //        print_composition_axes(stream);
    //        stream << "*/" << std::endl;
    //      }
    //      ptree comp_ptree = calc_composition_ptree();
    //      write_json(stream, comp_ptree);
    //    }

    void print_sublattice_map(std::ostream &stream) const {
      stream << "SUBLATTICE MAP" << std::endl << "-------------" << std::endl;
      stream << m_sublattice_map << std::endl;
    }
    void print_prim_end_members(std::ostream &stream) const {
      stream << "Number of end members: " << m_prim_end_members.rows() << std::endl;
      stream << m_components << std::endl << "------" << std::endl;
      stream << m_prim_end_members << std::endl;
    }
    void print_components(std::ostream &stream) const {
      stream << "Components: " << m_components << std::endl;
    };

    void print_composition_axes(std::ostream &stream) const;
    void print_curr_composition_axes(std::ostream &stream) const;
    void print_end_member_formula(const int &end_member_index, std::ostream &stream, const int &stream_width) const;
    void print_member_formula(const Eigen::VectorXd &member, std::ostream &stream, const int &stream_width) const;
    void print_origin_formula(std::ostream &stream, const int &stream_width) const;
    void print_composition_formula(std::ostream &stream, const int &stream_width) const;

    void print_composition_matrices(std::ostream &stream) const {
      stream << "components: ";
      print_components(stream);
      stream << "comp[PARAM_COMP] " << m_comp[PARAM_COMP] << std::endl;
      stream << "comp[NUMBER_ATOMS] " << m_comp[NUMBER_ATOMS] << std::endl;
      stream << "origin: " << m_origin << std::endl;
    }

    //*************************************************************
    //READ
    //    void read(const std::string &comp_filename);
    //    void read(std::istream &stream);
    //    void read(ptree comp_ptree);

    //*************************************************************
    //MISCELLANEOUS
    void max_out(const int &component_index, Eigen::MatrixXi &sublat_comp) const;
    void select_composition_axes(const Index &choice);

    //*************************************************************
    //ACCESSORS

    const Structure &prim() const {
      return *m_prim_struc;
    }

    const Array< Eigen::VectorXd > &spanning_end_members() const {
      return m_spanning_end_members;
    };

    /// \brief Return all possible end members as row matrix
    Eigen::MatrixXd prim_end_members() const {
      return m_prim_end_members;
    }

    const Array< Eigen::MatrixXd > &comp() const {
      return m_comp;
    };

    const Eigen::VectorXd &origin() const {
      return m_origin;
    };

    const int &rank_of_space() const {
      return m_rank_of_space;
    };

    const int &number_of_references() const {
      return m_rank_of_space;
    };

    /// \brief Components are ordered as in Structure::struc_molecule
    const Array<std::string> &components() const {
      return m_components;
    };

    std::string composition_formula() const;

    const Array<ParamComposition> &allowed_list() const {
      return m_allowed_list;
    }

    //*************************************************************
    //TEST FUNCTIONS

    //A ParamComposition is defined as 'set' if origin is non-empty, the
    //prim pointer does not point to nullptr, rank_of_space is greater
    //than 0, the comp matrices are non-empty and have their matrices
    //such that they are square
    bool is_set() const {
      if(m_rank_of_space <= 0 || m_components.size() == 0 || m_origin.size() != m_components.size()) {
        return false;
      }
      if(m_comp.size() == 2) {
        return (m_comp[0].rows() == m_components.size() && m_comp[0].cols() == m_components.size() &&
                m_comp[1].rows() == m_components.size() && m_comp[1].cols() == m_components.size());
      }
      else {
        return false;
      }
      return true;
    };
  };

}
#endif
