#ifndef DOF_HH
#define DOF_HH
#include<vector>
#include<functional>
#include<boost/algorithm/string.hpp>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymGroupRepID.hh"

namespace CASM {
  class Molecule;
  class jsonParser;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  namespace DoF_impl {
    /// \brief A class to manage dynamic evaluation of BasisFunctions

    /// A RemoteHandle can be initialized with either a double or integer reference and then passed to a
    /// DoF (or to a BasisFunction, which contains DoFs), in order for the DoF to access the remote value

    // For now this isn't used, but may be in the future if we need a dynamic clexulator class for prototyping
    class RemoteHandle {
    public:
      // \brief Initialize pointer to double
      explicit RemoteHandle(double const &d) : m_d_ptr(&d), m_i_ptr(nullptr) {}

      // \brief Initialize pointer to int
      explicit RemoteHandle(int const &i) : m_d_ptr(nullptr), m_i_ptr(&i) {}

      // \brief prevent construction by non-addressable integer
      RemoteHandle(int &&i) = delete;

      // \brief prevent construction by non-addressable double
      RemoteHandle(double &&d) = delete;

      // \brief const access of integer pointer
      const int *i_ptr() const {
        return m_i_ptr;
      }

      // \brief const access of double pointer
      const double *d_ptr() const {
        return m_d_ptr;
      }

    private:
      double const *const m_d_ptr;
      int const *const m_i_ptr;
      // In future, may need following:
      // std::string m_type_name;
      // Index m_dof_ID;
    };
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// \brief A class that represents an individual Degree of Freedom

  /// Base DoF class associates three pieces of data that identify the degree of freedom
  ///   1. type_name:     The type of DoF (e.g., occupation ("occ"), displacement ("disp"), etc.
  ///   2. variable_name: The name that will be associated with the individual DoF.  This depends
  ///                     on usage case (e.g., "x", "y", or "z" for displacement)
  ///   3. ID:            Integer subscript that identifies variables of the same name and type.
  ///                     For local DoF's, this is the site index; for global DoFs it denotes
  ///                     a subscript (e.g., for a strain component).
  ///
  /// The type_name and variable_name must be set at construction. For local DoF's, the ID can be
  /// updated as needed, but IDs of global DoFs persist throughout their lifetime.
  /// Optionally, the ID can be 'locked' which prevents it from being changed
  class DoF {
  public:
    using DoF_impl::Traits;
    using DoF_impl::RemoteHandle;

    /// DoFs are initialized using a TypeFunc
    /// e.g., DerivedDoF my_dof(DoFType::occupation);
    typedef std::function<Traits()> TypeFunc;

    static Traits const &traits(std::string const &_type_name) {
      auto it = _traits_map().find(_type_name);
      if(it == _traits_map().end()) {
        throw std::runtime_error("Could not find DoF Traits for DoF type" + _type_name);
      }
      return it->second;
    }

    /*
    DoF() :
      m_type_name("EMPTY"),
      m_var_name("EMPTY"),
      m_dof_ID(-1),
      m_ID_lock(false) {}
    */

    DoF(DoFType::TypeFunc _type_func,
        std::string const &_var_name,
        Index _ID) :
      m_type_name(_type_func().name()),
      m_var_name(_var_name),
      m_dof_ID(_ID),
      m_ID_lock(false) {
      auto it = _traits_map().find(type_name());
      if(it == _traits_map().end()) {
        _traits_map()[type_name()] = _type_func();
      }
      else {
        if(!_type_func().is_identical(it->second)) {
          throw std::runtime_error("Attempting to initialize multiple DoFs of type " + type_name() + " but having different traits!\n");
        }
      }
    }

    /// \brief allow destruction through base pointer
    /// (even though DoF shouldn't be used polymorphically)
    virtual
    ~DoF() {}

    /// \brief Const access of DoF type name
    std::string type_name() const {
      return m_type_name;
    }

    /// \brief Const access of variable name
    std::string var_name() const {
      return m_var_name;
    }

    /// \brief Const access of integer ID
    Index ID() const {
      return m_dof_ID;
    }

    /// \brief true if ID is locked
    bool is_locked() const {
      return m_ID_lock;
    }

    /// \brief mutator to set integer ID if it is unlocked
    void set_ID(Index new_ID) {
      if(m_ID_lock)
        return;
      //std::cout << "DoF " << this << " ID updated from " << ID() << " to " << new_ID << '\n';
      m_dof_ID = new_ID;
    }

    /// \brief mutator to lock integer ID
    void lock_ID() {
      m_ID_lock = true;
    }

    /// \brief mutator to unlock integer ID
    void unlock_ID() {
      m_ID_lock = false;
    }

  private:
    typedef std::map<std::string, Traits> TraitsMap;
    static TraitsMap &_traits_map();

    std::string m_type_name;
    std::string m_var_name;

    /// dof_ID is a way to distinguish between DoFs with the same name but different identities
    /// dof_ID for now usually refers to the site index of a cluster (e.g., 0, 1, 2 of a triplet)
    /// or an index into the primitive cell neighbor list. Other usage cases may be introduced later
    Index m_dof_ID;

    /// Locks the ID so that it can't be updated.  Is used for global DoF's
    bool m_ID_lock;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DiscreteDoF : public DoF {
  public:
    DiscreteDoF(): DoF(), m_current_state(0), m_remote_state(nullptr), m_sym_rep_ID(SymGroupRepID::identity(0)) {}

    DiscreteDoF(DoFType::TypeFunc _type_func,
                std::string const &_var_name,
                int _current_state = 0,
                SymGroupRepID _id = SymGroupRepID::identity(0)):
      DoF(_type_func, _var_name),
      m_current_state(_current_state),
      m_remote_state(nullptr),
      m_sym_rep_ID(_id) {

    }

    virtual
    ~DiscreteDoF() {}

    SymGroupRepID sym_rep_ID()const {
      return m_sym_rep_ID;
    }

    void set_sym_rep_ID(SymGroupRepID _id) const {
      m_sym_rep_ID = _id;
    }

    bool is_specified() const {
      return valid_index(m_current_state);
    }

    void invalidate() {
      m_current_state = -1;
    }


    int value() const {
      return m_current_state;
    }

    int remote_value() const {
      assert(m_remote_state && "In DiscreteDoF::remote_value() m_remote_state is nullptr.\n");
      return *m_remote_state;
    }

    int const *remote_ptr() const {
      return m_remote_state;
    }

    void register_remote(RemoteHandle const &handle) {
      assert(handle.i_ptr() && "In DiscreteDoF::register_remote(), attempting to register to nullptr!");
      m_remote_state = handle.i_ptr();
    }

    virtual
    DiscreteDoF *clone() const = 0;
    //      return new DiscreteDoF(*this);
    //}

    virtual
    Index size() const = 0;

    //John G 050513
    //Make an operator for this
    virtual
    void set_value(int new_state) {
      m_current_state = new_state;
      return;
    }

    virtual
    bool operator==(DiscreteDoF const &RHS) const {
      if(value() != RHS.value())
        return false;

      if(size() != RHS.size())
        return false;

      return true;
    }

    void print(std::ostream &out) const {
      out << value();
    }

    virtual
    jsonParser &to_json(jsonParser &json) const = 0;

  protected:
    /// index into domain of the current state, -1 if unspecified
    int m_current_state;

    /// Allows DoF to point to a remote value for faster/easier evaluation
    const int *m_remote_state;

    /// ID for the permutation representation for occupants
    mutable SymGroupRepID m_sym_rep_ID;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template< typename T>
  class OccupantDoF : public DiscreteDoF {
  public:
    OccupantDoF() : DiscreteDoF(DoFType::occupation,
                                  "s") { }
    OccupantDoF(std::string const &_var_name,
                std::vector<T> const &_domain,
                int _current_state = 0) :
      DiscreteDoF(DoFType::occupation,
                  _var_name,
                  _current_state,
                  SymGroupRepID::identity(_domain.size())),
      m_domain(_domain) { }

    OccupantDoF(std::vector<T> const &_domain, int _current_state = 0) :
      DiscreteDoF(DoFType::occupation,
                  "s",
                  _current_state,
                  SymGroupRepID::identity(_domain.size())),
      m_domain(_domain) { }

    OccupantDoF(std::initializer_list<T> _domain,
                int _current_state = 0) :
      DiscreteDoF(DoFType::occupation,
                  "s",
                  _current_state,
                  SymGroupRepID::identity(_domain.size())),
      m_domain(_domain.begin(),
               _domain.end()) {    }

    const T &occ() const {
      return m_domain[m_current_state];
    }

    DiscreteDoF *clone() const override {
      return new OccupantDoF(*this);
    }

    /// \brief Set occupant by index
    void set_value(int occ_index) {
      if(Index(new_state) >= m_domain.size()) {
        throw std::runtime_error("In OccupantDoF::set_value(): Bad Assignment, new_state>=size.\n");
      }
      m_current_state = new_state;
      return;
    }

    /// \brief Set occupant by value
    void set_current_state(T const &new_state) {
      Index new_index = find_index(m_domain.begin(),
                                   m_domain.end(),
                                   new_state);
      if(new_index >= m_domain.size()) {
        throw std::runtime_error("In OccupantDoF::set_current_state(): Bad Assignment, new_state not found in domain.\n");
      }
      m_current_state = new_index;
      return;
    }

    /// \brief Returns true if (*this) is equivalent to RHS, @param compare_value specifies
    /// whether to compare occupation values
    bool compare(OccupantDoF const &RHS, bool compare_value) const {
      if(compare_value && (value() != RHS.value()))
        return false;

      if(size() != RHS.size())
        return false;

      for(Index i = 0; i < size(); i++) {
        if((*this)[i] != RHS[i])
          return false;
      }

      return true;
    }

    /// \brief Equivalent to (*this).compare(RHS,true)
    bool operator==(OccupantDoF const &RHS) const {
      return compare(RHS, true);
    }

    void set_domain(std::vector<T> const &new_dom) {
      m_domain = new_dom;
      if(new_dom.size() == 1)
        m_current_state = 0;
      else
        m_current_state = -1;
      set_sym_rep_ID(SymGroupRepID::identity(new_dom.size()));
    }

    const std::vector<T> &domain() const {
      return m_domain;
    }

    T &operator[](Index i) {
      return m_domain[i];
    }

    const T &operator[](Index i)const {
      return m_domain[i];
    }

    Index size() const override {
      return m_domain.size();
    }

    void print(std::ostream &out) const {
      for(Index i = 0; i < size(); i++) {
        if(i == 0)
          out << m_domain[i].name();
        else
          out << ' ' << m_domain[i].name();
      }
    }

    void print_occ(std::ostream &out) const {
      if(valid_index(m_current_state))
        out << occ().name();
      else
        out << '?';
    }

    jsonParser &to_json(jsonParser &json) const override;

  private:
    /**** Inherited from DiscreteDoF ****
    // index into domain of the current state, -1 if unspecified
    //  int m_current_state;

    // Allows DoF to point to a remote value for faster/easier evaluation
    //  const int *m_remote_state;
    **/

    /// Allowed values, assume class T has method 'T.name()'
    std::vector<T> m_domain;
  };


  // int version
  jsonParser &to_json(OccupantDoF<int> const &dof, jsonParser &json);

  void from_json(OccupantDoF<int> &dof, jsonParser const &json);

  // molecule version
  jsonParser &to_json(OccupantDoF<Molecule> const &dof, jsonParser &json, Eigen::Matrix3d const &c2f_mat = Eigen::Matrix3d::Identity());

  void from_json(OccupantDoF<Molecule> &dof, jsonParser const &json, Eigen::Matrix3d const &f2c_mat = Eigen::Matrix3d::Identity());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class ContinuousDoF : public DoF {
  public:
    ContinuousDoF() : min_val(NAN), max_val(NAN), current_val(NAN),
      current_min(NAN), current_max(NAN), m_remote_ptr(nullptr) {}

    ContinuousDoF(DoFType::TypeFunc _type_func,
                  std::string const &_var_name,
                  double _min,
                  double _max) :
      DoF(_type_func, _var_name),
      min_val(_min),
      max_val(_max),
      current_val(NAN),
      current_min(min_val),
      current_max(max_val),
      m_remote_ptr(nullptr) {}

    ContinuousDoF(DoFType::TypeFunc _type_func,
                  std::string const &_var_name,
                  Index _ID,
                  double _min,
                  double _max) :
      DoF(_type_func, _var_name, _ID),
      min_val(_min),
      max_val(_max),
      current_val(NAN),
      current_min(min_val),
      current_max(max_val),
      m_remote_ptr(nullptr) {}

    double value() const {
      return current_val;
    }

    bool compare(ContinuousDoF const &RHS, bool compare_value) const {
      if(compare_value && !almost_equal(value(), RHS.value()))
        return false;

      return(type_name() == RHS.type_name()
             && var_name() == RHS.var_name()
             && (ID() == RHS.ID()));
    }

    bool operator==(ContinuousDoF const &RHS) const {
      return((type_name() == RHS.type_name())
             && var_name() == RHS.var_name()
             && (ID() == RHS.ID())
             && almost_equal(value(), RHS.value()));
    }

    double remote_value() const {
      assert(m_remote_ptr && "In ContinuousDoF::remote_value() m_remote_val is nullptr.\n");
      return *m_remote_ptr;
    }

    double const *remote_ptr() const {
      return m_remote_ptr;
    }

    void register_remote(RemoteHandle const &handle) {
      assert(handle.d_ptr() && "In ContinuousDoF::register_remote(), attempting to register to nullptr!");
      m_remote_ptr = handle.d_ptr();
    }

    ContinuousDoF *clone() const {
      return new ContinuousDoF(*this);
    }

    jsonParser &to_json(jsonParser &json) const {
      json.put_obj();
      json["DoF_type"] = "ContinuousDoF";
      json["type_name"] = type_name();
      json["var_name"] = var_name();
      json["min"] = min_val;
      json["max"] = max_val;
      json["current"] = current_val;
      return json;
    }

  private:
    double min_val, max_val, current_val;
    double current_min, current_max;

    /// Allows DoF to point to a remote value for faster/easier evaluation
    double const *m_remote_ptr;
  };

  void from_json(ContinuousDoF &dof, jsonParser const &json);

  jsonParser &to_json(ContinuousDoF const &dof, jsonParser &json);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DoFSet {
  public:
    using DoF_impl::Traits;

    typedef std::function<Traits()> TypeFunc;

    DoFSet(TypeFunc const &_type) :
      m_type_name(_type.type_name()) {}

    bool is_excluded_occ(std::string const &_occ_name) const {
      return m_excluded_occs.count(_occ_name);
    }

    /// \brief Matrix that relates DoFSet variables to a conventional coordiante system

    /// columns of coordinate_space() matrix are directions in conventional coordinate system
    /// so that  conventional_coord = DoFSet.coordinate_space()*DoFSet.values()
    /// coordinate_space() matrix has dimensions (N x size()), where N >= size()
    Eigen::MatrixXd const &coordinate_space() const {
      return m_coordinate_space;
    }

    /// \brief Return values of DoFs as a vector
    Eigen::VectorXd values() const;

    /// \brief Equivalent to coordinate_space()*values()
    Eigen::VectorXd conventional_values() const {
      return coordinate_space() * values();
    }

    void from_json(jsonParser const &json);

    jsonParser &to_json(jsonParser &json);


  private:
    std::string m_type_name;
    mutable SymGroupRepID m_sym_rep_ID;
    std::vector<ContinuousDoF> m_components;
    std::set<std::string> m_excluded_occs;
    Eigen::MatrixXd m_coordinate_space;
  };
}
#endif
