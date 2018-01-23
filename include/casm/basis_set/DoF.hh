#ifndef DOF_HH
#define DOF_HH
#include<vector>
#include<set>
#include<functional>
#include<boost/algorithm/string.hpp>
#include "casm/CASM_global_definitions.hh"
#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {
  class Molecule;
  class jsonParser;

  template<typename OccType>
  class OccupantDoF;
  using MoleculeOccupant = OccupantDoF<Molecule>;

  class ContinuousDoF;
  class DoFSet;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  namespace DoF_impl {
    enum DOF_DOMAIN {DISCRETE, CONTINUOUS};
    enum DOF_MODE {LOCAL, GLOBAL};

    /// \brief Base class for defining a collection of traits shared by a specific DoF type
    /// The BasicTraits class maintains only the traits that do not depend on other CASM classes

    /// In future, may include function pointers (wrapped in std::function<>) for controlling certain parts
    /// of program execution
    class BasicTraits {
    public:
      BasicTraits(std::string const &_type_name,
                  DOF_DOMAIN _domain,
                  DOF_MODE _mode) :
        m_type_name(_type_name),
        m_domain(_domain),
        m_mode(_mode) {
      }

      /// \brief Allow destruction through base pointer
      virtual ~BasicTraits() {}

      /// \brief const access of type_name
      std::string const &type_name() const {
        return m_type_name;
      }

      /// \brief returns true if DoF is global
      bool global()const {
        return m_mode == GLOBAL;
      }

      /// \brief returns true if DoF is discrete
      bool discrete() const {
        return m_domain == DISCRETE;
      }

      /// \brief equality comparison of type_name
      bool operator==(std::string const &other_name) const {
        return type_name() == other_name;
      }

      /// \brief lexicographic comparison of type_name
      bool operator<(std::string const &other_name) const {
        return type_name() < other_name;
      }

      /// \brief comparison of type_name, domain (discrete/continuous) and mode (local/global)
      bool identical(BasicTraits const &other) const {
        return type_name() == other.type_name()
               && m_domain == other.m_domain
               && m_mode == other.m_mode;
      }

      /// \brief allow implicit conversion to std::string (type_name)
      operator std::string const &() const {
        return type_name();
      }

      /// \brief returns true if time-reversal changes the DoF value
      virtual bool time_reversal_active() const {
        return false;
      }

      /// \brief returns true if DoF tracks a BasicTraits (specified by @param attr_name)
      virtual bool obscures_molecule_attribute(std::string const &attr_name) const {
        return false;
      }

      /// \brief returns true if DoF tracks the orientation of the occupying molecule (not typical)
      virtual bool obscures_occupant_orientation() const {
        return false;
      }

      /// \brief returns true if DoF tracks the chirality of the occupying molecule (not typical)
      virtual bool obscures_occupant_chirality() const {
        return false;
      }

      /// \brief implements json parsing of a specialized DoF.
      // In future, we may need to add another inheritance layer to handle DiscreteDoF types
      virtual void from_json(ContinuousDoF &_dof, jsonParser const &json) const;


      /// \brief implements json parsing of a specialized DoF.
      // In future, we may need to add another inheritance layer to handle DiscreteDoF types
      virtual void to_json(ContinuousDoF const &_dof, jsonParser &json) const;

      /// \brief implements json parsing of a specialized DoFSet.
      virtual void from_json(DoFSet &_dof, jsonParser const &json) const;


      /// \brief implements json parsing of a specialized DoFSet.
      virtual void to_json(DoFSet const &_dof, jsonParser &json) const;

      std::unique_ptr<BasicTraits> clone() const {
        return std::unique_ptr<BasicTraits>(_clone());
      }

    protected:
      virtual BasicTraits *_clone() const = 0;

      std::string m_type_name;
      DOF_DOMAIN m_domain;
      DOF_MODE m_mode;
    };


    //struct TraitsConverter {
    inline
    notstd::cloneable_ptr<BasicTraits> traits2cloneable_ptr(const BasicTraits &value) {
      return notstd::cloneable_ptr<BasicTraits>(value.clone().release());
    }
    //};


    /// \brief A class to manage dynamic evaluation of BasisFunctions

    /// A RemoteHandle can be initialized with either a double or integer reference and then passed to a
    /// DoF (or to a BasisFunction, which contains DoFs), in order for the DoF to access the remote value

    // For now this isn't used, but may be in the future if we need a dynamic clexulator class for prototyping
    class RemoteHandle {
    public:
      // \brief Initialize pointer to double
      RemoteHandle(std::string const &_type_name,
                   std::string const &_var_name,
                   Index _dof_ID) :
        m_type_name(_type_name),
        m_var_name(_var_name),
        m_dof_ID(_dof_ID) {

      }

      // \brief Initialize pointer to double
      RemoteHandle &operator=(double const &d) {
        m_d_ptr = &d;
        m_i_ptr = nullptr;
        return *this;
      }

      // \brief Initialize pointer to int
      RemoteHandle &operator=(int const &i) {
        m_d_ptr = nullptr;
        m_i_ptr = &i;
        return *this;
      }

      // \brief prevent construction by non-addressable integer
      RemoteHandle &operator=(double &&d) = delete;

      // \brief prevent construction by non-addressable double
      RemoteHandle &operator=(int &&i) = delete;

      // \brief Less-than comparison to allow use of STL algorithms
      bool operator<(RemoteHandle const &_rhs)const {
        return m_dof_ID < _rhs.m_dof_ID
               || (m_dof_ID == _rhs.m_dof_ID
                   && (m_type_name < _rhs.m_type_name
                       || (m_type_name == _rhs.m_type_name
                           && m_var_name < _rhs.m_var_name)));
      }

      // \brief Equality comparison, irrespective of internal pointers (compares name and ID only)
      bool operator==(RemoteHandle const &_rhs)const {
        return m_dof_ID == _rhs.m_dof_ID
               && m_type_name == _rhs.m_type_name
               && m_var_name == _rhs.m_var_name;
      }

      // \brief const access of integer pointer
      const int *i_ptr() const {
        return m_i_ptr;
      }

      // \brief const access of double pointer
      const double *d_ptr() const {
        return m_d_ptr;
      }

    private:
      double const *m_d_ptr;
      int const *m_i_ptr;

      std::string m_type_name;
      std::string m_var_name;
      Index m_dof_ID;
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
    using BasicTraits = DoF_impl::BasicTraits;
    using RemoteHandle = DoF_impl::RemoteHandle;

    /// DoFs are initialized using a TypeFunc
    /// e.g., DerivedDoF my_dof(DoFType::occupation);
    typedef std::function<notstd::cloneable_ptr<BasicTraits>()> TypeFunc;

    static BasicTraits const &traits(std::string const &_type_name);

    DoF() :
      m_type_name("EMPTY"),
      m_var_name("EMPTY"),
      m_dof_ID(-1),
      m_ID_lock(false) {}

    DoF(BasicTraits const &traits,
        std::string const &_var_name,
        Index _ID);

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

    /// \brief Create a RemoteHandle that refers to this DoF
    RemoteHandle handle()const {
      return RemoteHandle(type_name(), var_name(), ID());
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

    std::unique_ptr<DoF> clone() const {
      return std::unique_ptr<DoF>(this->_clone());
    }

  protected:
    void _set_type_name(std::string _type_name) {
      std::swap(m_type_name, _type_name);
    }

    void _set_var_name(std::string _var_name) {
      std::swap(m_var_name, _var_name);
    }
  private:

    virtual DoF *_clone() const = 0;

    typedef notstd::unique_cloneable_map<std::string, BasicTraits> TraitsMap;
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
    DiscreteDoF():
      DoF(),
      m_current_state(0),
      m_remote_state(nullptr),
      m_sym_rep_ID(SymGroupRepID::identity(0)) {}


    DiscreteDoF(BasicTraits const &_traits,
                std::string const &_var_name,
                Index _dof_ID = -1,
                int _current_state = 0,
                SymGroupRepID _id = SymGroupRepID::identity(0)) :
      DoF(_traits, _var_name, _dof_ID),
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

    void register_remote(RemoteHandle const &handle) const {
      assert(handle.i_ptr() && "In DiscreteDoF::register_remote(), attempting to register to nullptr!");
      m_remote_state = handle.i_ptr();
    }

    std::unique_ptr<DiscreteDoF> clone() const {
      return std::unique_ptr<DiscreteDoF>(this->_clone());
    }

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

  protected:

    /// index into domain of the current state, -1 if unspecified
    int m_current_state;

    /// Allows DoF to point to a remote value for faster/easier evaluation
    mutable int const *m_remote_state;

    /// ID for the permutation representation for occupants
    mutable SymGroupRepID m_sym_rep_ID;

  private:
    virtual DiscreteDoF *_clone() const = 0;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template< typename T>
  class OccupantDoF : public DiscreteDoF {
  public:
    OccupantDoF(BasicTraits const &_traits) :
      DiscreteDoF(_traits, "s") { }

    OccupantDoF(TypeFunc _func) :
      OccupantDoF(_func, "s") {}

    OccupantDoF(BasicTraits const &_traits,
                std::string const &_var_name,
                std::vector<T> const &_domain,
                int _current_state = 0) :
      DiscreteDoF(_traits,
                  _var_name,
                  -1,
                  _current_state,
                  SymGroupRepID::identity(_domain.size())),
      m_domain(_domain) { }

    OccupantDoF(TypeFunc _func,
                std::string const &_var_name,
                std::vector<T> const &_domain,
                int _current_state = 0) :
      OccupantDoF(*_func(),
                  _var_name,
                  _domain,
                  _current_state) {}

    OccupantDoF(BasicTraits const &_traits,
                std::vector<T> const &_domain,
                int _current_state = 0) :
      DiscreteDoF(_traits,
                  "s",
                  -1,
                  _current_state,
                  SymGroupRepID::identity(_domain.size())),
      m_domain(_domain) { }

    OccupantDoF(TypeFunc _func,
                std::vector<T> const &_domain,
                int _current_state = 0) :
      OccupantDoF(*_func(),
                  _domain,
                  _current_state) {}

    /*
    OccupantDoF(TypeFunc _func,
                std::initializer_list<T> _domain,
                int _current_state = 0) :
      DiscreteDoF(*_func(),
                  "s",
                  -1,
                  _current_state,
                  SymGroupRepID::identity(_domain.size())),
      m_domain(_domain.begin(),
               _domain.end()) {    }
    */
    const T &occ() const {
      return m_domain[m_current_state];
    }

    std::unique_ptr<OccupantDoF<T> > clone() const {
      return std::unique_ptr<OccupantDoF<T> >(this->_clone());
    }

    /// \brief Set occupant by index
    void set_value(int occ_index) override {
      if(Index(occ_index) >= m_domain.size()) {
        throw std::runtime_error("In OccupantDoF::set_value(): Bad Assignment, occ_index>=size.\n");
      }
      m_current_state = occ_index;
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
        if(!((*this)[i] == RHS[i]))
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

    template<typename...Args>
    void from_json(jsonParser const &json, Args &&... args);

  private:

    virtual OccupantDoF *_clone() const override {
      return new OccupantDoF(*this);
    }
    /**** Inherited from DiscreteDoF ****
    // index into domain of the current state, -1 if unspecified
    //  int m_current_state;

    // Allows DoF to point to a remote value for faster/easier evaluation
    //  const int *m_remote_state;
    **/

    /// Allowed values, assume class T has method 'T.name()()'
    std::vector<T> m_domain;
  };

  template<typename OccType, typename...Args>
  void from_json(OccupantDoF<OccType> &dof, const jsonParser &json, Args &&... args);

  template<typename OccType, typename...Args>
  jsonParser &to_json(OccupantDoF<int> const &dof, jsonParser &json, Args &&... args);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class ContinuousDoF : public DoF {
  public:
    ContinuousDoF(BasicTraits const &_traits,
                  std::string const &_var_name,
                  Index _ID,
                  double _min,
                  double _max) :
      DoF(_traits, _var_name, _ID),
      min_val(_min),
      max_val(_max),
      current_val(NAN),
      current_min(min_val),
      current_max(max_val),
      m_remote_ptr(nullptr) {}

    ContinuousDoF(BasicTraits const &_traits)
      : ContinuousDoF(_traits, "", -1, NAN, NAN) {}

    ContinuousDoF(TypeFunc _func,
                  std::string const &_var_name,
                  Index _ID,
                  double _min,
                  double _max) :
      ContinuousDoF(*_func(), _var_name, _ID, _min, _max) {}

    double value() const {
      return current_val;
    }

    bool compare(ContinuousDoF const &RHS, bool compare_value = false) const {
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

    void register_remote(RemoteHandle const &handle) const {
      assert(handle.d_ptr() && "In ContinuousDoF::register_remote(), attempting to register to nullptr!");
      m_remote_ptr = handle.d_ptr();
    }

    std::unique_ptr<ContinuousDoF> clone() const {
      return std::unique_ptr<ContinuousDoF>(this->_clone());
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
    virtual ContinuousDoF *_clone() const override {
      return new ContinuousDoF(*this);
    }

    double min_val, max_val, current_val;
    double current_min, current_max;

    /// Allows DoF to point to a remote value for faster/easier evaluation
    mutable double const *m_remote_ptr;
  };

  inline
  bool compare_no_value(ContinuousDoF const &A, ContinuousDoF const &B) {
    return A.compare(B);
  }

  void from_json(ContinuousDoF &dof, jsonParser const &json);

  jsonParser &to_json(ContinuousDoF const &dof, jsonParser &json);

  template<>
  struct jsonConstructor<ContinuousDoF> {

    /// \brief Allows ContinuousDoF to be constructed properly from json input
    static ContinuousDoF from_json(const jsonParser &json, DoF::BasicTraits const &_traits);
  };


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DoFSet {
  public:
    using BasicTraits = DoF_impl::BasicTraits;

    using TypeFunc =  std::function<notstd::cloneable_ptr<BasicTraits>()>;

    using Container = std::vector<ContinuousDoF>;

    using const_iterator = std::vector<ContinuousDoF>::const_iterator;

    DoFSet(TypeFunc const &_type) :
      m_type_name(_type()->type_name()) {}

    Index size() const {
      return m_components.size();
    }

    std::string const &type_name() const {
      return m_type_name;
    }

    ContinuousDoF const &operator[](Index i) const {
      return m_components[i];
    }

    const_iterator begin() const {
      return m_components.cbegin();
    }

    const_iterator end() const {
      return m_components.cend();
    }

    const_iterator cbegin() const {
      return m_components.cbegin();
    }

    const_iterator cend() const {
      return m_components.cend();
    }

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

    SymGroupRepID const &sym_rep_ID() const {
      return m_sym_rep_ID;
    }

    /// \brief Return values of DoFs as a vector
    Eigen::VectorXd values() const;

    /// \brief Equivalent to coordinate_space()*values()
    Eigen::VectorXd conventional_values() const {
      return coordinate_space() * values();
    }

    bool identical(DoFSet const &rhs) const;

    bool update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs);

    void from_json(jsonParser const &json);

    jsonParser &to_json(jsonParser &json) const;


  private:
    std::string m_type_name;
    mutable SymGroupRepID m_sym_rep_ID;
    std::vector<ContinuousDoF> m_components;
    std::set<std::string> m_excluded_occs;
    Eigen::MatrixXd m_coordinate_space;
  };

  //********************************************************************

  template<typename OccType> template<typename...Args>
  void OccupantDoF<OccType>::from_json(const jsonParser &json, Args &&... args) {
    _set_type_name(json["type_name"].get<std::string>());
    _set_var_name(json["var_name"].get<std::string>());
    set_ID(json["ID"].get<Index>());
    m_domain.clear();

    CASM::from_json(m_domain, json["domain"], std::forward<Args>(args)...);

    json.get_else(m_current_state, "value", int(-1));

  }

  //********************************************************************

  inline
  bool operator==(DoFSet const &A, DoFSet const &B) {
    return A.identical(B);
  }

  //********************************************************************

  inline
  bool operator!=(DoFSet const &A, DoFSet const &B) {
    return !A.identical(B);
  }

  //********************************************************************

  template<typename OccType, typename...Args>
  void from_json(OccupantDoF<OccType> &_dof, const jsonParser &json, Args &&... args) {
    _dof.from_json(json, std::forward<Args>(args)...);
  }

  //********************************************************************
  // molecule version
  template<typename OccType, typename...Args>
  jsonParser &to_json(OccupantDoF<OccType> const &_dof, jsonParser &json, Args &&... args) {
    json.put_obj();
    json["type_name"] = _dof.type_name();
    json["var_name"] = _dof.var_name();
    json["ID"] = _dof.ID();

    to_json(_dof.domain(), json["domain"], std::forward<Args>(args)...);
    json["value"] = _dof.value();

    return json;
  }

}
#endif
