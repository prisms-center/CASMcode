#ifndef DOF_HH
#define DOF_HH

#include<vector>
#include<set>
#include<functional>
#include<string>
#include<boost/algorithm/string.hpp>
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/basis_set/DoFDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/symmetry/SymGroupRepID.hh"

namespace CASM {

  class AnisoValTraits;

  class SymGroup;
  class jsonParser;

  template<typename OccType>
  class OccupantDoF;

  class ContinuousDoF;
  class DoFSet;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  namespace DoF {
    using BasicTraits = AnisoValTraits;

    /// A RemoteHandle can be initialized with either a double or integer reference and then passed to a
    /// DoF (or to a BasisFunction, which contains DoFs), in order for the DoF to access the remote value

    // For now this isn't used, but may be in the future if we need a dynamic clexulator class for prototyping
    class RemoteHandle {
    public:
      // \brief Construct with basic identification info
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
    class Base {
    public:
      using RemoteHandle = DoF::RemoteHandle;

      Base() :
        m_traits(BasicTraits::null()),
        m_var_name("NULL"),
        m_dof_ID(-1),
        m_ID_lock(false) {}

      Base(BasicTraits const &traits,
           std::string const &_var_name,
           Index _ID);

      BasicTraits const &traits() const {
        return m_traits;
      }

      /// \brief Const access of DoF type name
      std::string type_name() const {
        return m_traits.name();
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

    private:

      BasicTraits m_traits;
      std::string m_var_name;

      /// dof_ID is a way to distinguish between DoFs with the same name but different identities
      /// dof_ID for now usually refers to the site index of a cluster (e.g., 0, 1, 2 of a triplet)
      /// or an index into the primitive cell neighbor list. Other usage cases may be introduced later
      Index m_dof_ID;

      /// Locks the ID so that it can't be updated.  Is used for global DoF's
      bool m_ID_lock;
    };
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DiscreteDoF : public DoF::Base {
  public:
    using Base = DoF::Base;
    using BasicTraits = DoF::BasicTraits;
    DiscreteDoF():
      Base(),
      m_current_state(0),
      m_remote_state(nullptr),
      m_symrep_ID(SymGroupRepID::identity(0)) {}


    DiscreteDoF(BasicTraits const &_traits,
                std::string const &_var_name,
                Index _dof_ID = -1,
                int _current_state = 0,
                SymGroupRepID _id = SymGroupRepID::identity(0)) :
      Base(_traits, _var_name, _dof_ID),
      m_current_state(_current_state),
      m_remote_state(nullptr),
      m_symrep_ID(_id) {

    }

    virtual
    ~DiscreteDoF() {}

    SymGroupRepID symrep_ID()const {
      return m_symrep_ID;
    }

    void set_symrep_ID(SymGroupRepID _id) const {
      m_symrep_ID = _id;
    }

    /// \brief Allocates an empty symmetry representation in \param _group and records its SymGroupRepID
    /// This representation becomes THE representation for this DiscreteDoF, to be initialized accordingly
    void allocate_symrep(SymGroup const &_group) const;

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
    mutable SymGroupRepID m_symrep_ID;

  private:
    virtual DiscreteDoF *_clone() const = 0;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template< typename T>
  class OccupantDoF : public DiscreteDoF {
  public:
    OccupantDoF(BasicTraits const &_traits) :
      DiscreteDoF(_traits, "s") { }

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


    const T &occ() const {
      return m_domain[m_current_state];
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

    void set_domain(std::vector<T> new_dom) {
      m_domain = std::move(new_dom);
      if(m_domain.size() == 1)
        m_current_state = 0;
      else
        m_current_state = -1;
      set_symrep_ID(SymGroupRepID::identity(m_domain.size()));
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

    std::unique_ptr<OccupantDoF> clone() const {
      return std::unique_ptr<OccupantDoF>(static_cast<OccupantDoF *>(this->_clone()));
    }

  private:

    DiscreteDoF *_clone() const override {
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

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class ContinuousDoF : public DoF::Base {
  public:
    using Base = DoF::Base;
    using BasicTraits = DoF::BasicTraits;

    ContinuousDoF(BasicTraits const &_traits,
                  std::string const &_var_name,
                  Index _ID,
                  double _min,
                  double _max) :
      Base(_traits, _var_name, _ID),
      min_val(_min),
      max_val(_max),
      current_val(NAN),
      current_min(min_val),
      current_max(max_val),
      m_remote_ptr(nullptr) {}

    ContinuousDoF(BasicTraits const &_traits)
      : ContinuousDoF(_traits, "", -1, NAN, NAN) {}

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
      return notstd::make_unique<ContinuousDoF>(*this);
    }

  private:
    double min_val, max_val, current_val;
    double current_min, current_max;

    /// Allows DoF to point to a remote value for faster/easier evaluation
    mutable double const *m_remote_ptr;
  };

  inline
  bool compare_no_value(ContinuousDoF const &A, ContinuousDoF const &B) {
    return A.compare(B);
  }

}
#endif
