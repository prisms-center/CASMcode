#ifndef DOFTRAITS_HH
#define DOFTRAITS_HH
namespace CASM {
  class jsonParser;
  class MasterSymGroup;

  namespace DoF_impl {
    enum DOF_DOMAIN {DISCRETE, CONTINUOUS};
    enum DOF_MODE {LOCAL, GLOBAL};

    /// \brief Collection of all the traits specific to a DoF type

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
      virtual bool time_reversal_active() const = 0;

      /// \brief returns true if DoF tracks a BasicTraits (specified by @param attr_name)
      virtual bool obscures_molecule_attribute(std::string const &attr_name) const = 0;

      /// \brief returns true if DoF tracks the orientation of the occupying molecule (not typical)
      virtual bool obscures_occupant_orientation() const = 0;

      /// \brief returns true if DoF tracks the chirality of the occupying molecule (not typical)
      virtual bool obscures_occupant_chirality() const = 0;

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      virtual BasisSet construct_site_basis(Site const &_site) const = 0;

      /// \brief Populate @param _in from JSON
      virtual void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const = 0;

      /// \brief Output @param _in to JSON
      virtual void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const = 0;

      /// \brief Generate a symmetry representation for this DoF
      virtual SymGroupRepID generate_symrep(MasterSymGroup const &_group) const = 0;

      /// \brief non-virtual method to obtain copy through BasicTraits pointer
      std::unique_ptr<BasicTraits> clone() const {
        return _clone();
      }
    private:
      virtual BasicTraits *_clone() const = 0;

      std::string m_type_name;
      DOF_DOMAIN m_domain;
      DOF_MODE m_mode;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief Conversion Functor for inserting BasicTraits into unique_cloneable_map
    struct DictionaryConverter {
      // performs a static_cast of value.clone().unique().release()
      notstd::cloneable_ptr<BasicTraits> operator()(const BasicTraits &_traits) {
        return notstd::cloneable_ptr<BasicTraits>(static_cast<BasicTraits *>(_traits.clone().release()));
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief  Parsing dictionary for obtaining the correct BasicTraits given a name
    class TraitsDictionary :
      public notstd::unique_cloneable_map<std::string, BasicTraits> {

    public:
      typedef notstd::unique_cloneable_map<std::string, BasicTraits> UniqueMapType;
      typedef typename UniqueMapType::key_type key_type;
      typedef typename UniqueMapType::value_type value_type;
      typedef typename UniqueMapType::size_type size_type;
      typedef typename UniqueMapType::iterator iterator;
      typedef typename UniqueMapType::const_iterator const_iterator;

      TraitsDictionary() :
        UniqueMapType([](const value_type & value)->std::string {
        return value.name();
      },
      DictionaryConverter()) {}

      using UniqueMapType::insert;

      /// \brief Equivalent to find, but throw error with suggestion if @param _name not found
      notstd::cloneable_ptr<BasicTraits> lookup(const key_type &_name) const;

      /// \brief True if dictionary contains entry for @param _name
      bool contains(const key_type &_name) const {
        return this->find(_name) != this->end();
      }

      //void print_help(std::ostream &_stream) const;
    };

    /// This will eventually be managed by ProjectSettings
    TraitsDictionary const &traits_dict();
  }

}
#endif
