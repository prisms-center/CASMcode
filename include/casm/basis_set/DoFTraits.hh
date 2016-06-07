#ifndef DOFTRAITS_HH
#define DOFTRAITS_HH
namespace CASM {
  class jsonParser;
  class MasterSymGroup;

  namespace DoF_impl {
    /// \brief Collection of all the traits specific to a DoF type

    class Traits : public BasicTraits {
    public:
      Traits(std::string const &_type_name,
             DOF_DOMAIN _domain,
             DOF_MODE _mode) :
        BasicTraits(m_type_name(_type_name),
                    m_domain(_domain),
                    m_mode(_mode)) {

      }

      /// \brief Allow destruction through base pointer
      virtual ~Traits() {}

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      virtual BasisSet construct_site_basis(Site const &_site) const = 0;

      /// \brief Populate @param _in from JSON
      virtual void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const = 0;

      /// \brief Output @param _in to JSON
      virtual void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const = 0;

      /// \brief Generate a symmetry representation for this DoF
      virtual SymGroupRepID generate_symrep(MasterSymGroup const &_group) const = 0;

      /// \brief non-virtual method to obtain copy through Traits pointer
      std::unique_ptr<Traits> clone() const {
        return static_cast<Traits *>(_clone());
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class OccupationDoFTraits : public Traits {
    public:
      OccupationDoFTraits():
        Traits("occupation",
               DISCRETE,
               LOCAL) {
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      BasisSet construct_site_basis(Site const &_site) const override {
        throw std::runtime_error("OccupationDoFTraits::construct_site_basis not implemented!");
      }

      /// \brief Populate @param _in from JSON
      void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const override {
        throw std::runtime_error("OccupationDoFTraits::from_json not implemented!");
      }

      /// \brief Output @param _in to JSON
      void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const override {
        throw std::runtime_error("OccupationDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for this DoF
      SymGroupRepID generate_symrep(MasterSymGroup const &_group) const override {
        throw std::runtime_error("OccupationDoFTraits::generate_symrep not implemented!");
      }

    protected:
      BasicTraits *_clone() const override {
        return OccupationDoFTraits(*this);
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    namespace DoFType {
      DoF_impl::OccupationDoFTraits occupation() {
        return DoF_impl::OccupationDoFTraits();
      }
    }

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
