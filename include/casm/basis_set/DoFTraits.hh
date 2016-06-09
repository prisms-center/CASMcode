#ifndef DOFTRAITS_HH
#define DOFTRAITS_HH

#include "casm/basis_set/DoF.hh"

// remove these once implementation of derived classes gets moved out of this file
#include "casm/crystallography/Site.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {
  class jsonParser;
  class MasterSymGroup;
  class Structure;

  namespace DoF_impl {
    class Traits;
  }

  namespace DoFType {
    DoF_impl::Traits const &traits(std::string const &dof_key);

    DoF_impl::BasicTraits const &basic_traits(std::string const &dof_key);

    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> occupation();
  }

  // This is a hack to forward-declare IntegralCluster.  Forward declarations
  // of typedefs should probably get their own *.hh files, without any dependencies
  template<typename CoordType>
  class CoordCluster;
  class UnitCellCoord;
  typedef CoordCluster<UnitCellCoord> IntegralCluster;

  template <typename ClustType, typename SymCompareType>
  class Orbit;

  template <typename ClustType>
  class PrimPeriodicSymCompare;

  namespace DoF_impl {
    /// \brief Collection of all the traits specific to a DoF type

    class Traits : public BasicTraits {
    public:
      Traits(std::string const &_type_name,
             DOF_DOMAIN _domain,
             DOF_MODE _mode) :
        BasicTraits(_type_name,
                    _domain,
                    _mode) {

      }

      /// \brief Allow destruction through base pointer
      virtual ~Traits() {}

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      virtual std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                         std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                         jsonParser const &_bspecs) const = 0;

      /// \brief Populate @param _in from JSON
      virtual void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const = 0;

      /// \brief Output @param _in to JSON
      virtual void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const = 0;

      /// \brief Generate a symmetry representation for this DoF
      virtual SymGroupRepID generate_symrep(MasterSymGroup const &_group) const = 0;

      virtual std::string site_basis_description(BasisSet site_bset, Site site) const = 0;

      /// \brief non-virtual method to obtain copy through Traits pointer
      std::unique_ptr<Traits> clone() const {
        return std::unique_ptr<Traits>(static_cast<Traits *>(_clone()));
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

      std::string site_basis_description(BasisSet site_bset, Site site) const override {
        std::stringstream ss;
        if(site_bset.size() == 0)
          ss << "        [No site basis functions]\n\n";
        std::vector<DoF::RemoteHandle> handles(1, site.site_occupant().handle());
        int s;
        handles[0] = s;
        site_bset.register_remotes(handles);
        for(Index f = 0; f < site_bset.size(); f++) {
          for(s = 0; s < site.site_occupant().size(); s++) {
            if(s == 0)
              ss << "    ";
            ss << "    \\phi_" << site.basis_ind() << '_' << f << '[' << site.site_occupant()[s].name() << "] = "
               << site_bset[f]->remote_eval();
            if(s + 1 == site.site_occupant().size())
              ss << "\n";
            else
              ss << ",   ";
          }
        }
        return ss.str();
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;
    protected:
      BasicTraits *_clone() const override {
        return new OccupationDoFTraits(*this);
      }
    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief Conversion Functor for inserting BasicTraits into unique_cloneable_map
    struct DictionaryConverter {
      // performs a static_cast of value.clone().unique().release()
      notstd::cloneable_ptr<BasicTraits> operator()(const BasicTraits &_traits) {
        return notstd::cloneable_ptr<BasicTraits>(_traits.clone().release());
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
        return value.type_name();
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


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  namespace DoFType {
    DoF_impl::Traits const &traits(std::string const &dof_key) {
      return static_cast<DoF_impl::Traits const &>(DoF::traits(dof_key));
    }

    DoF_impl::BasicTraits const &basic_traits(std::string const &dof_key) {
      return DoF::traits(dof_key);
    }
    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> occupation() {
      return DoF_impl::OccupationDoFTraits().clone();
    }
  }

}
#endif
