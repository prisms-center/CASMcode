#ifndef CASM_DiffTransConfigDatabase
#define CASM_DiffTransConfigDatabase

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace boost {

  template<typename T>
  class iterator_range;

}

namespace CASM {

  namespace DB {

    /// Derived DiffTransConfigDatabase must implement public methods:
    /// - iterator update(const DiffTransConfiguration &config)
    /// - boost::iterator_range<iterator> orbit_scel_range(const name_type& difftransname,const name_type& scelname) const
    /// - boost::iterator_range<iterator> orbit_range(const name_type& difftransname) const
    ///
    /// Database insert methods by convention do not enforce canonical forms.
    /// That logic is included in DiffTransConfiguration::insert, which is
    /// the safest way to insert new DiffTransConfiguration in the database. But in cases
    /// where it is known that a DiffTransConfiguration is generated in
    /// canonical form, the Database insert methods may be used directly, for
    /// instance as in the method insert_unique_canon_diff_trans_configs.
    ///
    template<>
    class Database<DiffTransConfiguration> : public ValDatabase<DiffTransConfiguration> {

    public:

      Database(const PrimClex &_primclex) :
        ValDatabase<DiffTransConfiguration>(_primclex) {}

      virtual ~Database() {}


      /// Update record
      virtual iterator update(const DiffTransConfiguration &diff_trans_config) = 0;

      /// Range of DiffTransConfiguration in a particular supercell within an orbit
      ///
      /// - Should return range {end(), end()} if no DiffTransConfiguration in specified Supercell within an orbit
      /// - Note: boost::iterator_range<iterator>::size is not valid for
      ///   DatabaseIterator.  Use boost::distance instead.
      virtual boost::iterator_range<iterator> orbit_scel_range(const std::string &diff_trans_name, const std::string &scelname) const = 0;

      /// Number of DiffTransConfiguration in a particular supercell within an orbit
      Index orbit_scel_range_size(const std::string &diff_trans_name, const std::string &scelname) const;

      /// Range of DiffTransConfiguration in a particular diff_trans
      ///
      /// - Should return range {end(), end()} if no DiffTransConfiguration in specified diff_trans_name
      /// - Note: boost::iterator_range<iterator>::size is not valid for
      ///   DatabaseIterator.  Use boost::distance instead.
      virtual boost::iterator_range<iterator> orbit_range(const std::string &diff_trans_name) const = 0;

      /// Number of DiffTransConfiguration in a particular diff_trans
      Index orbit_range_size(const std::string &diff_trans_name) const;

    };

  }
}

#endif
