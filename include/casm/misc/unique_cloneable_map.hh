#ifndef CASM_unique_cloneable_map_HH
#define CASM_unique_cloneable_map_HH

#include "casm/misc/unique_map.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace notstd {

  template<typename MapType>
  struct GetSecondDereferenced {

    typedef typename MapType::iterator::reference reference;
    typedef typename MapType::iterator::value_type::second_type::reference result_type;

    GetSecondDereferenced() {}

    result_type operator()(reference pair) const {
      return *(pair.second);
    }

  };

  template<typename MapType>
  struct GetSecondDereferencedConst {

    typedef typename MapType::const_iterator::reference reference;
    typedef const typename MapType::const_iterator::value_type::second_type::reference result_type;

    GetSecondDereferencedConst() {}

    result_type operator()(reference pair) const {
      return *(pair.second);
    }

  };



  /// \brief Enable use of unique_map of cloneable_ptr<ValueType> via references
  ///        to ValueType
  template<typename KeyType, typename ValueType, typename MapType = std::map<KeyType, cloneable_ptr<ValueType> > >
  using unique_cloneable_map = unique_map<KeyType,
        ValueType,
        MapType,
        GetSecondDereferenced<MapType>,
        GetSecondDereferencedConst<MapType> >;


}

#endif

