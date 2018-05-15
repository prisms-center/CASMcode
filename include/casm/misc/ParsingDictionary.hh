#ifndef CASM_ParsingDictionary
#define CASM_ParsingDictionary

namespace CASM {
  namespace ParsingDictionary_impl {

    /// \brief Conversion Functor for inserting BasicTraits into unique_cloneable_map
    template<typename T>
    struct DictionaryConverter {
      // performs a static_cast of value.clone().unique().release()
      notstd::cloneable_ptr<T> operator()(const T &_t) {
        return notstd::cloneable_ptr<T>(static_cast<T *>(_t.clone().release()));
      }
    };
  }

  /// \brief  Parsing dictionary for obtaining the correct MoleculeAttribute given a name
  template<typename T>
  class ParsingDictionary :
    public notstd::unique_cloneable_map<std::string, T> {

  public:
    typedef notstd::unique_cloneable_map<std::string, T> Base;

    typedef typename Base::key_type key_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::size_type size_type;
    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

    using Base::find;
    using Base::begin;
    using Base::end;

    ParsingDictionary() :
      Base([](const value_type & value)->std::string {
      return value.name();
    },
    ParsingDictionary_impl::DictionaryConverter<value_type>()) {}

    using Base::insert;

    /// \brief Equivalent to find, but set 'home' and throws error with
    /// suggestion if @param _name not found
    std::unique_ptr<T> lookup(const key_type &_name) const;

    /// \brief True if dictionary contains entry for @param _name
    bool contains(const key_type &_name) const {
      return this->find(_name) != this->end();
    }

    void print_help(std::ostream &_stream,
                    int width = 60,
                    int separation = 8) const;

  };


  /// \brief Equivalent to find, but throw error with suggestion if _name not found
  template<typename T>
  std::unique_ptr<T> ParsingDictionary<T>::lookup(const ParsingDictionary<T>::key_type &_name) const {

    auto res = find(_name);
    if(res != end()) {
      return res->clone();
    }
    else {

      // If no match, try to use demerescau-levenshtein distance to make a helpful suggestion
      auto it = begin();
      int min_dist(-1);
      for(; it != end(); ++it) {
        int dist = dl_string_dist(_name, it->name());
        if(min_dist < 0 || dist < min_dist) {
          min_dist = dist;

          res = it;
        }
      }

      throw std::runtime_error("CRITICAL ERROR: Invalid " + T::class_desc() + " \"" + _name + "\" specified.\n"
                               + "                Did you mean \"" + res->name() + "\"?\n");

    }

  }
}
#endif
