#ifndef CASM_ParsingDictionary
#define CASM_ParsingDictionary

#include "casm/misc/CASM_math.hh"
#include "casm/misc/unique_cloneable_map.hh"
namespace CASM {
namespace ParsingDictionary_impl {

/// \brief Conversion Functor for inserting BasicTraits into
/// unique_cloneable_map
template <typename T>
struct DictionaryConverter {
  // performs a static_cast of value.clone().unique().release()
  notstd::cloneable_ptr<T> operator()(const T &_t) {
    return notstd::cloneable_ptr<T>(static_cast<T *>(_t.clone().release()));
  }
};
}  // namespace ParsingDictionary_impl

/// \brief  Parsing dictionary for obtaining the correct MoleculeAttribute given
/// a name
template <typename T>
class ParsingDictionary : public notstd::unique_cloneable_map<std::string, T> {
 public:
  typedef notstd::unique_cloneable_map<std::string, T> Base;

  typedef typename Base::key_type key_type;
  typedef typename Base::value_type value_type;
  typedef typename Base::size_type size_type;
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;

  using Base::begin;
  using Base::end;
  using Base::find;
  using Base::insert;

  ParsingDictionary()
      : Base(
            [](const value_type &value) -> std::string { return value.name(); },
            ParsingDictionary_impl::DictionaryConverter<value_type>()) {}

  std::unique_ptr<ParsingDictionary<value_type> > clone() const {
    return notstd::make_unique<ParsingDictionary<value_type> >(*this);
  }

  /// \brief Equivalent to find, but set 'home' and throws error with
  /// suggestion if @param _name not found
  T const &lookup(const key_type &_name) const;

  /// \brief True if dictionary contains entry for @param _name
  bool contains(const key_type &_name) const {
    return this->find(_name) != this->end();
  }

  void print_help(std::ostream &_stream, std::function<bool(T const &)> filter,
                  int width = 60, int separation = 8) const;
};

template <typename T>
ParsingDictionary<T> make_parsing_dictionary();

/// \brief Equivalent to find, but throw error with suggestion if _name not
/// found
template <typename T>
T const &ParsingDictionary<T>::lookup(
    const ParsingDictionary<T>::key_type &_name) const {
  auto res = find(_name);

  if (res == end()) {
    // If no match, try to use demerescau-levenshtein distance to make a helpful
    // suggestion
    auto it = begin();
    int min_dist(-1);
    for (; it != end(); ++it) {
      int dist = dl_string_dist(_name, it->name());
      if (min_dist < 0 || dist < min_dist) {
        min_dist = dist;

        res = it;
      }
    }

    throw std::runtime_error(
        "Invalid " + T::class_desc() + " \"" + _name + "\" specified.\n" +
        "                Did you mean \"" + res->name() + "\"?\n");
  }
  // std::cout << "Returning result of 'lookup' for " << _name << ". result is "
  // << res->name() << std::endl;
  return *res;  //->clone();
}

template <typename T>
void ParsingDictionary<T>::print_help(std::ostream &_stream,
                                      std::function<bool(T const &)> filter,
                                      int width, int separation) const {
  const_iterator it_begin(this->cbegin()), it_end(this->cend());
  std::string::size_type len(0);
  for (auto it = it_begin; it != it_end; ++it) {
    if (filter(*it)) len = max(len, it->name().size());
  }
  for (auto it = it_begin; it != it_end; ++it) {
    if (!filter(*it)) continue;
    _stream << std::string(5, ' ') << it->name()
            << std::string(len - it->name().size() + separation, ' ');
    std::string::size_type wcount(0);
    std::string::const_iterator str_end(it->description().cend());
    for (std::string::const_iterator str_it = it->description().cbegin();
         str_it != str_end; ++str_it) {
      if (wcount >= width && isspace(*str_it)) {
        _stream << std::endl << std::string(5 + len + separation, ' ');
        wcount = 0;
      } else {
        _stream << *str_it;
        wcount++;
      }
    }
    _stream << std::endl << std::endl;
  }
}

}  // namespace CASM
#endif
