#ifndef CASM_Selection
#define CASM_Selection

#include <map>
#include "casm/external/boost.hh"
#include "casm/database/Database.hh"

namespace CASM {

  namespace DB {

    template<typename ObjType>
    class Selection;

    template<typename ObjType>
    class SelectionIterator :
      public boost::iterator_facade <
      SelectionIterator<ObjType>,
      ObjType,
      std::forward_iterator_tag,
      const ObjType &,
      long > {

    public:

      /// Default constructor (equals end)
      SelectionIterator() {}

    private:

      friend boost::iterator_core_access;

      friend Selection<ObjType>;

      typedef std::map<std::string, bool>::iterator base_iterator;

      /// Construct iterator
      SelectionIterator(const Selection<ObjType> &_list, base_iterator _it, bool _selected_only) :
        m_list(&_list),
        m_it(_it),
        m_selected_only(_selected_only) {}

      /// boost::iterator_facade implementation
      void increment() {
        ++m_it;
        while(m_selected_only && m_it != m_list->data().end() && m_it->second == false) {
          ++m_it;
        }
      }

      /// boost::iterator_facade implementation
      const ObjType &dereference() const {
        return *(m_list->db()->find(m_it->first));
      }

      /// boost::iterator_facade implementation
      bool equal(const SelectionIterator &B) const {
        return m_it == B.m_it;
      }

      const Selection<ObjType> *m_list;
      base_iterator m_it;
      bool m_selected_only;
    };


    /// Provides a means of selecting and iterating over a subset of objects in
    /// a database
    ///
    /// - Stores a std::map of ObjType name_or_alias -> bool (true if selected)
    /// - Provides iterators over either all objects in list, or only selected
    ///   objects
    /// - All iterators provide const references to ObjType
    ///
    /// - ObjType must have:
    ///   - std::string name()
    ///   - std::string alias()
    ///
    template<typename ObjType>
    class Selection {

    public:

      typedef SelectionIterator<ObjType> iterator;
      typedef Index size_type;

      Selection(const Database<ObjType> &_db) :
        m_db(&_db) {}

      Selection(const Database<ObjType> &_db, fs::path selection_path = "MASTER");


      Database<ObjType> &db() {
        return *m_db;
      }
      const Database<ObjType> &db() const {
        return *m_db;
      }

      boost::iterator_range<iterator> all() const {
        return boost::make_iterator_range(
                 iterator(*this, m_data.begin(), false),
                 iterator(*this, m_data.end(), false));
      }

      boost::iterator_range<iterator> selected() const {
        return boost::make_iterator_range(
                 _iterator(*this, m_data.begin(), true),
                 _iterator(*this, m_data.end(), true));
      }


      std::map<std::string, bool> &data() {
        return m_data;
      }

      const std::map<std::string, bool> &data() const {
        return m_data;
      }


      const std::vector<std::string> &col_headers() const {
        return m_col_headers;
      }

      const std::string &name() const {
        return m_name;
      }


      /// \brief True if obj is in Selection and is selected; false otherwise
      bool selected(const ObjType &obj) const;

      /// \brief Ensure obj is in Selection and set selected to specified value
      void set_selected(const ObjType &obj, bool selected) const;

      void read(std::istream &_input);

      void print(const DataFormatterDictionary<ObjType> &_dict,
                 std::ostream &_out,
                 bool only_selected = false) const;


      const jsonParser &from_json(const jsonParser &_json);

      jsonParser &to_json(const DataFormatterDictionary<ObjType> &_dict,
                          jsonParser &_json,
                          bool only_selected = false) const;

    private:

      const Database<ObjType> *m_db;

      // first may be 'name' or 'alias'
      std::map<std::string, bool> m_data;
      std::vector<std::string> m_col_headers;
      std::string m_name;

    };

  }

}

#endif
