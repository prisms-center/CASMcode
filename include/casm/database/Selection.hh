#ifndef CASM_Selection
#define CASM_Selection

#include <map>
#include "casm/external/boost.hh"
#include "casm/database/Database.hh"

namespace CASM {

  namespace DB {

    template<typename ObjType>
    class SelectionIterator :
      public boost::iterator_facade <
      SelectionIterator<ObjType>,
      ObjType,
      boost::forward_iterator_tag,
      const ObjType &,
      long > {

    public:

      /// Default constructor (equals end)
      DatabaseIterator() {}

    private:

      friend boost::iterator_core_access;

      friend Selection<ObjType>;

      typedef std::map<name_type, bool>::iterator base_iterator;

      /// Construct iterator
      DatabaseIterator(const Selection<ObjType> &_list, base_iterator _it, bool _selected_only) :
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
      reference dereference() const {
        return *(m_list->db()->find(m_it->first));
      }

      /// boost::iterator_facade implementation
      bool equal(const SelectionIterator &B) const {

        bool A_is_end = (m_db == nullptr) || (m_it == m_list->data().end());
        bool B_is_end = (B.m_db == nullptr) || (B.m_it == B.m_list->data().end());

        if(A_is_end != B_is_end) {
          return false;
        }

        if(A_is_end) {
          return true;
        }

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
    template<typename ObjType>
    class Selection {

    public:

      typedef SelectionIterator<ObjType> iterator;
      typedef Database<ObjType>::name_type name_type;
      typedef Index size_type;

      Selection(Database<ObjType> &_db) :
        m_db(&_db) {}

      Selection(Database<ObjType> &_db, fs::path selection_path = "MASTER");


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


      std::map<name_type, bool> &data() {
        return m_data;
      }

      const std::map<name_type, bool> &data() const {
        return m_data;
      }


      const std::vector<std::string> &col_headers() const {
        return m_col_headers;
      }

      const std::string &name() const {
        return m_name;
      }


      void read(std::istream &_input);

      void print(const DataFormatterDictionary<ObjType> &_dict,
                 std::ostream &_out,
                 bool only_selected = false) const;


      const jsonParser &from_json(const jsonParser &_json);

      jsonParser &to_json(const DataFormatterDictionary<ObjType> &_dict,
                          jsonParser &_json,
                          bool only_selected = false) const;

    private:

      Database<ObjType> *m_db;
      std::map<name_type, bool> m_data;
      std::vector<std::string> m_col_headers;
      std::string m_name;

    };

  }

}

#endif
