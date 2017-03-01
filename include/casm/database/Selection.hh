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

      typedef std::map<std::string, bool>::const_iterator base_iterator;

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
        return *(m_list->db().find(m_it->first));
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
    class Selection : public Logging {

    public:

      typedef SelectionIterator<ObjType> iterator;
      typedef Index size_type;

      Selection(Database<ObjType> &_db) :
        Logging(_db.primclex()),
        m_db(&_db),
        m_primclex(&m_db->primclex()) {}

      Selection(Database<ObjType> &_db, fs::path selection_path = "MASTER");

      const PrimClex &primclex() const {
        return *m_primclex;
      }

      Database<ObjType> &db() const {
        return *m_db;
      }

      boost::iterator_range<iterator> all() const {
        return boost::make_iterator_range(
                 iterator(*this, m_data.begin(), false),
                 iterator(*this, m_data.end(), false));
      }

      boost::iterator_range<iterator> selected() const {
        return boost::make_iterator_range(
                 iterator(*this, m_data.begin(), true),
                 iterator(*this, m_data.end(), true));
      }


      std::map<std::string, bool> &data() {
        return m_data;
      }

      const std::map<std::string, bool> &data() const {
        return m_data;
      }

      Index size() const {
        return m_data.size();
      }

      Index selected_size() const;



      const std::vector<std::string> &col_headers() const {
        return m_col_headers;
      }

      const std::string &name() const {
        return m_name;
      }


      /// \brief True if obj is in Selection and is selected; false otherwise
      bool selected(const std::string &name_or_alias) const;

      /// \brief If obj is in Selection, set selected to specified value
      void set_selected(const std::string &name_or_alias, bool value) const;

      /// \brief Set selected objects to value of criteria
      void set(const DataFormatterDictionary<ObjType> &dict,
               const std::string &criteria);

      /// \brief Set selected objects to value, if criteria true
      void set(const DataFormatterDictionary<ObjType> &dict,
               const std::string &criteria,
               bool selected);

      void read(std::istream &_input);

      void print(const DataFormatterDictionary<ObjType> &_dict,
                 std::ostream &_out,
                 bool only_selected = false) const;

      const jsonParser &from_json(const jsonParser &_json);

      jsonParser &to_json(const DataFormatterDictionary<ObjType> &_dict,
                          jsonParser &_json,
                          bool only_selected = false) const;

      bool write(const DataFormatterDictionary<ObjType> &dict,
                 bool force,
                 const fs::path &out_path,
                 bool write_json,
                 bool only_selected) const;


    private:

      Database<ObjType> *m_db;
      const PrimClex *m_primclex;

      // first will only be 'name', no matter whether 'name' or 'alias' is
      // written in the selection file
      std::map<std::string, bool> m_data;

      std::vector<std::string> m_col_headers;
      std::string m_name;

    };

  }

}

#endif
