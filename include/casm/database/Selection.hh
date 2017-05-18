#ifndef CASM_Selection
#define CASM_Selection

#include <map>
#include "casm/external/boost.hh"
#include "casm/database/Database.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {
  namespace DB {

    template<typename ObjType>
    class Selection;

    template<typename ObjType, typename BaseIterator>
    class SelectionIterator :
      public boost::iterator_facade <
      SelectionIterator<ObjType, BaseIterator>,
      ObjType,
      std::bidirectional_iterator_tag,
      const ObjType &,
      long > {

    public:

      typedef typename CASM_TMP::ConstSwitch <
      std::is_same<BaseIterator, std::map<std::string, bool>::const_iterator>::value,
          bool > bool_type;

      /// Default constructor (equals end)
      SelectionIterator() {}

      /// \brief Name of object the iterator points at
      std::string name() const {
        return m_it->first;
      }

      /// \brief Reference to value 'is_selected'
      bool_type &is_selected() {
        return m_it->second;
      }

      /// \brief Reference to value 'is_selected'
      bool is_selected() const {
        return m_it->second;
      }

    private:

      friend boost::iterator_core_access;

      friend Selection<ObjType>;

      /// Construct iterator
      SelectionIterator(const Selection<ObjType> &_list, BaseIterator _it, bool _selected_only) :
        m_list(&_list),
        m_it(_it),
        m_selected_only(_selected_only) {
        if(m_selected_only && m_it != m_list->data().end() && m_it->second == false) {
          increment();
        }
      }

      /// boost::iterator_facade implementation
      void increment() {
        ++m_it;
        while(m_selected_only && m_it != m_list->data().end() && m_it->second == false) {
          ++m_it;
        }
      }

      /// boost::iterator_facade implementation
      void decrement() {
        --m_it;
        while(m_selected_only && m_it != m_list->data().begin() && m_it->second == false) {
          --m_it;
        }
      }

      /// boost::iterator_facade implementation
      const ObjType &dereference() const;

      /// boost::iterator_facade implementation
      bool equal(const SelectionIterator &B) const {
        return m_it == B.m_it;
      }

      /*
      long distance_to(const SelectionIterator &B) const {
        Do not define this or boost::iterator_range<T>::size will compile but
        cause runtime errors. Use boost::distance instead.
      }
      */

      const Selection<ObjType> *m_list;
      BaseIterator m_it;
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

      struct Compare {
        bool operator()(std::string A, std::string B) const {
          return traits<ObjType>::name_compare(A, B);
        }
      };

      typedef std::map<std::string, bool, Compare> map_type;
      typedef typename map_type::iterator base_iterator;
      typedef typename map_type::const_iterator base_const_iterator;
      typedef SelectionIterator<ObjType, base_iterator> iterator;
      typedef SelectionIterator<ObjType, base_const_iterator> const_iterator;
      typedef Index size_type;

      /// \brief Default construct into invalid state
      Selection() {};

      /// \brief Use default ObjType database
      Selection(const PrimClex &_primclex, fs::path selection_path = "MASTER");

      /// \brief Use specified ObjType database
      Selection(Database<ObjType> &_db, fs::path selection_path = "MASTER");

      const PrimClex &primclex() const {
        return *m_primclex;
      }

      Database<ObjType> &db() const {
        return *m_db;
      }

      boost::iterator_range<iterator> all() {
        return boost::make_iterator_range(
                 iterator(*this, m_data.begin(), false),
                 iterator(*this, m_data.end(), false));
      }

      boost::iterator_range<const_iterator> all() const {
        return boost::make_iterator_range(
                 const_iterator(*this, m_data.begin(), false),
                 const_iterator(*this, m_data.end(), false));
      }

      boost::iterator_range<iterator> selected() {
        return boost::make_iterator_range(
                 iterator(*this, m_data.begin(), true),
                 iterator(*this, m_data.end(), true));
      }

      boost::iterator_range<const_iterator> selected() const {
        return boost::make_iterator_range(
                 const_iterator(*this, m_data.begin(), true),
                 const_iterator(*this, m_data.end(), true));
      }


      map_type &data() {
        return m_data;
      }

      const map_type &data() const {
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
      bool is_selected(const std::string &name_or_alias) const;

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
      map_type m_data;

      std::vector<std::string> m_col_headers;
      std::string m_name;

    };

  }
}

#endif
