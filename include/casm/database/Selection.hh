#ifndef CASM_Selection
#define CASM_Selection

#include <map>
#include <boost/filesystem.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>
#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"
#include "casm/casm_io/enum/stream_io.hh"
#include "casm/casm_io/enum/json_io.hh"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {
  namespace DB {
    enum class SELECTION_TYPE {
      MASTER, ALL, NONE, EMPTY, CALCULATED
    };

    ENUM_IO_DECL(CASM::DB::SELECTION_TYPE)
  }
  template<> inline std::string singleline_enum_help<DB::SELECTION_TYPE>() {
    return standard_singleline_enum_help<DB::SELECTION_TYPE>(to_string(DB::SELECTION_TYPE::MASTER), "filename");
  }
  ENUM_TRAITS(DB::SELECTION_TYPE)
}

namespace CASM {

  class PrimClex;
  class jsonParser;

  namespace DB {
    template<typename ValueType>
    class Database;

    template<typename ValueType> class Database;

    template<typename ObjType> class Selection;

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
      std::string name() const;

      /// \brief Reference to value 'is_selected'
      bool_type &is_selected();

      /// \brief Reference to value 'is_selected'
      bool is_selected() const;

    private:

      friend boost::iterator_core_access;

      friend Selection<ObjType>;

      /// Construct iterator
      SelectionIterator(const Selection<ObjType> &_list, BaseIterator _it, bool _selected_only);

      /// boost::iterator_facade implementation
      void increment();

      /// boost::iterator_facade implementation
      void decrement();

      /// boost::iterator_facade implementation
      const ObjType &dereference() const;

      /// boost::iterator_facade implementation
      bool equal(const SelectionIterator &B) const;

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
    class Selection {

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
      Selection();

      /// \brief Use default ObjType database
      Selection(const PrimClex &_primclex, fs::path selection_path = "MASTER");

      /// \brief Use specified ObjType database
      Selection(Database<ObjType> &_db, fs::path selection_path = "MASTER");

      const PrimClex &primclex() const;

      Database<ObjType> &db() const;

      boost::iterator_range<iterator> all();

      boost::iterator_range<const_iterator> all() const;

      boost::iterator_range<iterator> selected();

      boost::iterator_range<const_iterator> selected() const;


      map_type &data();

      const map_type &data() const;

      Index size() const;

      Index selected_size() const;

      const std::vector<std::string> &col_headers() const;

      const std::string &name() const;


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
