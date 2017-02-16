

namespace CASM {

  namespace DB {

    template<typename ObjType>
    class List {

    public:

      typedef CASM_TMP::traits<ObjType>::db_iterator iterator;
      typedef CASM_TMP::traits<ObjType>::db_const_iterator const_iterator;
      typedef CASM_TMP::traits<ObjType>::db_reverse_iterator reverse_iterator;
      typedef CASM_TMP::traits<ObjType>::db_const_reverse_iterator const_reverse_iterator;
      typedef CASM_TMP::traits<ObjType>::db_name_type name_type;


      List() :
        m_all(this),
        m_selected(this) {}

      List(Database<ObjType> db) :
        m_db(db),
        m_all(this),
        m_selected(this) {}

      List(Database<ObjType> db, fs::path path);


      Database<ObjType> &db() {
        return m_db;
      }
      const Database<ObjType> &db() const {
        return m_db;
      }


      std::map<name_type, bool> &data() {
        return m_data;
      }

      const std::map<name_type, bool> &data() const {
        return m_data;
      }


      const All &all() {
        return m_all;
      }
      const All &all() const {
        return m_all;
      }

      const Selected &selected() {
        return m_selected;
      }
      const Selected &selected() const {
        return m_selected;
      }

      bool is_selected(const name_type &name) const {
        return data().find(name_type)->second;
      }

      void set_selected(const name_type &name_type, bool selected) {
        m_data.find(name_type)->second = selected;
      }

    private:

      Database<ObjType> m_db;
      std::map<name_type, bool> m_data;
      All m_all;
      Selected m_selected;

    };


    /// \brief Actions on only selected objects in this list
    template<typename ObjType>
    class Selected {

    public:

      Selected(List *list) :
        m_list(list) {}


      iterator begin() {
        return iterator(m_list, m_list->begin());
      }

      const_iterator begin() const {
        return iterator(m_list, m_list->begin());
      }

      const_iterator cbegin() const {
        return iterator(m_list, m_list->begin());
      }


      iterator end() {
        return iterator(m_list, m_list->end());
      }

      iterator end() const {
        return iterator(m_list, m_list->end());
      }

      const_iterator cend() const {
        return iterator(m_list, m_list->end());
      }


      reverse_iterator rbegin() {
        return reverse_iterator(end());
      }

      const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
      }

      const_reverse_iterator crbegin() const {
        return const_reverse_iterator(end());
      }


      reverse_iterator rend() {
        return reverse_iterator(begin());
      }

      const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
      }

      const_reverse_iterator crend() const {
        return const_reverse_iterator(begin());
      }

      /// \brief Return number of selected configurations
      ///
      /// - Sums number of selected configurations each time called
      size_type size() const {
        size_type s {0};
        for(auto it = m_list->begin(); it != m_list->end(); ++it) {
          s += (size_type) it->second;
        }
        return s;
      }

      iterator find(const name_type &name_type) const {
        auto pos = m_list->find(name_type);
        if(pos->second) {
          return iterator(m_list, pos);
        }
        return end();
      }

    private:

      List *m_list;

    };

    /// \brief Actions on all objects in this list (selected or unselected)
    template<typename ObjType>
    class All {

    public:

      All(List *list) :
        m_list(list) {}


      iterator begin() {
        return iterator(m_list, m_list->begin());
      }

      const_iterator begin() const {
        return iterator(m_list, m_list->begin());
      }

      const_iterator cbegin() const {
        return iterator(m_list, m_list->begin());
      }


      iterator end() {
        return iterator(m_list, m_list->end());
      }

      iterator end() const {
        return iterator(m_list, m_list->end());
      }

      const_iterator cend() const {
        return iterator(m_list, m_list->end());
      }


      reverse_iterator rbegin() {
        return reverse_iterator(end());
      }

      const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
      }

      const_reverse_iterator crbegin() const {
        return const_reverse_iterator(end());
      }


      reverse_iterator rend() {
        return reverse_iterator(begin());
      }

      const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
      }

      const_reverse_iterator crend() const {
        return const_reverse_iterator(begin());
      }


      size_type size() const {
        return m_list->size();
      }

      iterator find(const name_type &name_type) const {
        auto pos = m_list->find(config);
        if(pos != m_list->end()) {
          return iterator(m_list, pos);
        }
        return end();
      }

    private:

      List *m_list;

    };

  }

}
