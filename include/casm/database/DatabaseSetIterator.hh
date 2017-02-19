
namespace CASM {

  namespace DB {

    /// \brief DatabaseIterator for implementations using std::set<ValueType>
    ///
    /// DatabaseIterators must implement public methods:
    /// - Default constructor
    /// - std::unique_ptr<DatabaseIteratorBase<ValueType> > clone() const
    ///
    /// DatabaseIterators must implement private methods:
    /// - bool is_end() const
    /// - void increment()
    /// - reference dereference() const
    /// - DatabaseIteratorBase *_clone() const
    ///
    template<typename ValueType, typename DatabaseType>
    class DatabaseSetIterator : public DatabaseIteratorBase<ValueType> {

    public:

      DatabaseSetIterator() {}

      std::unique_ptr<DatabaseIteratorBase<ValueType> > clone() const {
        return std::unique_ptr<DatabaseIteratorBase<ValueType> >(this->_clone());
      }

    private:

      friend DatabaseType;

      typedef std::set<ValueType>::iterator base_iterator;

      DatabaseSetIterator(base_iterator _it, base_iterator _end) :
        m_it(_it),
        m_end(_end) {}

      base_iterator base() const {
        return m_it;
      }

      bool is_end() const override {
        return m_it == m_end;
      }

      void increment() override {
        ++m_it;
      }

      reference dereference() const override {
        return *m_it;
      }

      DatabaseSetIterator *_clone() const override {
        return new DatabaseSetIterator(*this);
      }

      base_iterator m_it;
      base_iterator m_end;
    };

  }
}

