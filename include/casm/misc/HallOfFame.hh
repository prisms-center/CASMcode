#ifndef CASM_HallOfFame
#define CASM_HallOfFame

#include <set>

#include "casm/misc/CASM_math.hh"

namespace CASM {

  /// \brief A container for storing best scoring objects
  ///
  /// Essentially a std::map with insert modified to calculate the 'score'
  /// of an object and insert while keeping only a particular number of best
  /// scoring objects.
  ///
  /// ObjectCompare implements comparison of ObjectType objects
  ///
  template<typename ObjectType, typename Metric, typename ObjectCompare = std::less<ObjectType> >
  class HallOfFame {

  public:

    typedef std::pair<double, ObjectType> PairType;

    /// Compare PairType lexicographically using FloatCompare for score,
    /// and ObjectCompare for objects
    class Compare {

    public:

      Compare(ObjectCompare _obj_compare, double _score_tol) :
        m_score_compare(_score_tol), m_obj_compare(_obj_compare) {}

      bool operator()(const PairType &A, const PairType &B) const {

        if(m_score_compare(A.first, B.first)) {
          return true;
        }
        if(m_score_compare(B.first, A.first)) {
          return false;
        }
        return m_obj_compare(A.second, B.second);
      }

    private:

      FloatCompare m_score_compare;
      ObjectCompare m_obj_compare;

    };

    typedef std::set<PairType, Compare> ContainerType;
    typedef typename ContainerType::const_iterator const_iterator;
    typedef typename ContainerType::const_reverse_iterator const_reverse_iterator;
    typedef Index size_type;

    /// Results data structure for HallOfFame::insert
    struct InsertResult {

      typedef HallOfFame::const_iterator const_iterator;

      InsertResult(const_iterator _pos, bool _success, double _score, bool _excluded, const_iterator _excluded_pos) :
        pos(_pos), success(_success), score(_score), excluded(_excluded), excluded_pos(_excluded_pos) {}

      InsertResult(std::pair<const_iterator, bool> _res, double _score, bool _excluded, const_iterator _excluded_pos) :
        pos(_res.first), success(_res.second), score(_score), excluded(_excluded), excluded_pos(_excluded_pos)  {}

      const_iterator pos;
      bool success;
      double score;
      bool excluded;
      const_iterator excluded_pos;
    };

    /// \brief Constructor
    ///
    /// \param _metric A function scoring 'ObjectType' objects
    /// \param _obj_compare Functor implementing comparison of 'ObjectType' objects
    /// \param _max_size The maximum number of 'ObjectType' objects to store in the HallOfFame
    /// \param _metric_tol The tolerance used for comparing scores, via FloatCompare
    ///
    HallOfFame(Metric _metric, ObjectCompare _obj_compare = std::less<ObjectType>(), Index _max_size = 100, double _score_tol = TOL):
      m_metric(_metric),
      m_obj_compare(_obj_compare),
      m_score_compare(_score_tol),
      m_halloffame(Compare(m_obj_compare, _score_tol)),
      m_check_exclude(false),
      m_exclude(Compare(m_obj_compare, _score_tol)),
      m_max_size(_max_size) {}

    /// \brief Add an object that should not be included in the hall of fame
    void exclude(const ObjectType &obj) {
      m_exclude.insert(PairType(m_metric(obj), obj));
      m_check_exclude = true;
    }

    /// \brief Add objects that should not be included in the hall of fame
    template<typename ObjectIterator>
    void exclude(ObjectIterator begin, ObjectIterator end) {
      for(auto it = begin; it != end; ++it) {
        exclude(*it);
      }
    }

    void clear_excluded() {
      m_exclude.clear();
    }

    Metric metric() const {
      return m_metric;
    }

    bool empty() const {
      return m_halloffame.empty();
    }

    size_type size() const {
      return m_halloffame.size();
    }

    size_type max_size() const {
      return m_max_size;
    }

    /// \brief Insert object in HallOfFame
    ///
    /// - Will score object, and insert into HallOfFame, erasing the worst scoring
    /// member if necessary to maintain the max_size specified at construction
    InsertResult insert(const ObjectType &obj) {

      auto excluded_pos = m_exclude.end();
      if(m_max_size <= 0) {
        return InsertResult(m_halloffame.end(), false, std::numeric_limits<double>::quiet_NaN(), false, excluded_pos);
      }

      double score = m_metric(obj);

      // if score is not good enough for hall of fame, do not insert
      if(m_halloffame.size() == m_max_size && m_score_compare(m_halloffame.rbegin()->first, score)) {
        return InsertResult(m_halloffame.end(), false, score, false, excluded_pos);
      }

      PairType test(score, obj);

      // insert
      auto res = m_halloffame.insert(test);

      // if not successful because it already exists in the hall of fame
      if(!res.second) {
        return InsertResult(res, score, false, excluded_pos);
      }

      // if successful, check if in exclude list
      if(m_check_exclude) {
        excluded_pos = m_exclude.find(test);
        if(excluded_pos != m_exclude.end()) {
          m_halloffame.erase(res.first);
          return InsertResult(m_halloffame.end(), false, score, true, excluded_pos);
        }
      }

      // remove any extras
      if(m_halloffame.size() > m_max_size) {
        m_halloffame.erase(std::prev(m_halloffame.end()));
      }

      return InsertResult(res, score, false, excluded_pos);
    }

    void clear() {
      m_halloffame.clear();
    }

    const_iterator begin() const {
      return m_halloffame.begin();
    }

    const_iterator cbegin() const {
      return m_halloffame.cbegin();
    }

    const_iterator end() const {
      return m_halloffame.end();
    }

    const_iterator cend() const {
      return m_halloffame.cend();
    }

    const_reverse_iterator rbegin() const {
      return m_halloffame.rbegin();
    }

    const_reverse_iterator crbegin() const {
      return m_halloffame.crbegin();
    }

    const_reverse_iterator rend() const {
      return m_halloffame.rend();
    }

    const_reverse_iterator crend() const {
      return m_halloffame.crend();
    }

    std::pair<const_iterator, const_iterator> equal_range(double key) const {
      return m_halloffame.equal_range(key);
    }

    size_type count(double key) const {
      return m_halloffame.count(key);
    }

    const_iterator find(double key) const {
      return m_halloffame.find(key);
    }

    const_iterator lower_bound(double key) const {
      return m_halloffame.lower_bound(key);
    }

    const_iterator upper_bound(double key) const {
      return m_halloffame.upper_bound(key);
    }


  private:

    Metric m_metric;
    ObjectCompare m_obj_compare;
    FloatCompare m_score_compare;
    ContainerType m_halloffame;
    bool m_check_exclude;
    ContainerType m_exclude;
    Index m_max_size;

  };

}

#endif
