#ifndef CASM_SubClusterGenerator
#define CASM_SubClusterGenerator


namespace CASM {


  /// \brief Generates subclusters of a cluster with an iterator-like interface
  ///
  /// - Includes the null cluster and the original cluster
  ///
  template<typename ClusterType>
  class SubClusterGenerator :
    public boost::iterator_facade <
    SubClusterIterator,
    ClusterType,
    boost::forward_traversal_tag > {

  public:

    /// \brief Default construtor
    SubClusterGenerator() {}

    /// \brief Construt with the cluster to find subclusters of
    explicit SubClusterGenerator(const ClusterType &clust) :
      m_cluster(std::make_unique<ClusterType>(clust)),
      m_current(std::make_unique<ClusterType>(clust)),
      m_current_valid(false),
      m_site_counter(std::vector<int>(m_cluster->size(), 0),
                     std::vector<int>(m_cluster->size(), 1),
                     std::vector<int>(m_cluster->size(), 1)) {}

  private:

    void increment() {
      ++m_site_counter;
      m_current_valid = false;
    }

    bool equal(const SubClusterGenerator &other) const {
      if(valid() && other.valid()) {
        return std::equal(m_site_counter.value_begin(),
                          m_site_counter.value_end(),
                          other.m_site_counter.value_begin());
      }
      return false;
    }

    ClusterType &dereference() const {
      // lazy construction of 'm_current'
      if(!m_current_valid) {
        m_current->elements().clear();
        for(Index i = 0; i < m_site_counter.size(); ++i) {
          if(m_site_counter()[i]) {
            m_current->elements().push_back(m_cluster->element(i));
          }
        }
      }
      return *m_current;
    }

    bool valid() const {
      if(!m_cluster) {
        return false;
      }
      return m_site_counter.valid();
    }

  private:

    /// the cluster we're finding subclusters of
    notstd::cloneable_ptr<ClusterType> m_cluster;

    /// for lazy construction of m_current
    bool m_current_valid;

    /// The current subcluster
    notstd::cloneable_ptr<ClusterType> m_current;

    /// Indicates which sites to include (1) or not include (0) in the subcluster
    Counter<std::vector<int> > m_site_counter;
  };
}

#endif
