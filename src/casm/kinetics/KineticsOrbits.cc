
#define LOCAL_CLUSTER_ORBITS_INST(INSERTER) \
  \
  template INSERTER make_local_orbits<INSERTER>( \
    const Kinetics::DiffusionTransformation &diff_trans, \
    const std::vector<double> &cutoff_radius, \
    const std::vector<double> &max_length, \
    const std::vector<IntegralCluster> &custom_generators, \
    const std::function<bool (Site)> &site_filter, \
    double xtal_tol, \
    INSERTER result, \
    std::ostream &status, \
    const SymGroup &generating_group); \
  \
  template INSERTER make_local_orbits<INSERTER>( \
    const Kinetics::DiffusionTransformation &diff_trans, \
    const jsonParser &bspecs, \
    const std::function<bool (Site)> &site_filter, \
    double xtal_tol, \
    INSERTER result, \
    std::ostream &status, \
    const SymGroup &generating_group); \
  \
