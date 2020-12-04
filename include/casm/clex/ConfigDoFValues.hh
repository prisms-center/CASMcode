#ifndef CASM_ConfigDoFValues
#define CASM_ConfigDoFValues

#include "casm/basis_set/DoFSet.hh"
#include "casm/crystallography/AnisoValTraits.hh"
namespace CASM {

  class ConfigDoFValues {
  public:
    ConfigDoFValues() : m_n_sublat(0), m_n_vol(0)
    {}
    ConfigDoFValues(DoF::BasicTraits const &_traits, Index _n_sublat, Index _n_vol) :
      m_type(_traits.name()),
      m_n_sublat(_n_sublat),
      m_n_vol(_n_vol) {
    }

    std::string const &type_name() const {
      return m_type;
    }

    Index n_vol() const {
      return m_n_vol;
    }

    Index n_sublat() const {
      return m_n_sublat;
    }

    void resize_vol(Index _n_vol) {
      m_n_vol = _n_vol;
      _resize();
    }

  protected:
    virtual void _resize() = 0;

  private:

    DoFKey m_type;
    Index m_n_sublat;
    Index m_n_vol;
  };

  class LocalDiscreteConfigDoFValues : public ConfigDoFValues {
  public:

    typedef Eigen::VectorXi ValueType;
    typedef Eigen::VectorXi &Reference;
    typedef const Eigen::VectorXi &ConstReference;

    typedef typename ValueType::Scalar SiteValueType;
    typedef int &SiteReference;
    typedef const int &ConstSiteReference;

    typedef ValueType SublatValueType;
    typedef typename ValueType::SegmentReturnType SublatReference;
    typedef typename ValueType::ConstSegmentReturnType ConstSublatReference;

    LocalDiscreteConfigDoFValues() {}

    LocalDiscreteConfigDoFValues(DoF::BasicTraits const &_traits,
                                 Index _n_sublat,
                                 Index _n_vol,
                                 Eigen::Ref< const ValueType > const &_vals,
                                 std::vector<SymGroupRepID> const &_symrep_IDs) :
      ConfigDoFValues(_traits, _n_sublat, _n_vol),
      m_vals(_vals),
      m_symrep_IDs(_symrep_IDs) {

    }

    /// Access occupation values (values are indices into Site::occupant_dof())
    Reference values() {
      return m_vals;
    }

    /// Const access occupation values (values are indices into Site::occupant_dof())
    ConstReference values() const {
      return m_vals;
    }

    /// Access vector block of values for all sites on one sublattice
    SublatReference sublat(Index b) {
      return m_vals.segment(b * n_vol(), n_vol());
    }

    /// Const access vector block of values for all sites on one sublattice
    ConstSublatReference sublat(Index b) const {
      return m_vals.segment(b * n_vol(), n_vol());
    }

    /// Provides the symmetry representations for transforming `values` (i.e. due to molecule
    /// orientation, not for permuting sites)
    std::vector<SymGroupRepID> const &symrep_IDs() const {
      return m_symrep_IDs;
    }

  protected:
    void _resize() override {
      m_vals.resize(n_vol()*n_sublat());
    }

  private:
    ValueType m_vals;
    std::vector<SymGroupRepID> m_symrep_IDs;
  };


  class LocalContinuousConfigDoFValues : public ConfigDoFValues {
  public:

    typedef Eigen::MatrixXd ValueType;
    typedef Eigen::MatrixXd &Reference;
    typedef const Eigen::MatrixXd &ConstReference;

    typedef Eigen::VectorXd SiteValueType;
    typedef typename ValueType::ColXpr SiteReference;
    typedef const typename ValueType::ConstColXpr ConstSiteReference;

    typedef ValueType SublatValueType;
    typedef typename Eigen::Block<ValueType> SublatReference;
    typedef const typename Eigen::Block<const ValueType> ConstSublatReference;

    LocalContinuousConfigDoFValues() {}

    LocalContinuousConfigDoFValues(DoF::BasicTraits const &_traits,
                                   Index _n_sublat,
                                   Index _n_vol,
                                   Eigen::Ref< const ValueType > const &_vals,
                                   std::vector<DoFSetInfo> const &_info) :
      ConfigDoFValues(_traits, _n_sublat, _n_vol),
      m_vals(_vals),
      m_info(_info) {

    }

    // /// DoF vector representation size
    // Index dim() const {
    //   return m_vals.rows(); // this is the standard basis dimension, not the prim basis dimension
    // }

    /// Access site DoF values (prim DoF basis, matrix representing all sites)
    ///
    /// Notes:
    /// - Matrix of size rows=standard basis dimension (this->info()[b].basis().rows()), cols=n_vol()*n_sublat(),
    /// - Each column represents a site DoF value in the prim DoF basis
    /// - The prim DoF basis can be accessed by `this->info()[b].basis()`, where `b` is the sublattice index,
    ///   `b = column_index / this->n_vol()`
    /// - If the prim DoF basis dimension (this->info()[b].basis().cols()) is less than the standard
    ///   DoF basis dimension, the matrix includes blocks of zeros for the corresponding
    ///   sublattice.
    Reference values() {
      return m_vals;
    }

    /// Const access DoF values (prim DoF basis, matrix representing all sites)
    ///
    /// Notes:
    /// - Matrix of size rows=dim(), cols=n_vol()*n_sublat(),
    /// - Each column represents a site DoF value in the prim DoF basis
    /// - The prim DoF basis can be accessed by `this->info()[b]`, where `b` is the sublattice index,
    ///   `b = column_index / this->n_vol()`
    ConstReference values() const {
      return m_vals;
    }

    /// Set local DoF values from standard DoF values
    ///
    /// Notes:
    /// - Standard DoF values are those expressed according to the standard DoF basis, i.e.
    ///   coordinates whose values corrspond to AnisoValTraits::standard_var_names().
    void from_standard_values(Eigen::Ref<const Eigen::MatrixXd> const &_standard_values);

    /// Get local DoF values as standard DoF values
    ///
    /// Notes:
    /// - Standard DoF values are those expressed according to the standard DoF basis, i.e.
    ///   coordinates whose values corrspond to AnisoValTraits::standard_var_names().
    Eigen::MatrixXd standard_values() const {
      Index rows = m_info[0].basis().rows();
      Eigen::MatrixXd result(rows, m_vals.cols());
      for(Index b = 0; b < n_sublat(); ++b) {
        result.block(0, b * n_vol(), rows, n_vol()) = info()[b].basis() * sublat(b).topRows(info()[b].dim());
      }
      return result;
    }

    /// Access site DoF value (prim DoF basis, vector associated with a single site)
    ///
    /// Note:
    /// - If the prim DoF basis dimension < standard DoF basis dimension, this includes a tail of
    ///   zeros (for rows >= this->info()[b].dim()) that should not be modified.
    SiteReference site_value(Index l) {
      return m_vals.col(l);
    }

    /// Const access site DoF value (prim DoF basis, vector associated with a single site)
    /// Note:
    /// - If the prim DoF basis dimension < standard DoF basis dimension, this includes a tail of
    ///   zeros (for rows >= this->info()[b].dim()).
    ConstSiteReference site_value(Index l) const {
      return m_vals.col(l);
    }

    /// Access matrix block of values for all sites on one sublattice
    SublatReference sublat(Index b) {
      return m_vals.block(0, b * n_vol(), m_vals.rows(), n_vol());
    }

    /// Const access matrix block of values for all sites on one sublattice
    ConstSublatReference sublat(Index b) const {
      return m_vals.block(0, b * n_vol(), m_vals.rows(), n_vol());
    }

    /// DoFSetInfo provides the basis and symmetry representations for `values`
    std::vector<DoFSetInfo> const &info() const {
      return m_info;
    }


  protected:
    void _resize() override {
      m_vals.resize(m_vals.rows(), n_vol()*n_sublat());
    }

  private:
    ValueType m_vals;
    std::vector<DoFSetInfo> m_info;
  };


  // TODO: It might be confusing that ValueType, Reference, and ConstReference are in terms of
  // Eigen::VectorXd, but dim(), from_standard_values(), and standard_values() are written using
  // Eigen::MatrixXd. Since Eigen::VectorXd is a 1 column Eigen::MatrixXd it could be all written
  // in terms of Eigen::VectorXd or Eigen::MatrixXd. Not sure if it is more useful to express
  // in terms of Eigen::MatrixXd and match the LocalContinuousConfigDoFValues interface or use
  // Eigen::VectorXd consistently to express that the value is 1 dimensional.

  class GlobalContinuousConfigDoFValues : public ConfigDoFValues {
  public:

    typedef Eigen::VectorXd ValueType;
    typedef Eigen::VectorXd &Reference;
    typedef const Eigen::VectorXd &ConstReference;

    typedef typename ValueType::Scalar SiteValueType;
    typedef int &SiteReference;
    typedef const int &ConstSiteReference;

    GlobalContinuousConfigDoFValues() :
      m_info(SymGroupRepID(), Eigen::MatrixXd::Zero(0, 0)) {}

    GlobalContinuousConfigDoFValues(DoF::BasicTraits const &_traits,
                                    Index _n_sublat,
                                    Index _n_vol,
                                    Eigen::Ref< const ValueType > const &_vals,
                                    DoFSetInfo const &_info) :
      ConfigDoFValues(_traits, _n_sublat, _n_vol),
      m_vals(_vals),
      m_info(_info) {

    }

    /// Global DoF vector representation dimension
    Index dim() const {
      return m_vals.rows();
    }

    /// Access values (coordinates defined by xtal::SiteDoFSet::basis() / DoFSetInfo::basis() )
    Reference values() {
      return m_vals;
    }

    /// Const access values (coordinates defined by xtal::SiteDoFSet::basis() / DoFSetInfo::basis() )
    ConstReference values() const {
      return m_vals;
    }

    /// Set global DoF values from standard DoF values
    ///
    /// Notes:
    /// - Standard DoF values are those expressed according to the standard DoF basis, i.e.
    ///   coordinates whose values corrspond to AnisoValTraits::standard_var_names().
    void from_standard_values(Eigen::Ref<const Eigen::MatrixXd> const &_standard_values);

    /// Get global DoF values as standard DoF values
    ///
    /// Notes:
    /// - Standard DoF values are those expressed according to the standard DoF basis, i.e.
    ///   coordinates whose values corrspond to AnisoValTraits::standard_var_names().
    Eigen::MatrixXd standard_values() const {
      return m_info.basis() * m_vals;
    }

    /// DoFSetInfo provides the basis and symmetry representations for `values`
    DoFSetInfo const &info() const {
      return m_info;
    }

  protected:
    void _resize() override { }

  private:
    ValueType m_vals;
    DoFSetInfo m_info;
  };

}

#endif
