#ifndef CASM_Norm
#define CASM_Norm

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  /** \ingroup Reference
   *  @{
   */

  template<typename DataObject>
  class Norm {

  public:

    virtual ~Norm() {}

    /// \brief Default normalization is 1.0
    virtual double operator()(const DataObject &obj) const {
      return 1.0;
    }

    std::unique_ptr<Norm> clone() const {
      return std::unique_ptr<Norm>(this->_clone());
    }

  private:

    virtual Norm *_clone() const {
      return new Norm(*this);
    }
  };

  class NormPerUnitCell : public Norm<Configuration> {

  public:

    /// \brief Return configuration supercell size
    double operator()(const Configuration &config) const override {
      return config.get_supercell().volume();
    }

    /// \brief Clone
    std::unique_ptr<NormPerUnitCell> clone() const {
      return notstd::make_unique<NormPerUnitCell>(*this);
    }

  private:

    NormPerUnitCell *_clone() const override {
      return new NormPerUnitCell(*this);
    }
  };

  class NormPerSpecies : public Norm<Configuration> {

  public:

    /// \brief Return number of non-Va species in configuration per unitcell
    double operator()(const Configuration &config) const override {
      return n_species(config);
    }

    /// \brief Clone
    std::unique_ptr<NormPerSpecies> clone() const {
      return notstd::make_unique<NormPerSpecies>(*this);
    }

  private:

    NormPerSpecies *_clone() const override {
      return new NormPerSpecies(*this);
    }
  };

  /** @} */
}

#endif

