#ifndef CASM_TestCompN
#define CASM_TestCompN

#include "casm/clex/ConfigIO.hh"

extern "C" {
  CASM::BaseDatumFormatter<CASM::Configuration> *make_TestCompN_formatter();
}

namespace CASM {

  namespace ConfigIO {

    /// \brief Calculate number of each species per unit cell
    ///
    /// \ingroup ConfigIO
    class TestCompN : public ConfigIO_impl::MolDependent {

    public:

      static const std::string Name;

      static const std::string Desc;


      TestCompN() :
        MolDependent(Name, Desc) {}


      // --- Required implementations -----------

      /// \brief Returns the parametric composition
      Eigen::VectorXd evaluate(const Configuration &config) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<TestCompN> clone() const {
        return std::unique_ptr<TestCompN>(this->_clone());
      }

    private:

      /// \brief Clone using copy constructor
      TestCompN *_clone() const override {
        return new TestCompN(*this);
      }

    };

  }

}

#endif
