#include "casm/kinetics/DoFTransformation.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  namespace Kinetics {

    // DoFTransformation

    DoFTransformation::DoFTransformation(const PrimType &prim) :
      m_prim(&prim) {};

    /// \brief Return a reference to the primitive Structure lattice
    const Lattice &DoFTransformation::lattice() const {
      return prim().lattice();
    }

    Configuration &DoFTransformation::apply_to(Configuration &config) const {
      return this->apply_to_impl(config);
    }

    Configuration &DoFTransformation::apply_reverse_to(Configuration &config) const {
      return this->apply_reverse_to_impl(config);
    }

    DoFTransformation &DoFTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(op);
      return *this;
    }

    void DoFTransformation::reverse() {
      this->reverse_impl();
    }

    std::unique_ptr<DoFTransformation> DoFTransformation::clone() const {
      return std::unique_ptr<DoFTransformation>(this->_clone());
    }


    // SiteDoFTransformation

    SiteDoFTransformation::SiteDoFTransformation(const UnitCellCoord &_uccoord) :
      DoFTransformation(_uccoord.unit()),
      uccoord(_uccoord) {}

    SiteDoFTransformation::~SiteDoFTransformation() {}

    SiteDoFTransformation &SiteDoFTransformation::operator+=(UnitCell frac) {
      uccoord += frac;
      return *this;
    }

    SiteDoFTransformation &SiteDoFTransformation::operator-=(UnitCell frac) {
      uccoord -= frac;
      return *this;
    }

    SiteDoFTransformation &SiteDoFTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(op);
      return *this;
    }

    std::unique_ptr<SiteDoFTransformation> SiteDoFTransformation::clone() const {
      return std::unique_ptr<SiteDoFTransformation>(this->_clone());
    }

    void SiteDoFTransformation::apply_sym_impl(const SymOp &op) {
      uccoord.apply_sym(op);
    }
  }
}

