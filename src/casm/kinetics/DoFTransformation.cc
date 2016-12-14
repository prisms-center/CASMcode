#include "casm/kinetics/DoFTransformation.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  namespace Kinetics {

    // DoFTransformation

    DoFTransformation::DoFTransformation(const PrimType &prim) :
      m_prim(&prim) {};

    virtual DoFTransformation::~DoFTransformation() {};

    Configuration &DoFTransformation::apply_to(Configuration &config) const {
      return this->apply_to_impl(config);
    }

    Configuration &DoFTransformation::apply_reverse_to(Configuration &config) const {
      return this->apply_reverse_to_impl(config);
    }

    DoFTransformation &DoFTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(config);
      return *this;
    }

    void DoFTransformation::reverse() {
      this->reverse_impl();
    }

    std::unique_ptr<DoFTransformation> DoFTransformation::clone() const {
      return std::unique_ptr<DoFTransformation>(this->_clone());
    }


    // SiteDoFTransformation

    SiteDoFTransformation::SiteDoFTransformation(const UnitCellCoord &uccoord) :
      DoFTransformation(uccoord.unit()),
      m_uccoord(uccoord) {}

    SiteDoFTransformation::~SiteDoFTransformation() {}

    UnitCellCoord &SiteDoFTransformation::uccoord() {
      return m_uccoord;
    }

    const UnitCellCoord &SiteDoFTransformation::uccoord() const {
      return m_uccoord;
    }

    SiteDoFTransformation &SiteDoFTransformation::operator+=(UnitCell frac) {
      m_uccoord += frac;
      return *this;
    }

    SiteDoFTransformation &SiteDoFTransformation::operator-=(UnitCell frac) {
      m_uccoord -= frac;
      return *this;
    }

    SiteDoFTransformation &SiteDoFTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(config);
      return *this;
    }

    std::unique_ptr<SiteDoFTransformation> SiteDoFTransformation::clone() const {
      return std::unique_ptr<SiteDoFTransformation>(this->_clone());
    }

    void SiteDoFTransformation::apply_sym_impl(const SymOp &op) {
      m_uccoord.apply_sym(op);
    }


    // OccupationTransformation

    OccupationTransformation::OccupationTransformation(
      const UnitCellCoord &uccoord,
      Index from_value,
      Index to_value) :
      SiteDoFTransformation(uccoord),
      m_from_value(from_value),
      m_to_value(to_value) {}

    Index &OccupationTransformation::from_value() {
      return m_from_value;
    }

    const Index &OccupationTransformation::from_value() const {
      return m_from_value;
    }

    Index &OccupationTransformation::to_value() {
      return m_to_value;
    }

    const Index &OccupationTransformation::to_value() const {
      return m_to_value;
    }

    OccupationTransformation &OccupationTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(config);
      return *this;
    }

    std::unique_ptr<OccupationTransformation> OccupationTransformation::clone() const {
      return std::unique_ptr<OccupationTransformation>(this->_clone());
    }

    Configuration &OccupationTransformation::apply_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord()), to_value());
      return config;
    }

    Configuration &OccupationTransformation::apply_reverse_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord()), from_value());
      return config;
    }

    void OccupationTransformation::apply_sym_impl(const SymOp &op) {
      static_cast<const SiteDoFTransformation &>(*this).apply_sym(op);

      //MOLECULE_SUPPORT: apply permutation to from_value & to_value
    }

    void OccupationTransformation::reverse_impl() {
      using std::swap;
      swap(m_from_value, m_to_value);
    }

    OccupationTransformation *OccupationTransformation::_clone() const {
      return new OccupationTransformation(*this);
    }


    // AtomTrajectory

    AtomTrajectory::AtomTrajectory(const UnitCellCoord &from_uccoord,
                                   Index from_value,
                                   Index from_atom_index,
                                   const UnitCellCoord &to_uccoord,
                                   Index to_value,
                                   Index to_atom_index) :
      m_from_uccoord(from_uccoord),
      m_from_value(from_value),
      m_from_atom_index(from_atom_index),
      m_to_uccoord(to_uccoord),
      m_to_value(to_value),
      m_to_atom_index(to_atom_index) {}


    UnitCellCoord &AtomTrajectory::from_uccoord() {
      return m_from_uccoord;
    }

    const UnitCellCoord &AtomTrajectory::from_uccoord() const {
      return m_from_uccoord;
    }

    UnitCellCoord &AtomTrajectory::to_uccoord() {
      return m_to_uccoord;
    }

    const UnitCellCoord &AtomTrajectory::to_uccoord() const {
      return m_to_uccoord;
    }


    Index &AtomTrajectory::from_value() {
      return m_from_value;
    }

    const Index &AtomTrajectory::from_value() const {
      return m_from_value;
    }

    Index &AtomTrajectory::to_value() {
      return m_to_value;
    }

    const Index &AtomTrajectory::to_value() const {
      return m_to_value;
    }


    Index &AtomTrajectory::from_atom_index() {
      return m_from_atom_index;
    }

    const Index &AtomTrajectory::from_atom_index() const {
      return m_from_atom_index;
    }

    Index &AtomTrajectory::to_atom_index() {
      return m_to_atom_index;
    }

    const Index &AtomTrajectory::to_atom_index() const {
      return m_to_atom_index;
    }

    AtomTrajectory &AtomTrajectory::operator+=(UnitCell frac) {
      m_from_uccoord += frac;
      m_to_uccoord += frac;
      return *this;
    }

    AtomTrajectory &AtomTrajectory::operator-=(UnitCell frac) {
      m_from_uccoord -= frac;
      m_to_uccoord -= frac;
      return *this;
    }

    std::unique_ptr<AtomTrajectory> AtomTrajectory::clone() const {
      return std::unique_ptr<AtomTrajectory>(this->_clone());
    }

    Configuration &AtomTrajectory::apply_to_impl(Configuration &config) const {
      auto id = config.atom_id(config.linear_index(from_uccoord()), from_atom_index());
      config.set_atom_id(config.linear_index(to_uccoord()), to_atom_index(), id);
      return config;
    }

    Configuration &AtomTrajectory::apply_reverse_to_impl(Configuration &config) const {
      auto id = config.atom_id(config.linear_index(to_uccoord()), to_atom_index());
      config.set_atom_id(config.linear_index(from_uccoord()), from_atom_index(), id);
      return config;
    }

    void AtomTrajectory::apply_sym_impl(const SymOp &op) {
      m_from_uccoord.apply_sym(op);
      m_to_uccoord.apply_sym(op);

      //MOLECULE_SUPPORT: apply permutation to to/from_value & to/from_atom_index
    }

    void AtomTrajectory::reverse_impl() {
      using std::swap;
      swap(m_from_uccoord, m_to_uccoord);
      swap(m_from_value, m_to_value);
      swap(m_from_atom_index, m_to_atom_index);
    }

    AtomTrajectory *AtomTrajectory::_clone() const {
      return new AtomTrajectory(*this);
    }


    // DiffusionTransformation

    DiffusionTransformation::DiffusionTransformation(const PrimType &prim) :
      DoFTransformation(prim) {}


    std::vector<OccupationTransformation> &DiffusionTransformation::occ_transform();

    const std::vector<OccupationTransformation> &DiffusionTransformation::occ_transform() const;


    std::vector<AtomTrajectory> &DiffusionTransformation::atom_traj();

    const std::vector<AtomTrajectory> &DiffusionTransformation::atom_traj() const;


    AtomTrajectory &DiffusionTransformation::operator+=(UnitCell frac) {
      for(auto &t : m_occ_transform) {
        t += frac;
      }

      for(auto &t : m_atom_traj) {
        t += frac;
      }
      return *this;
    }

    AtomTrajectory &DiffusionTransformation::operator-=(UnitCell frac) {
      for(auto &t : m_occ_transform) {
        t -= frac;
      }

      for(auto &t : m_atom_traj) {
        t -= frac;
      }
      return *this;
    }

    std::unique_ptr<DiffusionTransformation> DiffusionTransformation::clone() const {
      return std::unique_ptr<DiffusionTransformation>(this->_clone());
    }

    Configuration &DiffusionTransformation::apply_to_impl(Configuration &config) const {
      for(auto &t : m_occ_transform) {
        t.apply_to(config);
      }

      for(auto &t : m_atom_traj) {
        t.apply_to(config);
      }
    }

    Configuration &DiffusionTransformation::apply_reverse_to_impl(Configuration &config) const {
      for(auto &t : m_occ_transform) {
        t.apply_reverse_to(config);
      }

      for(auto &t : m_atom_traj) {
        t.apply_reverse_to(config);
      }
    }

    DoFTransformation &DiffusionTransformation::apply_sym_impl(const SymOp &op) {
      for(auto &t : m_occ_transform) {
        t.apply_sym(op);
      }

      for(auto &t : m_atom_traj) {
        t.apply_sym(op);
      }
    }

    void DiffusionTransformation::reverse_impl() {
      for(auto &t : m_occ_transform) {
        t.reverse();
      }

      for(auto &t : m_atom_traj) {
        t.reverse();
      }
    }

    DiffusionTransformation *DiffusionTransformation::_clone() const {
      return new DiffusionTransformation(*this);
    }

  }
}

#endif CASM_DoFTransformation
