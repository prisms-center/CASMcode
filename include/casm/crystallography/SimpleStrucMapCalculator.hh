#ifndef CASM_SimpleStrucMapCalculator
#define CASM_SimpleStrucMapCalculator

#include "casm/crystallography/StrucMapCalculatorInterface.hh"
#include "casm/crystallography/Adapter.hh"

namespace CASM {
  namespace xtal {
    class SimpleStructure;
    struct SymOp;
    typedef std::vector<SymOp> SymOpVector;

    // In this file:
    struct MappingNode;
    class StrucMapCalculatorInterface;
    class SimpleStrucMapCalculator;

    class SimpleStrucMapCalculator : public StrucMapCalculatorInterface {
    public:
      SimpleStrucMapCalculator(SimpleStructure _parent,
                               SymOpVector _point_group = {SymOp::identity()},
                               SimpleStructure::SpeciesMode species_mode = SimpleStructure::SpeciesMode::ATOM,
                               StrucMapping::AllowedSpecies allowed_species = {}) :
        StrucMapCalculatorInterface(std::move(_parent),
                                    std::move(_point_group),
                                    species_mode,
                                    std::move(allowed_species)) {

        for(Index i = 0; i < _allowed_species().size(); ++i) {
          for(std::string const &sp : _allowed_species()[i]) {
            if(sp == "Va" || sp == ("VA") || sp == "va") {
              m_va_allowed.insert(i);
              ++m_max_n_species["Va"];
            }
            else {
              ++m_max_n_species[sp];
            }
          }
        }
        _max_n_va() = m_va_allowed.size();
      }

      template <typename ExternSymOpVector>
      SimpleStrucMapCalculator(SimpleStructure _parent,
                               ExternSymOpVector _point_group = {SymOp::identity()},
                               SimpleStructure::SpeciesMode species_mode = SimpleStructure::SpeciesMode::ATOM,
                               StrucMapping::AllowedSpecies allowed_species = {}) :
        SimpleStrucMapCalculator(_parent, adapter::Adapter<SymOpVector, ExternSymOpVector>()(_point_group), species_mode, allowed_species) {
      }

      virtual ~SimpleStrucMapCalculator() {}

      std::map<std::string, Index> const &max_n_species() const {
        return m_max_n_species;
      }

      std::vector<Eigen::Vector3d> translations(MappingNode const &_node,
                                                SimpleStructure const &child_struc) const override;

      /// \brief Creates copy of _child_struc by applying isometry, lattice transformation, translation, and site permutation of _node
      virtual SimpleStructure resolve_setting(MappingNode const &_node,
                                              SimpleStructure const &_child_struc) const override;

      void finalize(MappingNode &_node,
                    SimpleStructure const &child_struc) const override;


      bool populate_cost_mat(MappingNode &_node,
                             SimpleStructure const &child_struc) const override;

      void populate_displacement(MappingNode &_node,
                                 SimpleStructure const &child_struc) const;

    private:
      bool _sublat_allows_va(Index b) const {
        return m_va_allowed.count(b);
      }

      /// \brief Make an exact copy of the calculator (including any initialized members)
      virtual StrucMapCalculatorInterface *_clone() const override {
        return new SimpleStrucMapCalculator(*this);
      }

      /// \brief Make an exact copy of the calculator (including any initialized members)
      virtual StrucMapCalculatorInterface *_quasi_clone(SimpleStructure _parent,
                                                        SymOpVector _point_group = {SymOp::identity()},
                                                        SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                        StrucMapping::AllowedSpecies _allowed_species = {}) const override {
        return new SimpleStrucMapCalculator(std::move(_parent), std::move(_point_group), _species_mode, std::move(_allowed_species));
      }

      std::unordered_set<Index> m_va_allowed;

      /// \brief maximum allowed number of each species
      std::map<std::string, Index> m_max_n_species;
    };
  }
}
#endif
