#ifndef CASM_StrucMapCalculatorInterface
#define CASM_StrucMapCalculatorInterface

#include <vector>
#include <unordered_set>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/Adapter.hh"

namespace CASM {
  namespace xtal {
    class SimpleStructure;

    // In this file:
    struct MappingNode;
    class StrucMapCalculatorInterface;

    namespace StrucMapping {

      // Denotes (name, number composition) of any species whose number-composition is fixed for the parent primitive cell
      using FixedSpecies = std::map<std::string, Index>;

      // List of species allowed at each site of primitive
      using AllowedSpecies = std::vector<std::unordered_set<std::string>>;

    }

    class StrucMapCalculatorInterface {
    public:

      StrucMapCalculatorInterface(SimpleStructure _parent,
                                  SymOpVector _point_group = {SymOp::identity()},
                                  SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                  StrucMapping::AllowedSpecies allowed_species = {}) :
        m_parent(std::move(_parent)),
        m_point_group(std::move(_point_group)),
        m_species_mode(_species_mode),
        m_allowed_species(std::move(allowed_species)) {

        if(m_allowed_species.empty()) {
          auto const &p_info(info(parent()));
          m_allowed_species.resize(p_info.size());
          for(Index i = 0; i < m_allowed_species.size(); ++i)
            m_allowed_species[i].insert(p_info.names[i]);
        }


        for(auto const &slist : m_allowed_species) {
          if(slist.size() != 1) {
            for(auto const &sp : slist) {
              _fixed_species()[sp] = 0;
            }
          }
          auto it = _fixed_species().find(*slist.begin());
          if(it == _fixed_species().end())
            _fixed_species()[*slist.begin()] = 1;
          else if((it->second) > 0)
            ++(it->second);
        }
        for(auto it = _fixed_species().begin(); it != _fixed_species().end(); ++it) {
          if((it->second) == 0)
            _fixed_species().erase(it);
        }
      }

      SimpleStructure::Info const &info(SimpleStructure const &_struc) const {
        return _struc.info(m_species_mode);
      }

      SimpleStructure const &parent() const {
        return m_parent;
      }

      SymOpVector const &point_group() const {
        return m_point_group;
      }

      void set_point_group(SymOpVector _point_group) {
        m_point_group = std::move(_point_group);
      }

      template <typename ExternSymOpVector>
      void set_point_group(ExternSymOpVector _point_group) {
        return this->set_point_group(adapter::Adapter<SymOpVector, ExternSymOpVector>()(_point_group));
      }

      StrucMapping::FixedSpecies const &fixed_species() const {
        return m_fixed_species;
      }

      /// \brief Return maximum possible number of vacancies in underlying primitive structure
      Index max_n_va() const {
        return m_max_n_va;
      }

      virtual ~StrucMapCalculatorInterface() {}

      /// \brief construct list of prospective mapping translations
      virtual std::vector<Eigen::Vector3d> translations(MappingNode const &_node,
                                                        SimpleStructure const &child_struc) const = 0;

      /// \brief Creates copy of _child_struc by applying isometry, lattice transformation, translation, and site permutation of _node
      virtual SimpleStructure resolve_setting(MappingNode const &_node,
                                              SimpleStructure const &_child_struc) const = 0;

      /// \brief Calculates final mapping score and sets _node.is_valid
      virtual void finalize(MappingNode &_node,
                            SimpleStructure const &child_struc) const = 0;

      virtual bool populate_cost_mat(MappingNode &_node,
                                     SimpleStructure const &child_struc) const = 0;

      /// \brief Make an exact copy of the calculator (including any initialized members)
      std::unique_ptr<StrucMapCalculatorInterface> clone() const {
        return std::unique_ptr<StrucMapCalculatorInterface>(this->_clone());
      }

      /// \brief Make an exact copy of the calculator (including any initialized members)
      std::unique_ptr<StrucMapCalculatorInterface> quasi_clone(SimpleStructure _parent,
                                                               SymOpVector _point_group = {SymOp::identity()},
                                                               SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                               StrucMapping::AllowedSpecies _allowed_species = {}) const {
        return std::unique_ptr<StrucMapCalculatorInterface>(this->_quasi_clone(std::move(_parent),
                                                                               std::move(_point_group),
                                                                               _species_mode,
                                                                               std::move(_allowed_species)));
      }

      template <typename ExternSymOpVector>
      std::unique_ptr<StrucMapCalculatorInterface> quasi_clone(SimpleStructure _parent,
                                                               ExternSymOpVector _point_group = {SymOp::identity()},
                                                               SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                               StrucMapping::AllowedSpecies _allowed_species = {}) const {
        return this->quasi_clone(_parent, adapter::Adapter<SymOpVector, ExternSymOpVector>()(_point_group), _species_mode, _allowed_species);
      }

    protected:
      StrucMapping::AllowedSpecies &_allowed_species() {
        return m_allowed_species;
      }

      StrucMapping::AllowedSpecies const &_allowed_species() const {
        return m_allowed_species;
      }

      StrucMapping::FixedSpecies &_fixed_species() {
        return m_fixed_species;
      }

      Index &_max_n_va() {
        return m_max_n_va;
      }
    private:
      SimpleStructure m_parent;

      SymOpVector m_point_group;

      SimpleStructure::SpeciesMode m_species_mode;

      StrucMapping::AllowedSpecies m_allowed_species;

      StrucMapping::FixedSpecies m_fixed_species;

      Index m_max_n_va;

      /// \brief Make an exact copy of the calculator (including any initialized members)
      virtual StrucMapCalculatorInterface *_clone() const = 0;

      /// \brief Make an exact copy of the calculator (including any initialized members)
      virtual StrucMapCalculatorInterface *_quasi_clone(SimpleStructure _parent,
                                                        SymOpVector _point_group = {SymOp::identity()},
                                                        SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                        StrucMapping::AllowedSpecies _allowed_species = {}) const = 0;

    };
  }
}

#endif
