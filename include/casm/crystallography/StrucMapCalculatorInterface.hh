#ifndef CASM_StrucMapCalculatorInterface
#define CASM_StrucMapCalculatorInterface

#include <vector>
#include <iostream>
#include <unordered_set>
#include "casm/global/definitions.hh"
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
      using AllowedSpecies = std::vector<std::vector<std::string>>;

    }

    class StrucMapCalculatorInterface {
    public:

      StrucMapCalculatorInterface(SimpleStructure _parent,
                                  SymOpVector _point_group = {SymOp::identity()},
                                  SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                  StrucMapping::AllowedSpecies allowed_species = {}) :
        m_parent(std::move(_parent)),
        m_point_group(std::move(_point_group)),
        m_species_mode(_species_mode) {

        //std::cout << "allowed_species:\n";
        if(allowed_species.empty()) {
          auto const &p_info(this->struc_info(parent()));
          allowed_species.resize(p_info.size());
          for(Index i = 0; i < p_info.size(); ++i) {
            allowed_species[i].push_back(p_info.names[i]);
            //std::cout << *(m_allowed_species[i].begin()) << "  ";
          }
          //std::cout << "\n";
        }
        set_allowed_species(allowed_species);
      }

      SimpleStructure::Info const &struc_info(SimpleStructure const &_struc) const {
        return _struc.info(m_species_mode);
      }

      SimpleStructure::Info &struc_info(SimpleStructure &_struc) const {
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

      void set_allowed_species(StrucMapping::AllowedSpecies allowed_species) {
        m_allowed_species = std::move(allowed_species);

        //Analyze allowed_species:
        for(Index i = 0; i < m_allowed_species.size(); ++i) {
          for(std::string &sp : m_allowed_species[i]) {
            if(sp == "Va" || sp == ("VA") || sp == "va") {
              //Is vacancy
              sp = "Va";
              m_va_allowed.insert(i);
            }
            else if(m_allowed_species[i].size() > 1) {
              //Not vacancy, but variable species
              m_fixed_species[sp] = 0;
            }
            else {
              //potentially fixed species
              auto it = m_fixed_species.find(sp);
              if(it == m_fixed_species.end())
                m_fixed_species[sp] = 1;
              else if((it->second) > 0)
                ++(it->second);
            }
            if(!m_max_n_species.count(sp))
              m_max_n_species[sp] = 0;
            ++m_max_n_species[sp];
          }
        }

        //std::cout << "Fixed_species: ";
        for(auto it = m_fixed_species.begin(); it != m_fixed_species.end();) {
          auto curr = it;
          ++it;
          if((curr->second) == 0)
            m_fixed_species.erase(curr);
          //std::cout << it->first << "  ";
        }

      }

      std::unordered_set<Index> const &va_allowed() const {
        return m_va_allowed;
      }

      StrucMapping::FixedSpecies const &fixed_species() const {
        return m_fixed_species;
      }

      /// \brief Return maximum possible number of vacancies in underlying primitive structure
      Index max_n_va() const {
        return va_allowed().size();
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

      StrucMapping::AllowedSpecies const &_allowed_species() const {
        return m_allowed_species;
      }

      StrucMapping::FixedSpecies const &_fixed_species() const {
        return m_fixed_species;
      }

      bool _sublat_allows_va(Index b) const {
        return m_va_allowed.count(b);
      }

      /// \brief maximum allowed number of each species
      std::map<std::string, Index> const &_max_n_species() const {
        return m_max_n_species;
      }

    private:
      SimpleStructure m_parent;

      SymOpVector m_point_group;

      SimpleStructure::SpeciesMode m_species_mode;

      StrucMapping::AllowedSpecies m_allowed_species;

      StrucMapping::FixedSpecies m_fixed_species;

      /// \brief maximum allowed number of each species
      std::map<std::string, Index> m_max_n_species;

      std::unordered_set<Index> m_va_allowed;

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
