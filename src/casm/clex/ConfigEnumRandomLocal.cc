#include "casm/clex/ConfigEnumRandomLocal.hh"
#include "casm/app/enum/EnumInterface_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/clex/FilteredConfigIterator.hh"


extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumRandomLocal_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumRandomLocal>();
  }
}

namespace CASM {

  const std::string ConfigEnumRandomLocal::enumerator_name = "ConfigEnumRandomLocal";

  std::string ConfigEnumRandomLocal::interface_help() {
    return
      "ConfigEnumRandomLocal: \n\n"

      "  confignames: Array of strings (optional) \n"
      "    Names of configurations to be used as initial state of enumeration. All \n"
      "    specified sublattices or sites will be enumerated on and all other DoFs will\n"
      "    maintain the values of the initial state.\n"
      "    Ex: \"confignames\" : [\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]\n\n"

      "  scelnames: Array of strings (optional) \n"
      "    Names of supercells used as initial state of enumeration. All site occupants\n"
      "    will be set to the first listed occupant, and all DoFs will be set to zero.\n"
      "    Ex: \"scelnames\" : [\"SCEL1_1_1_1_0_0_0\",\"SCEL2_2_1_1_0_0_0\"]\n\n"

      "  sublats: array of integers (optional, default none) \n"
      "    Restricts enumeration to specified sublattices. Each sublattice index corresponds\n"
      "    to a basis site in prim.json, indexed from 0.\n"
      "    Ex: \"sublats\" : [0,2]\n\n"

      "  sites: array of 4-entry integer arrays (optional, default none) \n"
      "    Restricts enumeration to specified sites. Sites are specified in [b,i,j,k] convention,\n"
      "    where 'b' is sublattice index and [i,j,k] specifies linear combinations of primitive-\n"
      "    cell lattice vectors.\n"
      "    Ex: \"sites\" : [[0,0,0,0],\n"
      "                   [2,0,0,0]]\n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
      "    Indicate supercells to use as initial states of enumeration in terms of size\n"
      "    and unit cell via a JSON object conforming to the format of 'ScelEnum' JSON\n"
      "    settings. \"scelnames\" will override \"supercells\", but if neither is specified\n"
      "    all existing supercells are used by default. See 'ScelEnum' description for details.\n\n"

      "  n_config: integer (optional, default=100) \n"
      "    How many random configurations to generate. Includes duplicate and pre-\n"
      "    existing configurations.                                                 \n\n"

      "  magnitude: positive number (optional, default=1) \n"
      "    Magnitude used to scale random vector at each site. If \"distribution\" == \"normal\",\n"
      "    magnitude specifies standard deviation of D-dimensional Gaussian (D is dimension of \n"
      "    site DoF value). If \"distribution\" == \"uniform\", magnitude is radius of D-dimensional\n"
      "    ball from which random vectors are chosen.\n\n"

      "  dof: string (required) \n"
      "    Name of site degree of freecom  for which normal coordinates are to be generated.\n"
      "    Must be one of the degrees of freedom under consideration in the current project,\n"
      "    as determined by prim.json\n\n"

      "  distribution: string (optional, default=\"normal\") \n"
      "    Distribution from which perturbation vectors are drawn. Options are \"uniform\"\n"
      "    (i.e., uniform on the unit sphere), or \"normal\" (i.e., zero-centered Gaussian).\n\n"

      "  primitive_only: bool (default=true)\n"
      "    If true, only the primitive form of a configuration is saved in the      \n"
      "    configuration list. Otherwise, both primitive and non-primitive          \n"
      "    configurations are saved. \n\n"

      "  Examples:\n"
      "    To enumerate 200 random occupations in supercells up to and including size 4:\n"
      "      casm enum --method ConfigEnumRandomLocal -i \n"
      "        '{\"supercell\":{\"max\":4}, \"n_config\": 200}' \n"
      "\n"
      "    To enumerate 200 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomLocal -i '{\"n_config\": 200}' \n"
      "\n"
      "    To enumerate 100 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomLocal' \n"
      "\n"
      "    To enumerate 200 random occupations in particular supercells:\n"
      "      casm enum --method ConfigEnumRandomLocal -i \n"
      "      '{ \n"
      "        \"supercells\": { \n"
      "          \"name\": [\n"
      "            \"SCEL1_1_1_1_0_0_0\",\n"
      "            \"SCEL2_1_2_1_0_0_0\",\n"
      "            \"SCEL4_1_4_1_0_0_0\"\n"
      "          ]\n"
      "        }, \n"
      "        \"n_config\": 200\n"
      "      }' \n\n";
  }

  int ConfigEnumRandomLocal::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt,
    EnumeratorMap const *interface_map) {

    std::vector<ConfigEnumInput> in_configs = make_enumerator_input_configs(primclex, _kwargs, enum_opt, interface_map);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);
    MTRand mtrand;

    DoFKey dof_key;
    if(_kwargs.contains("dof"))
      _kwargs["dof"].get(dof_key);
    else
      throw std::runtime_error("Keyword argument \"dof\" must be provided to run enumeration method \"ConfigEnumRandomLocal\".");

    double mag;
    if(_kwargs.contains("magnitude"))
      _kwargs["magnitude"].get(mag);
    else
      throw std::runtime_error("Keyword argument \"magnitude\" must be provided to run enumeration method \"ConfigEnumRandomLocal\".");

    bool normal = true;
    if(_kwargs.contains("distribution")) {
      std::string distribution;
      _kwargs["distribution"].get(distribution);
      if(distribution == "uniform")
        normal = false;
      else if(distribution != "normal")
        throw std::runtime_error("Keyword argument \"distribution\"ethod \"ConfigEnumRandomLocal\".");
    }

    Index n_config;
    _kwargs.get_else<Index>(n_config, "n_config", 100);

    bool primitive_only = true;
    _kwargs.get_if(primitive_only, "primitive_only");

    auto lambda = [&](const ConfigEnumInput & _input) {
      return notstd::make_unique<ConfigEnumRandomLocal>(_input, dof_key, n_config, mag, normal, mtrand);
    };

    int returncode = insert_configs(
                       enumerator_name,
                       primclex,
                       in_configs.begin(),
                       in_configs.end(),
                       lambda,
                       filter_expr,
                       primitive_only,
                       CASM::dry_run(_kwargs, enum_opt));

    return returncode;
  }

  /// \brief Constructor
  ///
  /// \param _initial,_final Initial and final configurations to interpolate between
  /// \param _size The total number of configurations to enumerate, including
  ///              the initial and final configurations
  ///
  /// - The `final` configuration is *not* pointed at by the end iterator,
  ///   which points past-the-final element, as is typical
  /// - `_size` will be equal to \code std::distance(this->begin(), this->end()) \endcode
  ConfigEnumRandomLocal::ConfigEnumRandomLocal(ConfigEnumInput const &_in_config,
                                               DoFKey const &_dof_key,
                                               Index _n_config,
                                               double _mag,
                                               bool _normal,
                                               MTRand &_mtrand):
    m_n_config(_n_config),
    m_mtrand(&_mtrand),
    m_mag(_mag),
    m_normal(_normal),
    m_unit_length(DoF::BasicTraits(_dof_key).unit_length()),
    m_site_selection(_in_config.sites().begin(), _in_config.sites().end()) {

    if(m_unit_length)
      m_normal = false;

    if(m_n_config < 0) {
      throw std::runtime_error("Error in ConfigEnumRandomLocal: n_config < 0");
    }
    if(m_n_config == 0) {
      this->_invalidate();
      return;
    }



    m_current = notstd::make_cloneable<Configuration>(_in_config.config());

    reset_properties(*m_current);
    this->_initialize(&(*m_current));

    m_dof_vals = &(m_current->configdof().local_dof(_dof_key));

    auto const &dof_info = m_dof_vals->info();
    for(Index l : m_site_selection)
      m_dof_dims.push_back(dof_info[m_current->sublat(l)].dim());

    // Make initial random config
    this->randomize();
    _set_step(0);
    m_current->set_source(this->source(step()));

    //std::cout << "Selection: " << m_site_selection <<"\n"
    //        << "dofkey: " << _dof_key << "\n"
    //        << "mag: " << m_mag << "\n"
    //        << "dof_dims: " << m_dof_dims << "\n"
    //        << "unit_length: " << m_unit_length << "\n"
    //        << "m_normal: " << m_normal << "\n";
  }

  /// Set m_current to correct value at specified step and return a reference to it
  void ConfigEnumRandomLocal::increment() {

    this->_increment_step();
    if(step() < m_n_config) {
      this->randomize();
      m_current->set_source(this->source(step()));
    }
    else {
      this->_invalidate();
    }

  }

  void ConfigEnumRandomLocal::randomize() {
    double tnorm;
    for(Index i = 0; i < m_site_selection.size(); ++i) {
      for(Index j = 0; j < m_dof_dims[i]; ++j) {
        m_dof_vals->site_value(m_site_selection[i])[j] = m_mtrand->randNorm(0., m_mag);
      }
      if(!m_normal) {
        tnorm = m_dof_vals->site_value(m_site_selection[i]).norm();
        if(!m_unit_length) {
          tnorm += -m_mag * std::log(m_mtrand->rand());
        }
        (m_dof_vals->site_value(m_site_selection[i])) /= tnorm;
      }
    }
    //std::cout << "Randomized: \n" << m_dof_vals->values() << "\n";
  }

}
