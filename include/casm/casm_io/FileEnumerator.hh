#ifndef CASM_FileEnumerator
#define CASM_FileEnumerator

#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"


namespace CASM {

  // -- Get lists of files -----------------------------

  /// \brief Lists all files in a CASM project, for use with 'casm files' command
  ///
  /// \ingroup casmIO
  ///
  class FileEnumerator {

  public:

    /// \brief A CASM project file enumerator
    FileEnumerator(
      const PrimClex &_primclex,
      bool _all_settings = false,
      bool _relative = false);

    /// \brief Enumerate all setting independent files
    template<typename OutputIterator>
    OutputIterator basic_files(OutputIterator result);

    /// \brief Enumerate bset files
    template<typename OutputIterator>
    OutputIterator bset_files(OutputIterator result);

    /// \brief Enumerate reference files
    template<typename OutputIterator>
    OutputIterator reference_files(OutputIterator result);

    /// \brief Enumerate eci files
    template<typename OutputIterator>
    OutputIterator eci_files(OutputIterator result);

    /// \brief Enumerate calculation settings files
    template<typename OutputIterator>
    OutputIterator calc_settings_files(OutputIterator result);

    /// \brief Enumerate calculation status files
    template<typename OutputIterator>
    OutputIterator calc_status_files(OutputIterator result);

    /// \brief Enumerate all training data files
    template<typename OutputIterator>
    OutputIterator all_calc_files(OutputIterator result);


  private:

    /// make paths relative to m_primclex.get_path() if m_relative
    fs::path _if_relative(fs::path path);

    /// output path if it exists
    template<typename OutputIterator>
    OutputIterator _if_exists(OutputIterator result, fs::path path);

    /// \brief Get all regular files that exist in directory 'location'
    template<typename OutputIterator>
    OutputIterator _all_that_exist(OutputIterator result, fs::path location);



    const PrimClex &m_primclex;
    DirectoryStructure m_dir;
    ProjectSettings m_set;
    bool m_all_settings;
    bool m_relative;

    std::vector<std::string> m_all_bset;
    std::vector<std::string> m_all_calctype;
    std::vector<std::string> m_all_property;

  };

  /// \brief A CASM project file enumerator
  ///
  /// \param _primclex The PrimClex representing the project to enumerate files
  /// \param _result OutputIterator taking list of fs::path
  /// \param _all_settings If true, output files for all settings combinations.
  /// \param _relative If true, enumerate paths relative to the project root
  ///         directory. If false, use absolute paths.
  ///
  /// \result Iterator to end of region containing the output file paths
  ///
  FileEnumerator::FileEnumerator(
    const PrimClex &_primclex,
    bool _all_settings,
    bool _relative) :

    m_primclex(_primclex),
    m_dir(m_primclex.dir()),
    m_set(m_primclex.settings()),
    m_all_settings(_all_settings),
    m_relative(_relative),
    m_all_bset(m_dir.all_bset()),
    m_all_calctype(m_dir.all_calctype()),
    m_all_property(m_dir.all_property()) {}


  /// make paths relative to CASM project root directory
  inline fs::path FileEnumerator::_if_relative(fs::path path) {
    if(m_relative) {
      auto a = m_primclex.get_path().string().size() + 1;
      auto b = path.string().size();
      return fs::path(path.string().substr(a, b));
    }
    return path;
  }

  /// output path if it exists
  template<typename OutputIterator>
  OutputIterator FileEnumerator::_if_exists(OutputIterator result, fs::path path) {
    if(fs::exists(path)) {
      *result++ = _if_relative(path);
    }
    return result;
  }

  /// \brief Get all regular files that exist in directory 'location'
  template<typename OutputIterator>
  OutputIterator FileEnumerator::_all_that_exist(OutputIterator result, fs::path location) {
    std::vector<fs::path> all;
    std::string dir;

    // get all
    if(!fs::exists(location)) {
      return result;
    }
    fs::directory_iterator it(location);
    fs::directory_iterator end_it;
    for(; it != end_it; ++it) {
      if(fs::is_regular_file(*it)) {
        *result++ = _if_relative(it->path());
      }
    }
    return result;
  }

  /// \brief Enumerate all setting independent files
  ///
  /// - prim.json
  /// - PRIM
  /// - project_settings.json
  /// - config_list.json
  /// - enumerator plugins
  /// - SCEL
  /// - lattice_point_group.json
  /// - factor_group.json
  /// - crystal_point_group.json
  template<typename OutputIterator>
  OutputIterator FileEnumerator::basic_files(OutputIterator result) {
    std::vector<fs::path> v {
      m_dir.prim(), m_dir.PRIM(),
      m_dir.project_settings(), m_dir.config_list(), m_dir.SCEL(),
      m_dir.lattice_point_group(), m_dir.factor_group(), m_dir.crystal_point_group()
    };
    for(auto it = v.begin(); it != v.end(); ++it) {
      result = _if_exists(result, *it);
    }
    result = _all_that_exist(result, m_dir.enumerator_plugins());
    return result;
  }

  /// \brief Enumerate bset files
  ///
  /// - bspecs.json
  /// - clust.json
  /// - basis.json
  /// - X_Clexulator.cc
  template<typename OutputIterator>
  OutputIterator FileEnumerator::bset_files(OutputIterator result) {

    // bset dependent:
    for(auto bset : m_all_bset) {
      if(!m_all_settings && bset != m_set.default_clex().bset) {
        continue;
      }
      result = _if_exists(result, m_dir.bspecs(bset));
      result = _if_exists(result, m_dir.clust(bset));
      result = _if_exists(result, m_dir.basis(bset));
      result = _if_exists(result, m_dir.clexulator_src(m_set.name(), bset));
    }
    return result;
  }

  /// \brief Enumerate reference files
  ///
  /// - composition_axes.json
  /// - chemical_reference.json
  template<typename OutputIterator>
  OutputIterator FileEnumerator::reference_files(OutputIterator result) {

    // calctype / ref dependent
    for(auto calctype : m_all_calctype) {
      if(!m_all_settings && calctype != m_set.default_clex().calctype) {
        continue;
      }
      auto all_ref = m_dir.all_ref(calctype);
      for(auto ref : all_ref) {
        if(!m_all_settings && ref != m_set.default_clex().ref) {
          continue;
        }
        result = _if_exists(result, m_dir.composition_axes());
        result = _if_exists(result, m_dir.chemical_reference(calctype, ref));
      }
    }
    return result;
  }

  /// \brief Enumerate eci files
  ///
  /// - eci.json
  /// - eci.out (deprecated)
  template<typename OutputIterator>
  OutputIterator FileEnumerator::eci_files(OutputIterator result) {

    // eci
    if(!m_all_settings) {
      result = _if_exists(result, m_dir.eci(m_set.default_clex().property, m_set.default_clex().calctype, m_set.default_clex().ref, m_set.default_clex().bset, m_set.default_clex().eci));
      result = _if_exists(result, m_dir.eci_out(m_set.default_clex().property, m_set.default_clex().calctype, m_set.default_clex().ref, m_set.default_clex().bset, m_set.default_clex().eci));
    }
    else {

      for(auto clex : m_all_property) {
        for(auto calctype : m_all_calctype) {
          auto all_ref = m_dir.all_ref(calctype);
          for(auto ref : all_ref) {
            for(auto bset : m_all_bset) {
              auto all_eci = m_dir.all_eci(clex, calctype, ref, bset);
              for(auto eci : all_eci) {
                result = _if_exists(result, m_dir.eci(clex, calctype, ref, bset, eci));
                result = _if_exists(result, m_dir.eci_out(clex, calctype, ref, bset, eci));
              }
            }
          }
        }
      }

    }
    return result;
  }

  /// \brief Enumerate calculation settings files
  ///
  /// - training_data/settings files (also at supercell and config level)
  template<typename OutputIterator>
  OutputIterator FileEnumerator::calc_settings_files(OutputIterator result) {

    for(auto calctype : m_all_calctype) {
      if(!m_all_settings && calctype != m_set.default_clex().calctype) {
        continue;
      }

      auto scel_begin = m_primclex.get_supercell_list().begin();
      auto scel_end = m_primclex.get_supercell_list().end();

      // calculation settings: global
      result = _all_that_exist(result, m_dir.calc_settings_dir(calctype));

      // supercell level
      for(auto scel_it = scel_begin; scel_it != scel_end; ++scel_it) {
        result = _all_that_exist(result, m_dir.supercell_calc_settings_dir(scel_it->get_name(), calctype));
      }

      auto begin = m_primclex.config_begin();
      auto end = m_primclex.config_end();

      // configuration level
      for(auto it = begin; it != end; ++it) {
        result = _all_that_exist(result, m_dir.configuration_calc_settings_dir(it->name(), calctype));
      }
    }
    return result;
  }

  /// \brief Enumerate calculation status files
  ///
  /// - properties.calc.json and status.json files
  template<typename OutputIterator>
  OutputIterator FileEnumerator::calc_status_files(OutputIterator result) {

    for(auto calctype : m_all_calctype) {
      if(!m_all_settings && calctype != m_set.default_clex().calctype) {
        continue;
      }

      auto begin = m_primclex.config_begin();
      auto end = m_primclex.config_end();

      // calculation summaries: properties.calc.json, status.json
      for(auto it = begin; it != end; ++it) {
        result = _if_exists(result, m_dir.calculated_properties(it->name(), calctype));
        result = _if_exists(result, m_dir.calc_status(it->name(), calctype));
      }

    }
    return result;
  }

  /// \brief Enumerate all training data files
  ///
  /// - all files in 'training_data' directory, recursively
  template<typename OutputIterator>
  OutputIterator FileEnumerator::all_calc_files(OutputIterator result) {

    // get all files in 'training_data', recursively
    fs::recursive_directory_iterator it(m_dir.training_data()), end;
    std::string pattern = "calctype.";
    for(; it != end; ++it) {

      // if not requesting all settings
      if(!m_all_settings) {
        // avoid recursing into other calctypes
        if(fs::is_directory(*it)) {
          std::string dir = it->path().filename().string();
          if(dir.substr(0, pattern.size()) == pattern) {
            if(dir.substr(pattern.size(), dir.size()) != m_set.default_clex().calctype) {
              it.no_push();
            }
          }
        }
      }

      if(fs::is_regular_file(*it)) {
        *result++ = _if_relative(it->path());
      }
    }
    return result;
  }

}

#endif
