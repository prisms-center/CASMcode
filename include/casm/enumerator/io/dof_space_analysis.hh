#ifndef CASM_enumerator_io_json_dof_space_analysis
#define CASM_enumerator_io_json_dof_space_analysis

#include "casm/global/definitions.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

  class ConfigEnumInput;
  class DirectoryStructure;
  class DoFSpace;
  class SymGroup;
  struct VectorSpaceSymReport;
  class jsonParser;

  namespace DoFSpaceIO {

    class OutputImpl {
    public:

      /// Write symmetry groups (lattice point group, factor_group, crystal_point_group)
      void write_symmetry(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        std::vector<PermuteIterator> const &group);

      /// Write symmetry groups (lattice point group, factor_group, crystal_point_group)
      ///
      /// - lattice_point_group: point group of supercell lattice
      /// - factor_group: subgroup of supercell factor group that keeps configuration DoF invariant
      ///   and selected sites invariant
      /// - crystal_point_group: factor_group excluding translations
      virtual void write_symmetry(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        SymGroup const &lattice_point_group,
        SymGroup const &factor_group,
        SymGroup const &crystal_point_group) = 0;

      // /// Write input state configuration
      // virtual void write_config(
      //   Index state_index,
      //   std::string const &identifier,
      //   ConfigEnumInput const &config_enum_input) = 0;

      /// Write input state structure
      virtual void write_structure(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) = 0;

      /// Write dof space analysis
      virtual void write_dof_space(
        Index state_index,
        DoFSpace const &dof_space,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        std::optional<VectorSpaceSymReport> const &sym_report) = 0;

    };

    class DirectoryOutput : public OutputImpl {
    public:

      using OutputImpl::write_symmetry;

      void write_symmetry(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        SymGroup const &lattice_point_group,
        SymGroup const &factor_group,
        SymGroup const &crystal_point_group) override;

      // void write_config(
      //   Index state_index,
      //   std::string const &identifier,
      //   ConfigEnumInput const &config_enum_input) override;

      void write_structure(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) override;

      void write_dof_space(
        Index state_index,
        DoFSpace const &dof_space,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        std::optional<VectorSpaceSymReport> const &sym_report) override;

    private:

      /// For SymmetryDirectoryOutput, configurations must exist in database and this will throw otherwise
      virtual void _check_config(Index state_index,
                                 std::string const &identifier,
                                 ConfigEnumInput const &config_enum_input) = 0;

      /// For SymmetryDirectoryOutput vs SequentialDirectoryOutput return output directory path
      virtual fs::path _output_dir(Index state_index,
                                   std::string const &identifier,
                                   ConfigEnumInput const &config_enum_input) = 0;

    };

    /// Implementation that outputs to <casm_project>/symmetry/analysis/<configname>
    class SymmetryDirectoryOutput : public DirectoryOutput {
    public:
      SymmetryDirectoryOutput(DirectoryStructure const &dir);

    private:

      /// For SymmetryDirectoryOutput, configurations must exist in database and this will throw otherwise
      void _check_config(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) override;

      /// For SymmetryDirectoryOutput return output directory path in symmetry dir
      fs::path _output_dir(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) override;

      DirectoryStructure const &m_dir;
    };

    /// Implementation that outputs to <output_dir>/dof_space/state.<index>
    class SequentialDirectoryOutput : public DirectoryOutput {
    public:
      SequentialDirectoryOutput(fs::path output_dir);

    private:

      /// For SequentialDirectoryOutput, any input state is allowed
      void _check_config(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) override;

      /// For SequentialDirectoryOutput return output directory path
      fs::path _output_dir(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) override;

      fs::path m_output_dir;

    };

    /// Implementation that outputs all results to one JSON file at <output_dir>/dof_space.json
    class CombinedJsonOutput : public OutputImpl {
    public:

      CombinedJsonOutput(fs::path output_dir);

      ~CombinedJsonOutput();

      using OutputImpl::write_symmetry;

      void write_symmetry(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        SymGroup const &lattice_point_group,
        SymGroup const &factor_group,
        SymGroup const &crystal_point_group) override;

      // void write_config(
      //   Index state_index,
      //   std::string const &identifier,
      //   ConfigEnumInput const &config_enum_input) override;

      void write_structure(
        Index state_index,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input) override;

      void write_dof_space(
        Index state_index,
        DoFSpace const &dof_space,
        std::string const &identifier,
        ConfigEnumInput const &config_enum_input,
        std::optional<VectorSpaceSymReport> const &sym_report) override;

    private:

      jsonParser &_output_json(Index state_index);

      jsonParser m_combined_json;

      fs::path m_output_dir;

    };


    struct DoFSpaceAnalysisOptions {

      std::vector<std::string> dofs;
      bool sym_axes = true;
      bool calc_wedge = false;
      bool write_symmetry = true;
      bool write_structure = true;
    };

    void output_dof_space(
      Index state_index,
      std::string const &identifier,
      ConfigEnumInput const &input_state,
      DoFSpaceAnalysisOptions const &options,
      OutputImpl &output);

    void dof_space_analysis(
      std::vector<std::pair<std::string, ConfigEnumInput>> const &named_inputs,
      DoFSpaceAnalysisOptions const &options,
      OutputImpl &output);
  }

}

#endif
