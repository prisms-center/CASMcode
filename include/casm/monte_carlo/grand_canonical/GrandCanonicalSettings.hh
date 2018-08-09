#ifndef CASM_GrandCanonicalSettings
#define CASM_GrandCanonicalSettings

#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/external/boost.hh"

namespace CASM {

  class GrandCanonicalConditions;

  class GrandCanonicalSettings : public EquilibriumMonteSettings {

  public:

    /// \brief Default constructor
    GrandCanonicalSettings() {}

    /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
    GrandCanonicalSettings(const PrimClex &primclex, const fs::path &read_path);


    // --- GrandCanonicalConditions settings ---------------------

    /// \brief Expects initial_conditions
    GrandCanonicalConditions initial_conditions() const;

    /// \brief Expects final_conditions
    GrandCanonicalConditions final_conditions() const;

    /// \brief Expects incremental_conditions
    GrandCanonicalConditions incremental_conditions() const;

    /// \brief Expects custom_conditions
    std::vector<GrandCanonicalConditions> custom_conditions() const;


    // --- Project settings ---------------------

    /// \brief Get formation energy cluster expansion
    ClexDescription formation_energy(const PrimClex &primclex) const;


    // --- Sampler settings ---------------------

    /// \brief Construct MonteSamplers as specified in the MonteSettings
    template<typename SamplerInsertIterator>
    SamplerInsertIterator samplers(const PrimClex &primclex, SamplerInsertIterator result) const;

    /// \brief Return true if all correlations should be sampled
    bool all_correlations() const;


  private:

    GrandCanonicalConditions _conditions(std::string name) const;
    GrandCanonicalConditions _conditions(const jsonParser &json) const;

    template<typename jsonParserIteratorType>
    std::tuple<bool, double> _get_precision(jsonParserIteratorType it) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_comp_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_comp_n_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_site_frac_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_atom_frac_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_all_correlations_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_non_zero_eci_correlations_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_query_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;

  };


  /// \brief Construct MonteSamplers as specified in the MonteSettings
  ///
  /// The requested MonteSamplers are inserted in 'result' as
  ///   std::pair<std::string, notstd::cloneable_ptr<MonteSampler> >
  ///
  template<typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::samplers(
    const PrimClex &primclex,
    SamplerInsertIterator result) const {

    if(method() == Monte::METHOD::LTE1) {//hack
      return result;
    }

    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler *ptr;

    // copy so we can add required measurements
    std::string level1 = "data";
    std::string level2 = "measurements";
    jsonParser t_measurements = (*this)[level1][level2];

    // find existing measurements
    std::set<std::string> input_measurements;
    for(auto it = t_measurements.cbegin(); it != t_measurements.cend(); it++) {
      input_measurements.insert((*it)["quantity"].get<std::string>());
    }

    std::vector<std::string> required = {
      "potential_energy",
      "formation_energy",
      "comp",
      "comp_n"
    };

    // add required if not already requested
    for(auto it = required.begin(); it != required.end(); ++it) {
      if(std::find(input_measurements.begin(), input_measurements.end(), *it) == input_measurements.end()) {
        jsonParser json;
        json["quantity"] = *it;
        t_measurements.push_back(json);
      }
    }



    try {

      for(auto it = t_measurements.cbegin(); it != t_measurements.cend(); it++) {

        prop_name = (*it)["quantity"].get<std::string>();

        // scalar quantities that we incrementally update
        std::vector<std::string> scalar_possible = {
          "formation_energy",
          "potential_energy"
        };

        // check if property found is in list of possible scalar properties
        if(std::find(scalar_possible.cbegin(), scalar_possible.cend(), prop_name) != scalar_possible.cend()) {

          std::tie(must_converge, prec) = _get_precision(it);

          // if 'must converge'
          if(must_converge) {
            ptr = new ScalarMonteSampler(prop_name, prop_name, prec, confidence(), data_maxlength);
          }
          else {
            ptr = new ScalarMonteSampler(prop_name, prop_name, confidence(), data_maxlength);
          }

          *result++ = std::make_pair(prop_name, notstd::cloneable_ptr<MonteSampler>(ptr));
          continue;
        }

        // scalar quantities that we incrementally update
        std::vector<std::string> vector_possible = {
          "comp",
          "comp_n",
          "site_frac",
          "atom_frac",
          "all_correlations",
          "non_zero_eci_correlations"
        };

        // check if property found is in list of possible vector properties
        if(std::find(vector_possible.cbegin(), vector_possible.cend(), prop_name) != vector_possible.cend()) {

          // construct MonteSamplers for 'comp'
          if(prop_name == "comp") {

            result = _make_comp_samplers(primclex, it, result);

          }

          // construct MonteSamplers for 'comp_n'
          else if(prop_name == "comp_n") {

            result = _make_comp_n_samplers(primclex, it, result);

          }

          // construct MonteSamplers for 'site_frac'
          else if(prop_name == "site_frac") {

            result = _make_site_frac_samplers(primclex, it, result);

          }

          // construct MonteSamplers for 'atom_frac'
          else if(prop_name == "atom_frac") {

            result = _make_atom_frac_samplers(primclex, it, result);

          }

          // construct MonteSamplers for 'all_correlations'
          else if(prop_name == "all_correlations") {

            result = _make_all_correlations_samplers(primclex, it, result);

          }

          // construct MonteSamplers for 'non_zero_eci_correlations'
          else if(prop_name == "non_zero_eci_correlations") {

            result = _make_non_zero_eci_correlations_samplers(primclex, it, result);

          }
          continue;
        }

        // custom query
        _make_query_samplers(primclex, it, result);

      }

    }
    catch(std::runtime_error &e) {
      Log &err_log = default_err_log();
      err_log.error<Log::standard>("'MonteSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result)'");
      err_log << "Error reading [\"" << level1 << "\"][\"" << level2 << "\"]\n" << std::endl;
      throw;
    }

    return result;
  }

  template<typename jsonParserIteratorType>
  std::tuple<bool, double> GrandCanonicalSettings::_get_precision(jsonParserIteratorType it) const {
    if(it->contains("precision")) {
      return std::make_tuple(true, (*it)["precision"]. template get<double>());
    }
    else {
      return std::make_tuple(false, 0.0);
    }
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_comp_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler *ptr;

    for(size_type i = 0; i < primclex.composition_axes().independent_compositions(); i++) {

      print_name = std::string("comp(") + std::string(1, (char)(i + ((int) 'a'))) + ")";

      std::tie(must_converge, prec) = _get_precision(it);

      // if 'must converge'
      if(must_converge) {
        ptr = new CompMonteSampler(i, primclex.composition_axes(), print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new CompMonteSampler(i, primclex.composition_axes(), print_name, confidence(), data_maxlength);
      }

      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

    }

    return result;
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_comp_n_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler *ptr;

    for(size_type i = 0; i < primclex.composition_axes().components().size(); i++) {

      prop_name = "comp_n";

      print_name = std::string("comp_n(") + primclex.composition_axes().components()[i] + ")";

      std::tie(must_converge, prec) = _get_precision(it);

      // if 'must converge'
      if(must_converge) {
        ptr = new VectorMonteSampler(prop_name, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new VectorMonteSampler(prop_name, i, print_name, confidence(), data_maxlength);
      }

      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

    }

    return result;
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_site_frac_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler *ptr;

    for(size_type i = 0; i < primclex.composition_axes().components().size(); i++) {

      // SiteFracMonteSampler uses 'comp_n' to calculate 'site_frac'
      print_name = std::string("site_frac(") + primclex.composition_axes().components()[i] + ")";

      std::tie(must_converge, prec) = _get_precision(it);

      // if 'must converge'
      if(must_converge) {
        ptr = new SiteFracMonteSampler(i, primclex.get_prim().basis.size(), print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new SiteFracMonteSampler(i, primclex.get_prim().basis.size(), print_name, confidence(), data_maxlength);
      }

      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

    }

    return result;
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_atom_frac_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler *ptr;

    // find the index for vacancies, if they are allowed,
    //  if not set to primclex.composition_axes().components().size()

    size_type vacancy_index = primclex.composition_axes().components().size();
    for(size_type i = 0; i < primclex.composition_axes().components().size(); i++) {

      // sample for non-vacancy components
      if(Specie(primclex.composition_axes().components()[i]).is_vacancy()) {
        vacancy_index = i;
        break;
      }
    }


    for(size_type i = 0; i < primclex.composition_axes().components().size(); i++) {

      // sample for non-vacancy components
      if(!Specie(primclex.composition_axes().components()[i]).is_vacancy()) {

        print_name = std::string("atom_frac(") + primclex.composition_axes().components()[i] + ")";

        std::tie(must_converge, prec) = _get_precision(it);

        // if 'must converge'
        if(must_converge) {
          ptr = new AtomFracMonteSampler(i, vacancy_index, print_name, prec, confidence(), data_maxlength);
        }
        else {
          ptr = new AtomFracMonteSampler(i, vacancy_index, print_name, confidence(), data_maxlength);
        }

        *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

      }
    }

    return result;
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_all_correlations_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler *ptr;

    for(size_type i = 0; i < primclex.clexulator(formation_energy(primclex)).corr_size(); i++) {

      prop_name = "corr";
      print_name = std::string("corr(") + std::to_string(i) + ")";

      std::tie(must_converge, prec) = _get_precision(it);

      // if 'must converge'
      if(must_converge) {
        ptr = new VectorMonteSampler(prop_name, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new VectorMonteSampler(prop_name, i, print_name, confidence(), data_maxlength);
      }

      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

    }

    return result;
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_non_zero_eci_correlations_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;

    MonteSampler *ptr;

    ECIContainer _eci = primclex.eci(formation_energy(primclex));

    for(size_type ii = 0; ii < _eci.index().size(); ii++) {

      prop_name = "corr";

      // store non-zero eci index in i
      size_type i = _eci.index()[ii];

      print_name = std::string("corr(") + std::to_string(i) + ")";

      std::tie(must_converge, prec) = _get_precision(it);

      // if 'must converge'
      if(must_converge) {
        ptr = new VectorMonteSampler(prop_name, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new VectorMonteSampler(prop_name, i, print_name, confidence(), data_maxlength);
      }

      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

    }

    return result;
  }

  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::_make_query_samplers(
    const PrimClex &primclex,
    jsonParserIteratorType it,
    SamplerInsertIterator result) const {

    size_type data_maxlength = max_data_length();
    double must_converge;
    double prec;
    std::string prop_name = (*it)["quantity"].template get<std::string>();
    MonteSampler *ptr;

    const auto &dict = primclex.settings().query_handler<Configuration>().dict();

    typedef QueryMonteSampler::Formatter FormatterType;
    std::shared_ptr<FormatterType> formatter = std::make_shared<FormatterType>(
                                                 dict.parse(prop_name));

    // make example config to test:
    Supercell tscel(const_cast<PrimClex *>(&primclex), simulation_cell_matrix());
    Configuration config(tscel);
    config.init_occupation();
    Eigen::VectorXd test = formatter->get().evaluate_as_matrix(config).row(0);
    auto col = formatter->get().col_header(config);

    if(test.size() != col.size()) {
      std::stringstream ss;
      ss << "Error constructing Monte Carlo samplers from query: '" << prop_name << "'";
      primclex.err_log() << ss.str();
      primclex.err_log() << "headers: " << col << std::endl;
      primclex.err_log() << "  Some queries may not be available for sampling at this time." << std::endl;
      throw std::runtime_error(ss.str());
    }

    for(int i = 0; i < col.size(); ++i) {

      std::string print_name = col[i];
      boost::algorithm::trim(print_name);

      std::tie(must_converge, prec) = _get_precision(it);

      // if 'must converge'
      if(must_converge) {
        ptr = new QueryMonteSampler(formatter, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new QueryMonteSampler(formatter, i, print_name, confidence(), data_maxlength);
      }

      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));

    }

    return result;

  }


}

#endif
