#ifndef CASM_GrandCanonicalSettings
#define CASM_GrandCanonicalSettings

#include "casm/monte_carlo/MonteSettings.hh"

namespace CASM {

  class GrandCanonicalConditions;

  class GrandCanonicalSettings : public EquilibriumMonteSettings {

  public:

    /// \brief Default constructor
    GrandCanonicalSettings() {}

    /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
    GrandCanonicalSettings(const fs::path &read_path);


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

    /// \brief Given a settings jsonParser figure out what the project clex settings to use are:
    std::string clex() const;

    /// \brief Given a settings jsonParser figure out what the project bset settings to use are:
    std::string bset() const;

    /// \brief Given a settings jsonParser figure out what the project calctype settings to use are:
    std::string calctype() const;

    /// \brief Given a settings jsonParser figure out what the project ref settings to use are:
    std::string ref() const;

    /// \brief Given a settings jsonParser figure out what the project eci settings to use are:
    std::string eci() const;



    // --- Sampler settings ---------------------

    /// \brief Construct MonteSamplers as specified in the MonteSettings
    template<typename SamplerInsertIterator>
    SamplerInsertIterator samplers(const PrimClex &primclex, SamplerInsertIterator result) const;

    /// \brief Return true if all correlations should be sampled
    bool all_correlations() const;


  private:

    CompositionConverter m_comp_converter;

    GrandCanonicalConditions _conditions(std::string name) const;
    GrandCanonicalConditions _conditions(const jsonParser &json) const;

    template<typename jsonParserIteratorType>
    std::tuple<bool, double> _get_precision(jsonParserIteratorType it, std::string input_name) const;

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

  };


  /// \brief Construct MonteSamplers as specified in the MonteSettings
  ///
  /// The requested MonteSamplers are inserted in 'result' as
  ///   std::pair<std::string, notstd::cloneable_ptr<MonteSampler> >
  ///
  template<typename SamplerInsertIterator>
  SamplerInsertIterator GrandCanonicalSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result) const {

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



    // construct scalar property samplers
    {
      std::vector<std::string> possible = {
        "formation_energy",
        "potential_energy"
      };


      try {

        for(auto it = t_measurements.cbegin(); it != t_measurements.cend(); it++) {
          prop_name = (*it)["quantity"].get<std::string>();

          // check if property found is in list of possible scalar properties
          if(std::find(possible.cbegin(), possible.cend(), prop_name) != possible.cend()) {

            std::tie(must_converge, prec) = _get_precision(it, prop_name);

            // if 'must converge'
            if(must_converge) {
              ptr = new ScalarMonteSampler(prop_name, prop_name, prec, confidence(), data_maxlength);
            }
            else {
              ptr = new ScalarMonteSampler(prop_name, prop_name, confidence(), data_maxlength);
            }

            *result++ = std::make_pair(prop_name, notstd::cloneable_ptr<MonteSampler>(ptr));

          }

        }

      }
      catch(std::runtime_error &e) {
        std::cerr << "ERROR in 'MonteSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result)'\n" << std::endl;
        std::cerr << "Error reading [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        throw;
      }
    }


    // construct vector properties to sample
    {
      std::vector<std::string> possible = {
        "comp",
        "comp_n",
        "site_frac",
        "atom_frac",
        "all_correlations",
        "non_zero_eci_correlations"
      };

      try {

        for(auto it = t_measurements.cbegin(); it != t_measurements.cend(); it++) {

          std::string input_name = (*it)["quantity"].get<std::string>();

          // check if property found is in list of possible vector properties
          if(std::find(possible.cbegin(), possible.cend(), input_name) != possible.cend()) {

            // construct MonteSamplers for 'comp'
            if(input_name == "comp") {

              result = _make_comp_samplers(primclex, it, result);

            }

            // construct MonteSamplers for 'comp_n'
            else if(input_name == "comp_n") {

              result = _make_comp_n_samplers(primclex, it, result);

            }

            // construct MonteSamplers for 'site_frac'
            else if(input_name == "site_frac") {

              result = _make_site_frac_samplers(primclex, it, result);

            }

            // construct MonteSamplers for 'atom_frac'
            else if(input_name == "atom_frac") {

              result = _make_atom_frac_samplers(primclex, it, result);

            }

            // construct MonteSamplers for 'all_correlations'
            else if(input_name == "all_correlations") {

              result = _make_all_correlations_samplers(primclex, it, result);

            }

            // construct MonteSamplers for 'non_zero_eci_correlations'
            else if(input_name == "non_zero_eci_correlations") {

              result = _make_non_zero_eci_correlations_samplers(primclex, it, result);

            }
          }
        }
      }
      catch(std::runtime_error &e) {
        std::cerr << "ERROR in 'MonteSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result)'\n" << std::endl;
        std::cerr << "Error reading [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        throw;
      }
    }

    return result;
  }

  template<typename jsonParserIteratorType>
  std::tuple<bool, double> GrandCanonicalSettings::_get_precision(jsonParserIteratorType it, std::string input_name) const {
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

      std::tie(must_converge, prec) = _get_precision(it, "comp");

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

      std::tie(must_converge, prec) = _get_precision(it, "comp_n");

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

      std::tie(must_converge, prec) = _get_precision(it, "site_frac");

      // if 'must converge'
      if(must_converge) {
        ptr = new SiteFracMonteSampler(i, primclex.prim().basis.size(), print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new SiteFracMonteSampler(i, primclex.prim().basis.size(), print_name, confidence(), data_maxlength);
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

      // sample for vacancy components
      if(is_vacancy(primclex.composition_axes().components()[i])) {
        vacancy_index = i;
        break;
      }
    }


    for(size_type i = 0; i < primclex.composition_axes().components().size(); i++) {

      // sample for non-vacancy components
      if(!is_vacancy(primclex.composition_axes().components()[i])) {

        print_name = std::string("atom_frac(") + primclex.composition_axes().components()[i] + ")";

        std::tie(must_converge, prec) = _get_precision(it, "atom_frac");

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

    for(size_type i = 0; i < primclex.global_clexulator().corr_size(); i++) {

      prop_name = "corr";
      print_name = std::string("corr(") + std::to_string(i) + ")";

      std::tie(must_converge, prec) = _get_precision(it, "all_correlations");

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

    ECIContainer _eci = read_eci(primclex.dir().eci(clex(), calctype(), ref(), bset(), eci()));

    for(size_type ii = 0; ii < _eci.index().size(); ii++) {

      prop_name = "corr";

      // store non-zero eci index in i
      size_type i = _eci.index()[ii];

      print_name = std::string("corr(") + std::to_string(i) + ")";

      std::tie(must_converge, prec) = _get_precision(it, "non_zero_eci_correlations");

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


}

#endif

