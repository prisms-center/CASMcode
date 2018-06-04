#include "casm/kinetics/DiffusionTransformationTraits.hh"
#include "casm/kinetics/PrimPeriodicDiffTransOrbitIO.hh"
#include "casm/database/Selected.hh"

#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/casm_io/DataFormatter_impl.hh"
#include "casm/casm_io/DataFormatterTools_impl.hh"

namespace CASM {

  template class BaseDatumFormatter<PrimPeriodicDiffTransOrbit>;
  template class DataFormatterOperator<bool, std::string, PrimPeriodicDiffTransOrbit>;
  template class DataFormatterOperator<bool, bool, PrimPeriodicDiffTransOrbit>;
  template class DataFormatterOperator<bool, double, PrimPeriodicDiffTransOrbit>;
  template class DataFormatterOperator<double, double, PrimPeriodicDiffTransOrbit>;
  template class DataFormatterOperator<Index, double, PrimPeriodicDiffTransOrbit>;
  template class DataFormatter<PrimPeriodicDiffTransOrbit>;
  //template std::string DataFormatterOperator<double, double, PrimPeriodicDiffTransOrbit >::short_header<std::string>(PrimPeriodicDiffTransOrbit const&) const;
  //template class DataFormatterOperator<double, double, PrimPeriodicDiffTransOrbit >::parse_args(std::string const&);
  template class DataFormatterDictionary<PrimPeriodicDiffTransOrbit>;


  namespace PrimPeriodicDiffTransOrbitIO {

    // Contains implementations

    const std::string Contains::Name = "contains";

    const std::string Contains::Desc =
      "Checks to see if all of the listed species are in the prototype"
      "of the orbit"
      "Requires an argument which is a string of the form:"
      "contains(species1,species2,...)";

    /// \brief Returns true if all of the species in m_search_list are in the prototype
    ///  of the orbit
    bool Contains::evaluate(const PrimPeriodicDiffTransOrbit &orbit) const {
      auto speciesmap = orbit.prototype().species_count();
      bool ret_val = true;
      for(auto it = m_search_list.begin(); it != m_search_list.end(); ++it) {
        bool tmp = false;
        for(auto it2 = speciesmap.begin(); it2 != speciesmap.end(); ++it2) {
          if(it2->first.name() == *it && it2->second != 0) {
            tmp = true;
          }
        }
        ret_val = (ret_val && tmp);
      }
      return ret_val;
    };

    /// \brief Expects 'contains("Specie1","Specie2",...)'
    bool Contains::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
      m_search_list = splt_vec;
      return true;
    };

    //Generic Formatter functions

    GenericDiffTransOrbitFormatter<Index> multiplicity() {
      return GenericDiffTransOrbitFormatter<Index>("multiplicity",
                                                   "Symmetric multiplicity of the diffusion transformation, excluding translational equivalents.",
      [](const PrimPeriodicDiffTransOrbit & orbit)->Index {
        return orbit.size();
      });
    }

    GenericDiffTransOrbitFormatter<Index> cluster_size() {
      return GenericDiffTransOrbitFormatter<Index>("cluster_size",
                                                   "Size of the cluster corresponding to the diffusion transformation.",
      [](const PrimPeriodicDiffTransOrbit & orbit)->Index {
        return orbit.prototype().size();
      });
    }

    GenericDiffTransOrbitFormatter<double> min_length() {
      return GenericDiffTransOrbitFormatter<double>("min_length",
                                                    "Returns the minimum length between two sites in the cluser of the diffusion transformation",
      [](const PrimPeriodicDiffTransOrbit & orbit)-> double {
        return orbit.prototype().min_length();
      });
    }

    GenericDiffTransOrbitFormatter<double> max_length() {
      return GenericDiffTransOrbitFormatter<double>("max_length",
                                                    "Returns the maximum length between two sites in the cluser of the diffusion transformation",
      [](const PrimPeriodicDiffTransOrbit & orbit)-> double {
        return orbit.prototype().max_length();
      });
    }

    GenericDiffTransOrbitFormatter<double> min_dist_to_path() {
      return GenericDiffTransOrbitFormatter<double>("min_dist_to_path",
                                                    "Returns the distance of the closest atom to the diffusion transformation",
      [](const PrimPeriodicDiffTransOrbit & orbit)-> double {
        return min_dist_to_path(orbit.prototype());
      });
    }

    GenericDiffTransOrbitFormatter<bool> path_collision() {
      return GenericDiffTransOrbitFormatter<bool>("path_collision",
                                                  "Returns true(1) if there is a chance that two moving species would collide on the designated path",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> bool {
        return path_collision(orbit.prototype());
      });
    }

    GenericDiffTransOrbitFormatter<std::string> species_list() {
      return GenericDiffTransOrbitFormatter<std::string>("species_list",
                                                         "Returns the species that move in this diffusion transformation",
      [](const PrimPeriodicDiffTransOrbit & orbit)-> std::string {

        auto speciesmap = orbit.prototype().species_count();
        std::string result = "";
        for(auto it = speciesmap.begin(); it != speciesmap.end(); ++it) {
          for(int i = 0; i < it->second; i++) {
            result += (it->first.name() + " ");
          }
        }
        return result;
      });
    }

    GenericDiffTransOrbitFormatter<std::string> diff_trans_orbitname() {
      return GenericDiffTransOrbitFormatter<std::string>("diff_trans_name",
                                                         "Returns the name of the diffusion transformation orbit",
      [](const PrimPeriodicDiffTransOrbit & orbit)-> std::string {
        return orbit.name();
      });
    }

    GenericDiffTransOrbitFormatter<double> min_activation_barrier() {
      return GenericDiffTransOrbitFormatter<double>("min_activation_barrier",
                                                    "Returns the locally cluster expanded activation barrier associated with this diffusion transformation",
      [](const PrimPeriodicDiffTransOrbit & orbit)-> double {
        //TEMPORARY FUNCTION UNTIL ACTIVATION BARRIER CAN BE CALCULATED
        return min_dist_to_path(orbit.prototype());
      }, [](const PrimPeriodicDiffTransOrbit & orbit)-> bool {
        return false;
      });
    }
  }



  template<>
  StringAttributeDictionary<PrimPeriodicDiffTransOrbit> make_string_dictionary<PrimPeriodicDiffTransOrbit>() {

    using namespace PrimPeriodicDiffTransOrbitIO;
    StringAttributeDictionary<PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      name<PrimPeriodicDiffTransOrbit>(),
      species_list(),
      diff_trans_orbitname()
    );

    return dict;
  }

  template<>
  BooleanAttributeDictionary<PrimPeriodicDiffTransOrbit> make_boolean_dictionary<PrimPeriodicDiffTransOrbit>() {

    using namespace PrimPeriodicDiffTransOrbitIO;
    BooleanAttributeDictionary<PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      Contains(),
      DB::Selected<PrimPeriodicDiffTransOrbit>(),
      path_collision()
    );

    return dict;
  }

  template<>
  IntegerAttributeDictionary<PrimPeriodicDiffTransOrbit> make_integer_dictionary<PrimPeriodicDiffTransOrbit>() {

    using namespace PrimPeriodicDiffTransOrbitIO;
    IntegerAttributeDictionary<PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      multiplicity(),
      cluster_size()
    );

    return dict;
  }

  template<>
  ScalarAttributeDictionary<PrimPeriodicDiffTransOrbit> make_scalar_dictionary<PrimPeriodicDiffTransOrbit>() {

    using namespace PrimPeriodicDiffTransOrbitIO;
    ScalarAttributeDictionary<PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      min_length(),
      max_length(),
      min_dist_to_path()
    );

    return dict;
  }

  template<>
  VectorXdAttributeDictionary<PrimPeriodicDiffTransOrbit> make_vectorxd_dictionary<PrimPeriodicDiffTransOrbit>() {
    return VectorXdAttributeDictionary<PrimPeriodicDiffTransOrbit>();
  }


  template<>
  VectorXiAttributeDictionary<PrimPeriodicDiffTransOrbit> make_vectorxi_dictionary<PrimPeriodicDiffTransOrbit>() {
    return VectorXiAttributeDictionary<PrimPeriodicDiffTransOrbit>();
  }

}
