#ifndef CONFIGIOHULL_HH
#define CONFIGIOHULL_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/hull/Hull.hh"
namespace CASM {

  class Configuration;

  namespace ConfigIO_impl {
    /*
     */
    class BaseHullConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      BaseHullConfigFormatter(const std::string &_name, const std::string &_desc, const std::string &_dependent_prop) :
        BaseDatumFormatter<Configuration>(_name, _desc), m_dependent_prop(_dependent_prop) {}

      /// Initialize the convex hull and determine on-hull configurations
      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;

      //void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      //void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      //jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    protected:
      const BP::Geo &_hull() const {
        return m_hull;
      }
      const DataFormatter<Configuration> &_format() const {
        return m_format;
      }
      const std::map<std::string, bool> &_on_hull() const {
        return m_on_hull;
      }
      const std::vector<std::string> &_independent_props() const {
        return m_independent_props;
      }
      const Eigen::MatrixXd &_projection() const {
        return m_projection;
      }
      //const std::string &_dependent_prop const{ return m_dependent_prop;}
      //void _parse_args(const std::string &args, const std::string &_dep_prop);
    private:
      // specifies the dependent property to use (e.g., formation_energy or clex(formation_energy) )
      const std::string m_dependent_prop;
      mutable BP::Geo m_hull;

      // Matrix that describes subspace spanned by data
      mutable Eigen::MatrixXd m_projection;
      // names of on-hull configurations, all the bools are true--just used for fast access
      mutable std::map<std::string, bool> m_on_hull;
      // used to extract data from configurations for determining hull distance
      mutable DataFormatter<Configuration> m_format;

      // Parsed arguments
      //  -- what selection to use for constructing
      std::string m_selection;
      //  -- what variables to use for "axes" of convex hull. can be list of compositions or a single composition descriptor
      //     e.g., {"parametric_composition"}
      std::vector<std::string> m_independent_props;

    };

    /*
     *  BaseHullConfigFormatter handles initialization of data members, OnHullConfigFormatter uses that info to
     *  report whether it is on the hull.
     */

    class OnHullConfigFormatter: public BaseHullConfigFormatter {
    public:
      OnHullConfigFormatter(const std::string &_name, const std::string &_desc, const std::string _dependent_prop) :
        BaseHullConfigFormatter(_name, _desc, _dependent_prop) {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new OnHullConfigFormatter(*this);
      }

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

    protected:
      //inherits:
      // const BP::Geo& _hull() const{ return m_hull;}
      // const std::map<std::string, bool> &_on_hull() const{ return m_on_hull;}
      // const std::vector<std::string> &_independent_props const{ return m_independent_props;}
    };
    /*
     */

    class HullDistConfigFormatter: public BaseHullConfigFormatter {
    public:
      HullDistConfigFormatter(const std::string &_name, const std::string &_desc, const std::string _dependent_prop) :
        BaseHullConfigFormatter(_name, _desc, _dependent_prop) {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new HullDistConfigFormatter(*this);
      }

      bool validate(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

    protected:
      //inherits:
      // const BP::Geo& _hull() const{ return m_hull;}
      // const std::map<std::string, bool> &_on_hull() const{ return m_on_hull;}
      // const std::vector<std::string> &_independent_props const{ return m_independent_props;}
    };
  }
}
#endif

