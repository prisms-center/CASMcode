#include "casm/casm_io/EigenDataStream.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  namespace ConfigIO {
    void Novelty::init(const Configuration &_tmplt) const {

      m_format.clear();
      m_format.push_back(Corr());

      // setup DataStream
      //    -- fills a matrix so that each row corresponds to a configuration and column corresponds to correlation value
      MatrixXdDataStream mat_wrapper(DataStream::skipfail);

      // Cases control which configurations to use for obtaining novelty data
      if(m_selection == "all")
        mat_wrapper << m_format(_tmplt.get_primclex().config_cbegin(), _tmplt.get_primclex().config_cend());
      else if(m_selection == "MASTER" || m_selection.empty())
        mat_wrapper << m_format(_tmplt.get_primclex().selected_config_cbegin(), _tmplt.get_primclex().selected_config_cend());
      else {
        ConstConfigSelection select(_tmplt.get_primclex(), m_selection);
        mat_wrapper << m_format(select.selected_config_cbegin(), select.selected_config_cend());
      }

      m_avg_corr = mat_wrapper.matrix().colwise().sum().transpose() / double(mat_wrapper.matrix().cols());
      Eigen::MatrixXd tcovar = mat_wrapper.matrix().transpose() * mat_wrapper.matrix() / double(mat_wrapper.matrix().cols()) - m_avg_corr * m_avg_corr.transpose();

      m_gram_mat = (double(tcovar.rows()) * tcovar + 1E-6 * Eigen::MatrixXd::Identity(tcovar.rows(), tcovar.cols())).inverse();

    }

    //****************************************************************************************
    bool Novelty::parse_args(const std::string &args) {
      if(m_selection.size())
        return false;
      if(args.empty())
        m_selection = "MASTER";
      else
        m_selection = args;
      return true;
    }

    // ****************************************************************************************

    std::string Novelty::short_header(const Configuration &_tmplt) const {
      return name() + "(" + (m_selection.empty() ? "MASTER" : m_selection) + ")";
    }

    // ****************************************************************************************
    /*
        void Novelty::inject(const Configuration &_config, DataStream &_stream, Index) const {
          _stream << _evaluate(_config);
        }

        // ****************************************************************************************

        void Novelty::print(const Configuration &_config, std::ostream &_stream, Index) const {
          _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
          _stream.precision(8);

          _stream << _evaluate(_config);
        }

        // ****************************************************************************************

        jsonParser &Novelty::to_json(const Configuration &_config, jsonParser &json)const {
          json = _evaluate(_config);
          return json;
        }
    */

    //****************************************************************************************

    double Novelty::evaluate(const Configuration &_config) const {
      VectorDataStream<double> tstr;
      tstr << m_format(_config);
      Eigen::Map<const Eigen::VectorXd> corr(tstr.vector().data(), tstr.vector().size());
      return sqrt(double((corr - m_avg_corr).transpose() * m_gram_mat * (corr - m_avg_corr)));
    }
  }
}

