#include <functional>
#include "casm/casm_io/EigenDataStream.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
//#include "casm/app/DirectoryStructure.hh"
//#include "casm/app/ProjectSettings.hh"

namespace CASM {

  namespace ConfigIO_impl {
    void BaseHullConfigFormatter::init(const Configuration &_tmplt) const {
      // initialize hull_props to look something like:
      //      {"comp", "formation_energy", "configname"}
      // or
      //      {"atom_frac(Au)", "atom_frac(Pt)", "clex(formation_energy)", "configname"}
      // in general:
      //      {$extensive_quantity1, $extensive_quantity2, ..., $energy_metric, "configname"}

      std::vector<std::string>  hull_props(m_independent_props);
      if(hull_props.size() == 0)
        hull_props.push_back("comp");
      hull_props.push_back(m_dependent_prop);
      // grabs the name for each row
      hull_props.push_back("configname");

      // setup DataStream
      //    -- fills a matrix so that each row corresponds to a configuration and pushes the confignames into an array
      //    -- columns of matrix correspond to {$extensive_quantity1, $extensive_quantity2, ..., $energy_metric}
      //    -- initialized with DataStram::skipfail -> skip configs that are missing any of the necessary data
      LabeledMatrixXdDataStream mat_wrapper(DataStream::skipfail);
      m_format = ConfigIOParser::parse(hull_props);

      // Cases control which configurations to use for obtaining hull data
      if(m_selection == "all")
        mat_wrapper << m_format(_tmplt.get_primclex().config_cbegin(), _tmplt.get_primclex().config_cend());
      else if(m_selection == "MASTER")
        mat_wrapper << m_format(_tmplt.get_primclex().selected_config_cbegin(), _tmplt.get_primclex().selected_config_cend());
      else {
        ConstConfigSelection select(_tmplt.get_primclex(), m_selection);
        mat_wrapper << m_format(select.selected_config_cbegin(), select.selected_config_cend());
      }

      EigenIndex numcomps(mat_wrapper.matrix().cols() - 1), rank(0);
      // Make the data matrix full-rank (hopefully)
      Eigen::MatrixXd compcovar(mat_wrapper.matrix().leftCols(numcomps).transpose()*mat_wrapper.matrix().leftCols(numcomps));
      Eigen::JacobiSVD<Eigen::MatrixXd> tsvd(compcovar, Eigen::ComputeFullU);
      // get rank manually (eigen's rank() method uses relative tolerances. absolute tolerances are better for our purposes...)
      for(Index i = 0; i < tsvd.singularValues().size(); i++) {
        if(!almost_zero(tsvd.singularValues()[i])) {
          rank++;
        }
        else//singular values are always positive and decreasing
          break;
      }
      //std::cout << "compcovar matrix is\n" << compcovar << "\n\n and U matrix is \n" << tsvd.matrixU() << "\n\n";
      m_projection.resize(rank + 1, numcomps + 1);
      m_projection << tsvd.matrixU().topRows(rank), Eigen::MatrixXd::Zero(rank, 1),
                   Eigen::MatrixXd::Zero(1, numcomps), 1.0;

      //std::cout << "Size before: " << mat_wrapper.matrix().rows() << ", " << mat_wrapper.matrix().cols() << "\n";

      // ***The following two lines turn off optimal projection***
      // ***Delete them after convex hull is fixed***
      m_projection = Eigen::MatrixXd::Identity(numcomps + 1, numcomps + 1);
      rank = numcomps;

      Eigen::MatrixXd reduced_mat(m_projection * mat_wrapper.matrix().transpose());
      //std::cout << "Size after: " << reduced_mat.rows() << ", " << reduced_mat.cols() << "\n";
      //std::cout << "reduced_mat is \n" << reduced_mat.transpose() << "\n\n and projection is \n" << m_projection << "\n\n";
      m_hull.set_verbosity(0);

      //(using check_for_repeats=false doesn't seem to be an improvement)
      m_hull.reset_points(reduced_mat, true);

      if(!m_hull.calc_CH()) { //calculates hull
        throw std::runtime_error("Failure to construct convex hull from selection " + m_selection
                                 + " for formatted output!\n");
      }

      // define which direction is down
      Eigen::VectorXd down = Eigen::VectorXd::Zero(rank + 1);
      down(rank) = -1;
      m_hull.CH_bottom(down);


      /*
      for(Index i=0; i<m_hull.size(); i++){
        std::cout << "  i" << i<<":  " << reduced_mat.col(i).transpose() << "   " << mat_wrapper.matrix().row(i) << "   " <<m_hull.pos(i).transpose() << "   " << m_hull.CH_dist_to_hull(m_hull.pos(i)) << "\n";
        }*/


      // record names of on-hull configs
      auto hull_inds = m_hull.CH_verts_indices();
      for(Index i = 0; i < hull_inds.size(); i++) {
        m_on_hull[mat_wrapper.labels()[hull_inds[i]]] = true;
      }

    }

    //****************************************************************************************
    bool BaseHullConfigFormatter::parse_args(const std::string &args) {
      if(m_independent_props.size() || m_selection.size())
        return false;
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
      if(splt_vec.size() > 2) {
        throw std::runtime_error("Attempted to initialize format tag " + name()
                                 + " with " + std::to_string(splt_vec.size()) + " arguments ("
                                 + args + "), but no more than 2 arguments allowed.\n");
        return false;
      }
      m_independent_props.push_back(splt_vec.size() < 2 ? "comp" : splt_vec[1]);
      m_selection = splt_vec.size() < 1 ? "all" : splt_vec[0];
      return true;
    }

    //****************************************************************************************

    std::string BaseHullConfigFormatter::short_header(const Configuration &_tmplt) const {
      std::stringstream t_ss;
      t_ss << name() << "(" << m_selection << ",";
      for(Index i = 0; i < m_independent_props.size(); i++)
        t_ss << m_independent_props[i];
      t_ss << ")";
      return t_ss.str();
    }

    //****************************************************************************************

    void OnHullConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      _stream << (_on_hull().find(_config.name()) != _on_hull().cend());
    }

    //****************************************************************************************

    void OnHullConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      _stream << (_on_hull().find(_config.name()) != _on_hull().cend());
    }

    //****************************************************************************************

    jsonParser &OnHullConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      json = _on_hull().find(_config.name()) != _on_hull().cend();
      return json;
    }

    //****************************************************************************************

    bool HullDistConfigFormatter::validate(const Configuration &_config)const {
      return _format().validate(_config);
    }

    //****************************************************************************************

    void HullDistConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      MatrixXdDataStream mat_wrapper;
      mat_wrapper << _format()(_config);
      if(mat_wrapper.fail())
        _stream << DataStream::failbit << double(NAN);
      else
        _stream << _hull().CH_dist_to_hull(Eigen::VectorXd(_projection()*mat_wrapper.matrix().transpose()));
    }

    //****************************************************************************************

    void HullDistConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      MatrixXdDataStream mat_wrapper;
      mat_wrapper << _format()(_config);
      //std::cout << "mat_wrapper.matrix() is " << mat_wrapper.matrix() << "\n";
      //std::cout << "rotated, is " << (_projection()*mat_wrapper.matrix().transpose()).transpose() << "\n";
      //std::cout << "projection is \n" << _projection() << "\n";
      if(mat_wrapper.fail())
        _stream << "unknown";
      else
        _stream << _hull().CH_dist_to_hull(Eigen::VectorXd(_projection()*mat_wrapper.matrix().transpose()));
    }

    //****************************************************************************************

    jsonParser &HullDistConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      MatrixXdDataStream mat_wrapper;
      mat_wrapper << _format()(_config);
      if(mat_wrapper.fail())
        json = "unknown";
      else
        json = _hull().CH_dist_to_hull(Eigen::VectorXd(_projection() * mat_wrapper.matrix().transpose()));
      return json;
    }
  }
}

