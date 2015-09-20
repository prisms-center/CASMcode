#ifndef CASM_Hull
#define CASM_Hull

#include "casm/BP_C++/BP_Geo.hh"
#include "casm/clex/ConfigSelection.hh"

namespace CASM {

  ///Run through configurations, and make a hull out of them, keeping track of the valid configurations.
  template<typename ConSelectOutputIterator>
  bool populate_convex_hull(BP::Geo &hull,
                            ConfigSelectionIterator<false, false> begin,
                            ConfigSelectionIterator<false, false> end,
                            ConSelectOutputIterator valid_config);

  ///Run through all the configurations and update hull variables into generated Properties. Also return hull in json.
  jsonParser update_hull_props(PrimClex &primclex,
                               ConfigSelectionIterator<false, false> begin,
                               ConfigSelectionIterator<false, false> end,
                               double geo_tol = 0.0005);



  /**
   * For every configuration that's selected, check to see if it has relaxed_energy. If it does, then include
   * its energy and composition into a list. Once we have all the energies, we use the list to generate the convex
   * hull. Because we may not use all the configurations, also keep track of the indices of supercells
   * and configurations that were used to make hull.
   */
  template<typename ConSelectOutputIterator>
  bool populate_convex_hull(BP::Geo &hull,
                            ConfigSelectionIterator<false, false> begin,
                            ConfigSelectionIterator<false, false> end,
                            ConSelectOutputIterator valid_config) {

    //As we go through, we store the energies and compositions. These will then get put into an Eigen matrix
    Array<double> energies;
    Array<Eigen::VectorXd> compositions;

    for(auto it = begin; it != end; ++it) {
      if(it.selected()) {
        if(it->delta_properties().contains("relaxed_energy")) {
          //If you made it in here then the configuration is turned on and has a relaxed energy, so we store values
          double good_config_energy = it->delta_properties()["relaxed_energy"].get<double>();
          Eigen::VectorXd good_config_composition = it->get_param_composition();

          energies.push_back(good_config_energy);
          compositions.push_back(good_config_composition);

          //Remember this configuration
          *valid_config++ = it;
        }
        else {
          std::cerr << "WARNING in PrimClex::populate_convex_hull" << std::endl
                    << "Requested to include " << it->get_path().string()
                    << " but there is no relaxed energy for this configuration." << std::endl
                    << "Skipping..." << std::endl;
        }
      }
    }

    if(energies.size() == 0 || compositions.size() == 0) {
      std::cerr << "ERROR in PrimClex::update_hull" << std::endl;
      std::cerr << "Best case scenario: none of your configurations have a relaxed energy. Be sure to update the properties" << std::endl;
      exit(323);
    }

    //We have all the energies and compositions now, which we load into a matrix (unsigned long seems overkill but w/e)
    Index rows = compositions[0].size() + 1; //we have number of compositions + energy rows
    Index cols = compositions.size();     //we have how ever many configurations columns
    Eigen::MatrixXd enercomps(rows, cols);

    //The matrix is sized correctly, but is empty. Now we fill in values
    //Loop over every configuration (columns)
    for(Index i = 0; i < cols; i++) {
      //First stick the energy (0th row value)
      enercomps(0, i) = energies[i];
      //Loop over the compositions (rows after energy, start at 1th row)
      for(Index j = 1; j < rows; j++) {
        //Fill in composition values (all the rows after the first energy row)
        enercomps(j, i) = compositions[i](j - 1);
      }
    }

    //Your matrix looks something like:
    //
    //    conf0 conf1 conf2 conf3  ...
    //E    ...   ...   ...   ...   ...
    //x0   ...   ...   ...   ...   ...
    //x1   ...   ...   ...   ...   ...
    //x2   ...   ...   ...   ...   ...
    //...  ...   ...   ...   ...   ...

    //Let the BP magic begin
    //We start by making a BP::Geo object for our hull
    hull.reset_points(enercomps, true); //Check for repeats is good?
    bool hull_found = hull.calc_CH();

    return hull_found;
  }

}

#endif
