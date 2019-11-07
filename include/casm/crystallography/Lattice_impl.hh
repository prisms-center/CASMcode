#include "casm/crystallography/Lattice.hh"
#include "casm/CASM_global_enum.hh"
#include "casm/container/Counter.hh"

namespace CASM {
  namespace xtal {

    /**
     * This function generates a grid of points between max_radius and
     * min_radius. Additionally, it also fills up the points with a basis
     */

    template<typename CoordType, typename CoordType2>
    std::vector<CoordType> Lattice::gridstruc_build(double max_radius, double min_radius, std::vector<CoordType> basis, CoordType2 lat_point) {
      Eigen::Vector3i dim;
      dim = enclose_sphere(max_radius);
      EigenCounter<Eigen::Vector3i > grid_count(-dim, dim, Eigen::Vector3i(1));
      double min_dist, dist;
      std::vector<CoordType> gridstruc;
      Eigen::Vector3i temp;

      do {
        lat_point(FRAC) = grid_count();

        for(Index i = 0; i < basis.size(); i++) {
          CoordType tatom(basis[i] + lat_point);
          //get distance to closest basis site in the unit cell at the origin

          min_dist = 1e20;
          for(Index j = 0; j < basis.size(); j++) {
            dist = tatom.dist(basis[j]);
            if(dist < min_dist)
              min_dist = dist;
          }
          if(min_dist < min_radius) {
            continue;
          }
          if(min_dist < max_radius) {
            gridstruc.push_back(tatom);
          }
        }
      }
      while(++grid_count);

      return gridstruc;
    }

  }
}

