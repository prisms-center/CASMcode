#ifndef TEMPLATE_ALGORITHMS_HH
#define TEMPLATE_ALGORITHMS_HH

#include <cmath>

//#include "casm/../CASM_global_definitions.cc"

namespace CASM {

  /*
  template<class Container, class T>
  class GaussJordan {
  public:
    static void eliminate(Array<Container> &rows, double elim_tol = TOL) {
      if(rows.size() == 0) return;


      int i, j, nc, nr, min_row_num = 0;
      Container trow;

      //Loop over columns
      for(nc = 0; nc < rows[0].size(); nc++) {

        //Loop over rows
        for(nr = min_row_num; nr < rows.size(); nr++) {
          if(std::abs(rows[nr][nc]) < TOL) {
            rows[nr][nc] = 0;
            continue;
          }

          trow = rows[nr];
          rows[nr] = rows[min_row_num];
          rows[min_row_num] = trow;

          rows[min_row_num] /= rows[min_row_num][nc];
          break;
        }

        if(nr >= rows.size()) continue;

        for(i = 0; i < rows.size(); i++) {
          if(i != min_row_num && elim_tol < std::abs(rows[i][nc]))
            rows[i] -= rows[i][nc] * rows[min_row_num];

        }

        min_row_num++;
      }

      for(nr = rows.size() - 1; nr >= 0; nr--) {
        for(nc = 0; nc < rows[nr].size(); nc++) {
          if(elim_tol < std::abs(rows[nr][nc]))
            break;
          rows[nr][nc] = 0;
        }
        if(nc == rows[nr].size())
          rows.pop_back();
      }
      return;
    };

  };
  */
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
  template<class Container, class T>
  class GaussElim {
  public:
    static void eliminate(Array<Container> &rows, double elim_tol = TOL) {
      if(rows.size() == 0) return;


      int i, j, nc, nr, min_row_num = 0;
      Container trow;
      for(nc = 0; nc < rows[0].size(); nc++) {

        for(nr = min_row_num; nr < rows.size(); nr++) {
          if(std::abs(rows[nr][nc]) < TOL) {
            rows[nr][nc] = 0;
            continue;
          }

          trow = rows[nr];
          rows[nr] = rows[min_row_num];
          rows[min_row_num] = trow;

          rows[min_row_num] /= rows[min_row_num][nc];
          break;
        }

        if(nr >= rows.size()) continue;

        for(i = min_row_num + 1; i < rows.size(); i++) {
          if(elim_tol < std::abs(rows[i][nc]))
            rows[i] -= rows[i][nc] * rows[min_row_num];

        }

        min_row_num++;
      }

      for(nr = rows.size() - 1; nr >= 0; nr--) {
        for(nc = 0; nc < rows[nr].size(); nc++) {
          if(elim_tol < std::abs(rows[nr][nc]))
            break;
          rows[nr][nc] = 0;
        }
        if(nc == rows[nr].size())
          rows.pop_back();
      }
      return;
    };


  };
  */
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<class Container, class T>
  class GramSchmidt {
  public:
    static void orthogonalize(Array<Container> &rows) {

      Index i, j;
      T tcoeff;
      for(i = 0; i < rows.size(); i++) {
        for(j = 0; j < i; j++)
          rows[i] -= dot(rows[i], rows[j]) * rows[j];

        tcoeff = norm(rows[i]);
        if(TOL < tcoeff)
          rows[i] /= tcoeff;

        else {
          rows.remove(i);
          i--;
        }

      }
      return;
    }
  };


};
#endif /* TEMPLATE_ALGORITHMS_HPP */
