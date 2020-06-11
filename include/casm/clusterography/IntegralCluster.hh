#ifndef CASM_IntegralCluster
#define CASM_IntegralCluster

#include <vector>

#include "casm/clusterography/GenericCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

  class Structure;

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /** \defgroup IntegralCluster

      \brief Functions and classes related to IntegralCluster
      \ingroup Clusterography
      \ingroup CoordCluster
  */

  /* -- IntegralCluster ------------------------------------- */

  class IntegralCluster;

  /// traits, required for GenericCluster
  template <>
  struct traits<IntegralCluster> {
    typedef xtal::UnitCellCoord Element;
    //typedef Index size_type;
    //static const std::string name;
  };

  class IntegralCluster : public xtal::Translatable<GenericCluster<CRTPBase<IntegralCluster>>> {

  public:

    typedef Structure PrimType;
    typedef xtal::Translatable<GenericCluster<CRTPBase<IntegralCluster>>> Base;
    using Base::Element;
    using Base::size_type;

    IntegralCluster(PrimType const &prim);

    template<typename Iterator>
    IntegralCluster(PrimType const &prim, Iterator begin, Iterator end);

    const PrimType &prim() const;

    /// \brief Access vector of elements
    std::vector<Element> &elements();

    /// \brief const Access vector of elements
    const std::vector<Element> &elements() const;

    /// \brief Return the coordinate corresponding to element(i)
    xtal::Coordinate coordinate(size_type i) const;

    /// \brief Translate the cluster by a UnitCell translation
    IntegralCluster &operator+=(xtal::UnitCell trans);

  protected:

    friend GenericCluster<CRTPBase<IntegralCluster>>;

  private:

    std::vector<xtal::UnitCellCoord> m_element;
    const PrimType *m_prim_ptr;

  };
}

// TODO: remove
// #include <iosfwd>
//
// namespace CASM {
//
//   class jsonParser;
//   template<typename T> struct jsonConstructor;
//
//   /// \brief Print IntegralCluster to stream, using default Printer<IntegralCluster>
//   std::ostream &operator<<(std::ostream &sout, IntegralCluster const &clust);
//
//   /// \brief Write IntegralCluster to JSON object
//   jsonParser &to_json(IntegralCluster const &clust, jsonParser &json);
//
//   /// \brief Read from JSON
//   void from_json(IntegralCluster &clust, jsonParser const &json, double xtal_tol);
//
//   template<>
//   struct jsonConstructor<IntegralCluster> {
//
//     /// \brief Construct from JSON
//     static IntegralCluster from_json(jsonParser const &json, Structure const &prim, double xtal_tol);
//   };
//
// }

#endif
