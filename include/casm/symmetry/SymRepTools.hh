#ifndef CASM_SymRepTools
#define CASM_SymRepTools

#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {
  namespace SymRepTools {
    struct IrrepWedge {
      IrrepWedge() {}
      IrrepWedge(Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                 std::vector<Index> _mult) :
        axes(_axes),
        mult(_mult) {}
      Eigen::MatrixXd axes;
      std::vector<Index> mult;
    };

    class SubWedge {
    public:
      SubWedge(std::vector<IrrepWedge> const &_iwedges) :
        m_iwedges(_iwedges),
        m_trans_mat(_subwedge_to_trans_mat(m_iwedges)) {

      }

      std::vector<IrrepWedge> const &irrep_wedges() const {
        return m_iwedges;
      }

      Eigen::MatrixXd const &trans_mat() const {
        return m_trans_mat;
      }

    private:
      std::vector<IrrepWedge> m_iwedges;
      Eigen::MatrixXd m_trans_mat;

      static Eigen::MatrixXd _subwedge_to_trans_mat(std::vector<IrrepWedge> const &_iwedges);
    };


    std::vector<IrrepWedge> irreducible_wedges(const SymGroup &head_group, SymGroupRepID id);
    std::vector<SubWedge> symrep_subwedges(SymGroup const &_group, SymGroupRepID id);
  }
}
#endif
