#ifndef EIGENDATASTREAM_HH
#define EIGENDATASTREAM_HH

#include "casm/casm_io/dataformatter/DataStream.hh"
#include "casm/external/Eigen/Dense"
#include "casm/global/eigen.hh"

namespace CASM {

/// \ingroup DataFormatter
///
class MatrixXdDataStream : public DataStream {
 public:
  MatrixXdDataStream(DataStreamTraits _traits = none)
      : DataStream(_traits), m_matrix(1, 0), m_row(0), m_col(0) {}
  DataStream &operator<<(double _d) {
    if (m_col == m_matrix.cols()) {
      if (m_row > 0 || m_matrix.rows() != 1)
        throw std::runtime_error(
            "Attempting to stream non-rectangular data to Eigen::MatrixXd "
            "using MatrixXdDataStream, at row=" +
            std::to_string(m_row) + ", col=" + std::to_string(m_row) + "\n");
      m_matrix.conservativeResize(Eigen::NoChange, m_col + 1);
    }
    if (m_row == m_matrix.rows()) {
      m_matrix.conservativeResize(m_row + 1, Eigen::NoChange);
    }
    m_matrix(m_row, m_col++) = _d;
    return *this;
  }

  DataStream &operator<<(long _l) { return operator<<((double)_l); }

  virtual DataStream &newline() {
    if (fail() && _skipfail()) {
      clear_fail();
    } else
      ++m_row;
    m_col = 0;
    return *this;
  }

  const Eigen::MatrixXd &matrix() { return m_matrix; }

 protected:
  EigenIndex row() const { return m_row; }

 private:
  Eigen::MatrixXd m_matrix;
  EigenIndex m_row, m_col;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// \ingroup DataFormatter
///
class LabeledMatrixXdDataStream : public MatrixXdDataStream {
 public:
  LabeledMatrixXdDataStream(DataStreamTraits _traits = none)
      : MatrixXdDataStream(_traits) {}
  const std::vector<std::string> &labels() const { return m_labels; }

  DataStream &operator<<(const std::string &_str) {
    if (labels().size() == row())
      m_labels.push_back(_str);
    else if (row() + 1 == labels().size()) {
      std::cerr << "WARNING: Attempting to collect labeled data from a "
                   "datastream, but too many labels exist for row "
                << row() << "\n";
      DataStream::operator<<(failbit);
    }
    return *this;
  }

  DataStream &newline() {
    if (!fail() && labels().size() < row() + 1) {
      std::cerr << "WARNING: Attempting to collect labeled data from a "
                   "datastream, but no label was available for row "
                << row() << "\n";
      DataStream::operator<<(failbit);
    } else if (_skipfail() && fail() && row() + 1 == labels().size()) {
      m_labels.pop_back();
    }

    return MatrixXdDataStream::newline();
  }

 private:
  std::vector<std::string> m_labels;
};
}  // namespace CASM

#endif
