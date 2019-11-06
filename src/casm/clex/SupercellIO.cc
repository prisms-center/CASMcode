#include "casm/clex/SupercellIO.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/Supercell_impl.hh"
#include "casm/database/Selected_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"

namespace CASM {

  template class BaseDatumFormatter<Supercell>;
  template class DataFormatterOperator<bool, std::string, Supercell>;
  template class DataFormatterOperator<bool, bool, Supercell>;
  template class DataFormatterOperator<bool, double, Supercell>;
  template class DataFormatterOperator<double, double, Supercell>;
  template class DataFormatterOperator<Index, double, Supercell>;
  template class DataFormatter<Supercell>;
  template class DataFormatterDictionary<Supercell>;

  namespace ScelIO {

    // --- template<typename Base> class SupercellCheckBase ---

    template<typename Base>
    SupercellCheckBase<Base>::SupercellCheckBase(std::string name, std::string desc) :
      Base(name, desc),
      m_refcell(nullptr),
      m_last_result(notstd::make_cloneable<result_type>()),
      m_last_scel(nullptr),
      m_last_unit(nullptr) {}

    /// \brief Expects arguments of the form 'is_supercell_of(scelname)'
    template<typename Base>
    bool SupercellCheckBase<Base>::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

      if(splt_vec.size() != 1) {
        std::stringstream ss;
        ss << this->name() << " expected 1 argument.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }

      m_refcell_name = args;
      return true;
    }

    /// \brief Set pointer to ref supercell
    template<typename Base>
    bool SupercellCheckBase<Base>::init(const Supercell &_tmplt) const {
      m_refcell = &*_tmplt.primclex().db<Supercell>().find(m_refcell_name);
      return true;
    }

    /// \brief col_header returns: {'name(refcell_name)'}
    template<typename Base>
    std::vector<std::string>
    SupercellCheckBase<Base>::col_header(const Supercell &_tmplt) const {
      return std::vector<std::string> {this->name() + "(" + m_refcell_name + ")"};
    }

    /// Call is_supercell using prim.factor_group() to try possible orientations
    ///
    /// Returns (bool, SymOp op, Eigen::MatrixXi T) with:
    /// - scel is supercell of unit?
    /// - If true: scel.lattice() == apply(op, unit.lattice()) * T
    ///
    template<typename Base>
    const typename SupercellCheckBase<Base>::result_type &
    SupercellCheckBase<Base>::_evaluate(const Supercell &scel, const Supercell &unit) const {
      if(&scel != m_last_scel || &unit != m_last_unit) {

        auto res = is_supercell(
                     scel.lattice(),
                     unit.lattice(),
                     unit.prim().factor_group().begin(),
                     unit.prim().factor_group().end(),
                     unit.crystallography_tol());

        *m_last_result = std::make_tuple(
                           res.first != unit.prim().factor_group().end(),
                           *res.first,
                           iround(res.second));
      }
      return *m_last_result;
    }


    // --- IsSupercellOf ---

    const std::string IsSupercellOf::Name = "is_supercell_of";
    const std::string IsSupercellOf::Desc =
      "Returns true for all supercells that are a supercell the specified supercell. "
      "All re-orientations allowed by the crystal point group are checked. "
      "Ex: 'is_supercell_of(SCELV_A_B_C_D_E_F)'";

    IsSupercellOf::IsSupercellOf() :
      SupercellCheckBase(Name, Desc) {}

    bool IsSupercellOf::evaluate(const Supercell &scel) const {
      return std::get<0>(_evaluate(scel, *m_refcell));
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<IsSupercellOf> IsSupercellOf::clone() const {
      return std::unique_ptr<IsSupercellOf>(this->_clone());
    }

    /// \brief Clone using copy constructor
    IsSupercellOf *IsSupercellOf::_clone() const {
      return new IsSupercellOf(*this);
    }


    // --- IsUnitcellOf ---

    const std::string IsUnitcellOf::Name = "is_unitcell_of";
    const std::string IsUnitcellOf::Desc =
      "Returns true for all supercells that can tile the specified supercell. "
      "All re-orientations allowed by the crystal point group are checked. "
      "Ex: 'is_unitcell_of(SCELV_A_B_C_D_E_F)'";

    IsUnitcellOf::IsUnitcellOf() :
      SupercellCheckBase(Name, Desc) {}

    bool IsUnitcellOf::evaluate(const Supercell &unit) const {
      return std::get<0>(_evaluate(*m_refcell, unit));
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<IsUnitcellOf> IsUnitcellOf::clone() const {
      return std::unique_ptr<IsUnitcellOf>(this->_clone());
    }

    /// \brief Clone using copy constructor
    IsUnitcellOf *IsUnitcellOf::_clone() const {
      return new IsUnitcellOf(*this);
    }


    // --- TransfMat ---

    const std::string TransfMat::Name = "transf_mat";
    const std::string TransfMat::Desc =
      "For all supercells, S, returns the transformation matrix, T, that can be "
      "used to create S from the specified unit cell, U, if possible, i.e. "
      "S.lat = (op*U.lat)*T, where 'op' is an element of the crystal "
      "point group, and lattices are represented by column vector matrices. "
      "T is returned in column-major form: (T00, T10, T20, T01, ...) "
      "If not specified, the primitive cell is used for the unit "
      "cell. Ex: 'transf_mat', 'transf_mat(SCELV_A_B_C_D_E_F)'";

    TransfMat::TransfMat() :
      SupercellCheckBase(Name, Desc) {}

    Eigen::VectorXi TransfMat::evaluate(const Supercell &scel) const {
      // should be column-major anyways, but let's ensure it
      Eigen::Matrix<int, 3, 3, Eigen::ColMajor> T = std::get<2>(_evaluate(scel, *m_refcell));
      return Eigen::Map<Eigen::VectorXi>(T.data(), T.size());
    }

    bool TransfMat::validate(const Supercell &scel) const {
      // should be column-major anyways, but let's ensure it
      return std::get<0>(_evaluate(scel, *m_refcell));
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<TransfMat> TransfMat::clone() const {
      return std::unique_ptr<TransfMat>(this->_clone());
    }

    /// \brief Expects arguments of the form 'transf_mat(unitcell_name)'
    bool TransfMat::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

      if(splt_vec.size() > 1) {
        std::stringstream ss;
        ss << this->name() << " expected 0 or 1 argument.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }
      else if(args.empty()) {
        m_refcell_name = "SCEL1_1_1_1_0_0_0";
      }
      else {
        m_refcell_name = args;
      }
      return true;
    }

    /// \brief Clone using copy constructor
    TransfMat *TransfMat::_clone() const {
      return new TransfMat(*this);
    }


    // --- ConfigCountBase ---

    ConfigCountBase::ConfigCountBase(std::string name, std::string desc) :
      IntegerAttribute<Supercell>(name, desc) {}

    /// \brief Expects arguments of the form 'is_supercell_of(scelname)'
    bool ConfigCountBase::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

      if(splt_vec.size() > 1) {
        std::stringstream ss;
        ss << this->name() << " expected 0 or 1 argument.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }

      m_type = args;
      return true;
    }

    /// \brief col_header returns: {'name(refcell_name)'}
    std::vector<std::string> ConfigCountBase::col_header(const Supercell &_tmplt) const {
      return std::vector<std::string> {this->name() + "(" + m_type + ")"};
    }


    // --- Nconfig ---

    const std::string Nconfig::Name = "Nconfig";
    const std::string Nconfig::Desc =
      "Number of enumerated configurations of (default) all types, or specified "
      "type, by supercell. Ex: 'Nconfig', 'Nconfig(config)'";

    Nconfig::Nconfig() :
      ConfigCountBase(Name, Desc) {}

    Index Nconfig::evaluate(const Supercell &scel) const {
      if(m_type.empty()) {
        return DB::config_count(scel.name(), scel.primclex());
      }
      else {
        return DB::config_count(m_type, scel.name(), scel.primclex());
      }
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<Nconfig> Nconfig::clone() const {
      return std::unique_ptr<Nconfig>(this->_clone());
    }

    /// \brief Clone using copy constructor
    Nconfig *Nconfig::_clone() const {
      return new Nconfig(*this);
    }


    // --- Ncalc ---

    const std::string Ncalc::Name = "Ncalc";
    const std::string Ncalc::Desc =
      "Number of configurations with completed calculations of (default) all "
      "types, or specified type, by supercell. Ex: 'Ncalc', 'Ncalc(config)'";

    Ncalc::Ncalc() :
      ConfigCountBase(Name, Desc) {}

    Index Ncalc::evaluate(const Supercell &scel) const {
      if(m_type.empty()) {
        return DB::config_calculated_count(scel.name(), scel.primclex());
      }
      else {
        return DB::config_calculated_count(m_type, scel.name(), scel.primclex());
      }
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<Ncalc> Ncalc::clone() const {
      return std::unique_ptr<Ncalc>(this->_clone());
    }

    /// \brief Clone using copy constructor
    Ncalc *Ncalc::_clone() const {
      return new Ncalc(*this);
    }


    // --- Ndata ---

    const std::string Ndata::Name = "Ndata";
    const std::string Ndata::Desc =
      "Number of configurations of (default) all types, or specified "
      "type, which have any data or files, by supercell. Ex: 'Ndata', 'Ndata(config)'";

    Ndata::Ndata() :
      ConfigCountBase(Name, Desc) {}

    Index Ndata::evaluate(const Supercell &scel) const {
      if(m_type.empty()) {
        return DB::config_data_count(scel.name(), scel.primclex());
      }
      else {
        return DB::config_data_count(m_type, scel.name(), scel.primclex());
      }
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<Ndata> Ndata::clone() const {
      return std::unique_ptr<Ndata>(this->_clone());
    }

    /// \brief Clone using copy constructor
    Ndata *Ndata::_clone() const {
      return new Ndata(*this);
    }


    // --- GenericDatumFormatter generating functions ---

    GenericScelFormatter<std::string> pointgroup_name() {
      return GenericScelFormatter<std::string>(
               "pointgroup_name",
               "Supercell point group name.",
      [](const Supercell & scel)->std::string {
        return scel.factor_group().get_name();
      });
    }

    GenericScelFormatter<Index> scel_size() {
      return GenericScelFormatter<Index>(
               "scel_size",
               "Supercell volume, given as the integer number of primitive cells",
      [](const Supercell & scel)->Index {
        return scel.volume();
      });
    }

    GenericScelFormatter<Index> multiplicity() {
      return GenericScelFormatter<Index>(
               "multiplicity",
               "Number of equivalent supercells",
      [](const Supercell & scel)->Index {
        return scel.prim().factor_group().size() / scel.factor_group().size();
      });
    }

    GenericScelFormatter<Index> factorgroup_size() {
      return GenericScelFormatter<Index>(
               "factorgroup_size",
               "Supercell factor group size",
      [](const Supercell & scel)->Index {
        return scel.factor_group().size();
      });
    }

    GenericScelFormatter<double> volume() {
      return GenericScelFormatter<double>(
               "volume",
               "Supercell volume (length^3)",
      [](const Supercell & scel)->double {
        return scel.volume() * scel.lattice().vol();
      });
    }

    GenericVectorXdScelFormatter lattice() {
      return GenericVectorXdScelFormatter(
               "lattice",
               "Lattice vectors, unrolled: (a0, a1, a2, b0, ...)",
      [](const Supercell & scel)->Eigen::VectorXd {
        Eigen::Matrix<double, 3, 3, Eigen::ColMajor> L = scel.lattice().lat_column_mat();
        return Eigen::Map<Eigen::VectorXd>(L.data(), L.size());
      });
    }

    GenericVectorXdScelFormatter lattice_params() {
      return GenericVectorXdScelFormatter(
               "lattice_params",
               "Lattice parameters, as: (a, b, c, alpha, beta, gamma)",
      [](const Supercell & scel)->Eigen::VectorXd {
        Eigen::VectorXd res;
        res << scel.lattice().length(0), scel.lattice().length(1), scel.lattice().length(2),
            scel.lattice().angle(0), scel.lattice().angle(1), scel.lattice().angle(2);
        return res;
      });
    }

  }

  template<>
  StringAttributeDictionary<Supercell> make_string_dictionary<Supercell>() {

    using namespace ScelIO;
    StringAttributeDictionary<Supercell> dict;

    dict.insert(
      name<Supercell>(),
      alias<Supercell>(),
      alias_or_name<Supercell>(),
      pointgroup_name()
    );

    return dict;
  }

  template<>
  BooleanAttributeDictionary<Supercell> make_boolean_dictionary<Supercell>() {

    using namespace ScelIO;
    BooleanAttributeDictionary<Supercell> dict;

    dict.insert(
      DB::Selected<Supercell>(),
      IsSupercellOf(),
      IsUnitcellOf()
    );

    return dict;
  }

  template<>
  IntegerAttributeDictionary<Supercell> make_integer_dictionary<Supercell>() {

    using namespace ScelIO;
    IntegerAttributeDictionary<Supercell> dict;

    dict.insert(
      scel_size(),
      multiplicity(),
      factorgroup_size(),
      Nconfig(),
      Ncalc(),
      Ndata()
    );

    return dict;
  }

  template<>
  ScalarAttributeDictionary<Supercell> make_scalar_dictionary<Supercell>() {

    using namespace ScelIO;
    ScalarAttributeDictionary<Supercell> dict;

    dict.insert(
      volume()
    );

    return dict;
  }

  template<>
  VectorXiAttributeDictionary<Supercell> make_vectorxi_dictionary<Supercell>() {

    using namespace ScelIO;
    VectorXiAttributeDictionary<Supercell> dict;

    dict.insert(
      TransfMat()
    );

    return dict;
  }

  template<>
  VectorXdAttributeDictionary<Supercell> make_vectorxd_dictionary<Supercell>() {

    using namespace ScelIO;
    VectorXdAttributeDictionary<Supercell> dict;

    dict.insert(
      lattice(),
      lattice_params()
    );

    return dict;
  }

}
