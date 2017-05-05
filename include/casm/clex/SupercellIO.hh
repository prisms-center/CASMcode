#ifndef CASM_SupercellIO
#define CASM_SupercellIO

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"


namespace CASM {
  class Supercell;
  class SymOp;

  namespace ScelIO {

    // --- BaseDatumFormatter-derived classes ---

    template<typename Base>
    class SupercellCheckBase : public Base {
    public:

      typedef std::tuple<bool, SymOp, Eigen::MatrixXi> result_type;

      SupercellCheckBase(std::string name, std::string desc);

      // --- Specialized implementation -----------

      /// \brief Expects arguments of the form 'is_supercell_of(scelname)'
      bool parse_args(const std::string &args) override;

      /// \brief Set pointer to ref supercell
      void init(const Supercell &_tmplt) const override;

      /// \brief col_header returns: {'short_name(refcell_name)'}
      std::vector<std::string> col_header(const Supercell &_tmplt) const override;

      // --- Shared logic -----------

      /// Call is_supercell using prim.factor_group() to try possible orientations
      ///
      /// Returns (bool, SymOp op, Eigen::MatrixXi T) with:
      /// - scel is supercell of unit?
      /// - If true: scel.lattice() == apply(op, unit.lattice()) * T
      ///
      const result_type &_evaluate(const Supercell &scel, const Supercell &unit) const;

    protected:

      /// Reference supercell name, given meaning by derived class
      std::string m_refcell_name;

      /// Reference supercell, given meaning by derived class
      mutable const Supercell *m_refcell;

    private:
      mutable notstd::cloneable_ptr<result_type> m_last_result;
      mutable const Supercell *m_last_scel;
      mutable const Supercell *m_last_unit;
    };

    class IsSupercellOf : public SupercellCheckBase<BooleanAttribute<Supercell> > {
    public:
      static const std::string Name;
      static const std::string Desc;

      IsSupercellOf();

      bool evaluate(const Supercell &scel) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<IsSupercellOf> clone() const;

    private:

      /// \brief Clone using copy constructor
      IsSupercellOf *_clone() const override;

    };

    class IsUnitcellOf : public SupercellCheckBase<BooleanAttribute<Supercell> > {
    public:
      static const std::string Name;
      static const std::string Desc;

      IsUnitcellOf();

      bool evaluate(const Supercell &unit) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<IsUnitcellOf> clone() const;

    private:

      /// \brief Clone using copy constructor
      IsUnitcellOf *_clone() const override;

    };

    class TransfMat : public SupercellCheckBase<VectorXiAttribute<Supercell> > {
    public:
      static const std::string Name;
      static const std::string Desc;

      TransfMat();

      Eigen::VectorXi evaluate(const Supercell &scel) const override;

      bool validate(const Supercell &scel) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<TransfMat> clone() const;

      /// \brief Expects arguments of the form 'transf_mat(unitcell_name)'
      bool parse_args(const std::string &args) override;

    private:

      /// \brief Clone using copy constructor
      TransfMat *_clone() const override;

    };

    class ConfigCountBase : public IntegerAttribute<Supercell> {
    public:

      ConfigCountBase(std::string name, std::string desc);

      // --- Specialized implementation -----------

      /// \brief Expects arguments of the form 'is_supercell_of(scelname)'
      bool parse_args(const std::string &args) override;

      /// \brief col_header returns: {'short_name(refcell_name)'}
      std::vector<std::string> col_header(const Supercell &_tmplt) const override;

    protected:

      /// Reference supercell name, given meaning by derived class
      std::string m_type;

    };

    class Nconfig : public ConfigCountBase {
    public:
      static const std::string Name;
      static const std::string Desc;

      Nconfig();

      Index evaluate(const Supercell &scel) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Nconfig> clone() const;

    private:

      /// \brief Clone using copy constructor
      Nconfig *_clone() const override;

    };

    class Ncalc : public ConfigCountBase {
    public:
      static const std::string Name;
      static const std::string Desc;

      Ncalc();

      Index evaluate(const Supercell &scel) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Ncalc> clone() const;

    private:

      /// \brief Clone using copy constructor
      Ncalc *_clone() const override;

    };

    class Ndata : public ConfigCountBase {
    public:
      static const std::string Name;
      static const std::string Desc;

      Ndata();

      Index evaluate(const Supercell &scel) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Ndata> clone() const;

    private:

      /// \brief Clone using copy constructor
      Ndata *_clone() const override;

    };

    // --- GenericDatumFormatter generating functions ---

    template<typename ValueType>
    using GenericScelFormatter = GenericDatumFormatter<ValueType, Supercell>;

    typedef Generic1DDatumFormatter<Eigen::VectorXd, Supercell> GenericVectorXdScelFormatter;

    GenericScelFormatter<std::string> name();
    GenericScelFormatter<std::string> alias();
    GenericScelFormatter<std::string> name_or_alias();
    GenericScelFormatter<std::string> scelname();
    GenericScelFormatter<std::string> pointgroup_name();
    GenericScelFormatter<Index> scel_size();
    GenericScelFormatter<Index> multiplicity();
    GenericScelFormatter<Index> factorgroup_size();
    GenericScelFormatter<double> volume();
    GenericVectorXdScelFormatter lattice();
    GenericVectorXdScelFormatter lattice_params();

  }

  template<>
  StringAttributeDictionary<Supercell> make_string_dictionary<Supercell>();

  template<>
  BooleanAttributeDictionary<Supercell> make_boolean_dictionary<Supercell>();

  template<>
  IntegerAttributeDictionary<Supercell> make_integer_dictionary<Supercell>();

  template<>
  ScalarAttributeDictionary<Supercell> make_scalar_dictionary<Supercell>();

  template<>
  VectorXiAttributeDictionary<Supercell> make_vectorxi_dictionary<Supercell>();

  template<>
  VectorXdAttributeDictionary<Supercell> make_vectorxd_dictionary<Supercell>();

}

#endif
