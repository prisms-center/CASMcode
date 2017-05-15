#ifndef CASM_Selected
#define CASM_Selected

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/database/Selection.hh"

namespace CASM {

  namespace DB {

    /// \brief Returns true if configuration is specified in given selection (default: MASTER)
    ///
    /// Ex: 'selected_in', 'selected_in(myselection.txt)'
    ///
    /// \ingroup ConfigIO
    ///
    template<typename ObjType>
    class Selected : public BooleanAttribute<ObjType> {
    public:
      Selected():
        BooleanAttribute<ObjType>(
          "selected_in",
          "Returns true if object is specified in given selection (default MASTER). "
          "Ex: 'selected_in(myselection.txt)'") {}

      explicit Selected(const Selection<ObjType> &_selection):
        BooleanAttribute<ObjType>(
          "selected_in",
          "Returns true if configuration is specified in given selection (default MASTER). "
          "Ex: 'selected_in(myselection.txt)'"),
        m_selection(notstd::clone(_selection)) {}

      explicit Selected(std::unique_ptr<Selection<ObjType> > _selection):
        BooleanAttribute<ObjType>(
          "selected_in",
          "Returns true if configuration is specified in given selection (default MASTER). "
          "Ex: 'selected_in(myselection.txt)'"),
        m_selection(std::move(_selection)) {}

      // --- Required implementations -----------

      std::unique_ptr<Selected> clone() const;

      /// \brief True if obj is in Selection and it is selected
      bool evaluate(const ObjType &_obj) const override;


      // --- Specialized implementation -----------

      void init(const ObjType &_tmplt) const override;

      std::string short_header(const ObjType &_config) const override;

      bool parse_args(const std::string &args) override;


    private:

      /// \brief Clone
      Selected *_clone() const override;

      mutable std::string m_selection_name;
      mutable notstd::cloneable_ptr<Selection<ObjType> > m_selection;
    };
  }
}

#endif
