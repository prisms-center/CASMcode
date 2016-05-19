#ifndef CASM_DoFMod_HH
#define CASM_DoFMod_HH

namespace CASM {

  /// \brief Describes the modification of a variable on a basis site
  ///
  /// SiteMod has information about the occupational change on a single site. It holds only
  /// an index corresponding to the site in the ConfigDoF and the final value a variable
  /// on that site will have after the modification takes place.
  ///
  template<class T>
  class SiteMod {
  public:

    typedef Index size_type;

    /// \brief Default constructor
    SiteMod() {}

    /// \brief Construct a SiteMod
    ///
    /// \param _site_linear_index Index into the unrolled array of sites in the ConfigDoF of site being modified
    /// \param _sublat Sublattice index of site being modified
    /// \param _to_value Value the variable on the site being modified will change to
    ///
    SiteMod(size_type _site_linear_index, size_type _sublat, const T &_to_value) {
      set(_site_linear_index, _sublat, _to_value);
    }

    /// \brief Returns the linear index corresponding to site in ConfigDoF
    size_type site_index() const {
      return m_site_index;
    }

    /// \brief Returns the sublattice index of site being modified
    size_type sublat() const {
      return m_sublat;
    }


    /// \brief Returns the value the variable on the site being modified will change to
    const T &to_value() const {
      return m_to_value;
    }

    /// \brief Set the values of a SiteMod
    ///
    /// \param _site_linear_index Index into the unrolled array of sites in the ConfigDoF of site being modified
    /// \param _sublat Sublattice index of site being modified
    /// \param _to_value Value the variable on the site being modified will change to
    ///
    void set(size_type _site_linear_index, size_type _sublat, const T &_to_value) {
      m_site_index = _site_linear_index;
      m_sublat = _sublat;
      m_to_value = _to_value;
      return;
    }

  private:

    size_type m_site_index;
    size_type m_sublat;
    T m_to_value;
  };

  /// \brief An OccMod describes the change in occupation variable on a site
  typedef SiteMod<int> OccMod;

}

#endif
