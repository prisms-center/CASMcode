#ifndef CONFIGENUM_HH
#define CONFIGENUM_HH

#include "casm/casm_io/jsonParser.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  template <typename ConfigType>
  class ConfigEnumIterator;

  template <typename ConfigType>
  class ConfigEnum {
  public:
    typedef long step_type;

    // ConfigType is either Configuration or ConfigDoF
    typedef ConfigType value_type;

    typedef ConfigEnumIterator<ConfigType> iterator;


  private:
    value_type m_initial, m_final;
    value_type m_current;
    step_type m_step, m_N_steps;

    /// member to store 'source' info that will be used to tag enumerated configurations
    jsonParser m_source;
  public:
    ConfigEnum() {};

    ConfigEnum(const value_type &_initial, const value_type &_final, Index _N_steps) :
      m_initial(_initial), m_final(_final), m_current(_initial), m_step(0), m_N_steps(_N_steps) {

    };

    void initialize(const value_type &_initial, const value_type &_final, Index _N_steps) {
      m_initial = _initial;
      m_final = _final;
      m_current = _initial;
      m_step = 0;
      m_N_steps = _N_steps;
    }

    virtual ~ConfigEnum() {};

    const value_type &current() const {
      return m_current;
    }
    const value_type &initial() const {
      return m_initial;
    }
    const value_type &final() const {
      return m_final;
    }
    const jsonParser &source() const {
      return m_source;
    }

    step_type step() const {
      return m_step;
    };

    // Number of steps separating m_initial from m_final, with m_final being a num_steps()-1
    step_type num_steps() const {
      return m_N_steps;
    };

    iterator begin() {
      return iterator(*this, 0);
    }

    iterator end() {
      return iterator(*this, num_steps());
    }

    // **** Mutators ****
    // increment m_current and return a reference to it
    virtual const value_type &increment() = 0;

    // set m_current to correct value at specified step and return a reference to it
    virtual const value_type &goto_step(step_type _step) = 0;

  protected:
    value_type &_current() {
      return m_current;
    }
    value_type &_initial() {
      return m_initial;
    }
    value_type &_final() {
      return m_final;
    }

    jsonParser &_source() {
      return m_source;
    }

    step_type &_step() {
      return m_step;
    };

    // Number of steps separating m_initial from m_final, with m_final being a num_steps()-1
    step_type &_num_steps() {
      return m_N_steps;
    };



  };


}
#endif
