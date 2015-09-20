
namespace CASM {
  template<typename ConfigType>
  ConfigEnumAllOccupations<ConfigType>::ConfigEnumAllOccupations(const ConfigEnumAllOccupations<ConfigType>::value_type &_initial,
                                                                 const ConfigEnumAllOccupations<ConfigType>::value_type &_final,
                                                                 PermuteIterator _perm_begin, PermuteIterator _perm_end) :
    ConfigEnum<ConfigType>(_initial, _final, -1),
    m_counter(_initial.occupation(), _final.occupation(), Array<int>(_initial.size(), 1)),
    m_perm_begin(_perm_begin), m_perm_end(_perm_end) {
    //std::cout << "INITIALIZING OCCUPATION ENUMERATOR\n";
    //std::cout << "starting at: "<< current().occupation() << "\n";
    //std::cout << "running to: " << final().occupation() << "\n";
    //std::cout << "current().is_primitive()=" << current().is_primitive(_perm_begin) << " and current().is_canonical()=" << current().is_canonical(_perm_begin, _perm_end)<<"\n";

    // set source to describe enumeration procedure -- just the algorithm description for now
    // TODO: Add information about 'seed' configuration (_perm_begin)?
    _source() = "occupation_enumeration";

    // Make sure that current() has primitive canonical config
    if(!(current().is_primitive(_perm_begin) && current().is_canonical(_perm_begin, _perm_end))) {
      //std::cout << "INITIAL ENUMERATION STATE IS NOT CANONICAL!\n";
      increment();
    }

    if(!m_counter.valid())
      _step() = -1;
    else
      _step() = 0;

  }

  //*******************************************************************************************
  // **** Mutators ****
  // increment m_current and return a reference to it
  template<typename ConfigType>
  const typename ConfigEnumAllOccupations<ConfigType>::value_type &ConfigEnumAllOccupations<ConfigType>::increment() {
    bool is_valid_config(false);
    //std::cout << "Incrementing...\n";
    while(!is_valid_config && ++m_counter) {

      _current().set_occupation(m_counter());
      is_valid_config = current().is_primitive(_perm_begin()) && current().is_canonical(_perm_begin(), _perm_end());
      //std::cout << "counter() is: " << m_counter() << ";  is_valid_config: " << is_valid_config
      //<< ";  is_valid_counter: " << m_counter.valid() << "\n";
    }

    if(m_counter.valid())
      _step()++;
    else {
      //std::cout << "REACHED END OF THE LINE!\n";
      _step() = -1;
    }
    //std::cout << "--FINISHED SEARCH " << _step()<< "--\n";
    return current();
  }

  //*******************************************************************************************
  // set m_current to correct value at specified step and return a reference to it
  template<typename ConfigType>
  const typename ConfigEnumAllOccupations<ConfigType>::value_type &ConfigEnumAllOccupations<ConfigType>::goto_step(step_type _step) {
    std::cerr << "CRITICAL ERROR: Class ConfigEnumAllOccupations does not implement a goto_step() method. \n"
              << "                You may be using a ConfigEnumIterator in an unsafe way!\n"
              << "                Exiting...\n";
    assert(0);
    exit(1);
    return current();
  };
}

