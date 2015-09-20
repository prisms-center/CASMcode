#ifndef GEOMETRYPIECES_HH
#define GEOMETRYPIECES_HH

namespace BP {
  class Geo;
  template<class T> class BP_Vec;
}

namespace CASM {
  class jsonParser;
  template<class T> class Array;

  ///Collect all the interesting things from your hull and give it back in json format
  jsonParser hull_data(BP::Geo &hull);
  ///Collect all the interesting things from a facet on the hull and give it back in json format
  jsonParser facet_data(BP::Geo &hull, int index);
  ///Collect all the interesting things from a index on the hull and give it back in json format
  jsonParser vertex_data(BP::Geo &hull, int index);

  ///This is horrible and should never be used
  Array<int> BP_to_CASM(const BP::BP_Vec<int> &value);
}
#endif
