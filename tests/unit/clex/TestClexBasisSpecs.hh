#ifndef CASMtest_TestClexBasisSpecs
#define CASMtest_TestClexBasisSpecs

#include <string>

namespace test {

  // site_basis_functions : "occupation"
  // cluster_specs:
  //   method: "periodic_max_length"
  //   orbit_branch_specs:
  //     2: max_length=4.01
  //     3: max_length=3.01
  //   orbit_specs: two custom quadruplet clusters
  std::string FCC_ternary_clex_basis_specs_str_ex0();

  inline std::string FCC_ternary_clex_basis_specs_str_ex0() {
    return R"({
"basis_function_specs" : {
"dof_specs": {
  "occ": {
    "site_basis_functions" : "occupation"
  }
}
},
"cluster_specs": {
"method": "periodic_max_length",
"params": {
  "orbit_branch_specs" : {
    "2" : {"max_length" : 4.01},
    "3" : {"max_length" : 3.01}
  },
  "orbit_specs" : [
    {
      "coordinate_mode" : "Direct",
      "prototype" : [
        [ 0.000000000000, 0.000000000000, 0.000000000000 ],
        [ 1.000000000000, 0.000000000000, 0.000000000000 ],
        [ 2.000000000000, 0.000000000000, 0.000000000000 ],
        [ 3.000000000000, 0.000000000000, 0.000000000000 ]
      ],
      "include_subclusters" : true
    },
    {
      "coordinate_mode" : "Direct",
      "prototype" : [
        [ 0.000000000000, 0.000000000000, 0.000000000000 ],
        [ 0.000000000000, 1.000000000000, 0.000000000000 ],
        [ 0.000000000000, 0.000000000000, 1.000000000000 ],
        [ 1.000000000000, 1.000000000000, 1.000000000000 ]
      ],
      "include_subclusters" : true
    }
  ]
}
}
})";
  }
}

#endif
