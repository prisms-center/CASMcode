#include "casm/global/enum/io_traits.hh"

namespace CASM {

const std::string traits<COORD_TYPE>::name = "coordinate_mode";

/// Excludes COORD_DEFAULT, which isn't used for IO
const std::multimap<COORD_TYPE, std::vector<std::string> >
    traits<COORD_TYPE>::strval = {
        {COORD_TYPE::FRAC,
         {"Fractional", "fractional", "Direct", "direct", "FRAC"}},
        {COORD_TYPE::CART, {"Cartesian", "cartesian", "CART"}},
        {COORD_TYPE::INTEGRAL, {"Integral", "integral", "INTEGRAL"}}};

const std::string traits<PERIODICITY_TYPE>::name = "periodicity_type";

/// Excludes PERIODICITY_DEFAULT, which isn't used for IO
const std::multimap<PERIODICITY_TYPE, std::vector<std::string> >
    traits<PERIODICITY_TYPE>::strval = {
        {PERIODICITY_TYPE::PERIODIC, {"PERIODIC"}},
        {PERIODICITY_TYPE::APERIODIC, {"APERIODIC", "LOCAL"}}};

const std::string traits<EQUIVALENCE_TYPE>::name = "equivalence_type";

const std::multimap<EQUIVALENCE_TYPE, std::vector<std::string> >
    traits<EQUIVALENCE_TYPE>::strval = {
        {EQUIVALENCE_TYPE::PRIM, {"PRIM", "prim"}},
        {EQUIVALENCE_TYPE::SCEL, {"SCEL", "scel"}},
        {EQUIVALENCE_TYPE::CONFIG, {"CONFIG", "config"}}};

const std::string traits<CELL_TYPE>::name = "cell_type";

const std::multimap<CELL_TYPE, std::vector<std::string> >
    traits<CELL_TYPE>::strval = {{CELL_TYPE::PRIM, {"PRIM"}},
                                 {CELL_TYPE::SCEL, {"SCEL"}}};

const std::string traits<OnError>::name = "on_error";

const std::multimap<OnError, std::vector<std::string> >
    traits<OnError>::strval = {{OnError::THROW, {"THROW", "throw"}},
                               {OnError::WARN, {"WARN", "warn"}},
                               {OnError::CONTINUE, {"CONTINUE", "continue"}}};

}  // namespace CASM
