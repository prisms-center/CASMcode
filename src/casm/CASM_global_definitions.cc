#include "casm/CASM_global_definitions.hh"

namespace CASM {

  bool valid_index(Index i) {
    return 0 <= i;
  }


  void print_splash(std::ostream &out) {

    out << "       .::::::::.        .:::::.         .::::::.          .:.       .:.     \n"
        << "     .:::::::::::.      .:::::::.      .::::::::::.       .:::.     .:::.    \n"
        << "    .:::'     ':::.    .:::' ':::.    .:::'   '::::      .:::::.   .:::::.   \n"
        << "    ::::       ::::   .:::'   ':::.   ::::     '::'     .:::::::. .:::::::.  \n"
        << "    ::::              ::::     ::::   '::::.            ::::'':::.:::''::::  \n"
        << "    ::::              ::::     ::::    '::::::.         ::::  ':::::'  ::::  \n"
        << "    ::::              :::::::::::::      ''::::::.      ::::   ':::'   ::::  \n"
        << "    ::::              :::::::::::::         '::::::     ::::    ':'    ::::  \n"
        << "    ::::      .::::   ::::     ::::            :::::    ::::           ::::  \n"
        << "    ':::.    .::::'   ::::     ::::   .::.     .::::    ::::           ::::  \n"
        << "     '::::::::::'     ::::     ::::   :::::...:::::'    ::::           ::::  \n"
        << "        ':::::'       '::'     '::'   ':::::::::::'     '::'           '::'  \n";
  }

}

namespace CASM {

  const std::string traits<COORD_TYPE>::name = "coordinate_mode";

  const std::multimap<COORD_TYPE, std::vector<std::string> > traits<COORD_TYPE>::strval = {
    {COORD_TYPE::FRAC, {"Direct", "direct", "Fractional", "fractional", "FRAC"} },
    {COORD_TYPE::CART, {"Cartesian", "cartesian", "CART"} },
    {COORD_TYPE::INTEGRAL, {"Integral", "integral", "INTEGRAL"} },
    {COORD_TYPE::COORD_DEFAULT, {"COORD_DEFAULT"} }
  };

  const std::string traits<PERIODICITY_TYPE>::name = "periodicity_type";

  const std::multimap<PERIODICITY_TYPE, std::vector<std::string> > traits<PERIODICITY_TYPE>::strval = {
    {PERIODICITY_TYPE::PERIODIC, {"PERIODIC"} },
    {PERIODICITY_TYPE::LOCAL, {"LOCAL"} },
    {PERIODICITY_TYPE::PERIODICITY_DEFAULT, {"PERIODICITY_DEFAULT"} }
  };

  const std::string traits<CELL_TYPE>::name = "cell_type";

  const std::multimap<CELL_TYPE, std::vector<std::string> > traits<CELL_TYPE>::strval = {
    {CELL_TYPE::PRIM, {"PRIM"} },
    {CELL_TYPE::SCEL, {"SCEL"} }
  };

}


