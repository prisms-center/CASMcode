#ifndef BP_StopWatch_HH
#define BP_StopWatch_HH

#include<sys/time.h>
#include <cstddef>
#include <sstream>
#include "casm/BP_C++/BP_Parse.hh"

// //////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////
//  BP_StopWatch class

namespace BP {

  /// \ingroup BP_useful
  class BP_StopWatch {
  private:
    double start;	//seconds
    double lap;		//seconds

  public:
    BP_StopWatch();

    void reset();

    void set_start();

    void set_lap();

    double gettime() const;

    double gettime_ms() const;

    double gettime_us() const;

    double get_us() const;

    double get_s() const;

    double total_time_s() const;

    double lap_time_s();

    /// return current date and time as string: YEAR-MM-DD HH:MM:SS
    std::string date_time() const;

  };

}

#endif // BP_StopWatch_HH
