#ifndef BP_StopWatch_CC
#define BP_StopWatch_CC

#include "casm/BP_C++/BP_StopWatch.hh"

namespace BP {

  BP_StopWatch::BP_StopWatch() {
    reset();
  };

  void BP_StopWatch::reset() {
    set_start();
    lap = start;
  };

  void BP_StopWatch::set_start() {
    start = gettime();
  };

  void BP_StopWatch::set_lap() {
    lap = gettime();
  };

  double BP_StopWatch::gettime() const {
    timeval tim;
    gettimeofday(&tim, NULL);
    return tim.tv_sec + (tim.tv_usec / 1000000.0);
  };

  double BP_StopWatch::gettime_ms() const {
    timeval tim;
    gettimeofday(&tim, NULL);
    return 1000.0 * (tim.tv_sec + (tim.tv_usec / 1000000.0));
  };

  double BP_StopWatch::gettime_us() const {
    timeval tim;
    gettimeofday(&tim, NULL);
    return 1000000.0 * (tim.tv_sec + (tim.tv_usec / 1000000.0));
  };

  double BP_StopWatch::get_us() const {
    timeval tim;
    gettimeofday(&tim, NULL);
    return tim.tv_usec;
  };

  double BP_StopWatch::get_s() const {
    timeval tim;
    gettimeofday(&tim, NULL);
    return tim.tv_sec;
  };

  double BP_StopWatch::total_time_s() const {
    return gettime() - start;
  };

  double BP_StopWatch::lap_time_s() {
    double new_lap = gettime();
    double lap_time = new_lap - lap;
    lap = new_lap;
    return lap_time;

  };

  /// return current date and time as string: YEAR-MM-DD HH:MM:SS
  std::string BP_StopWatch::date_time() const {
    time_t t = time(0);   // get time now
    struct tm *_now = localtime(& t);
    std::ostringstream ss;

    std::string smin = "";
    int min = _now->tm_min;
    if(min < 10)
      smin = "0" + BP::itos(min);
    else
      smin = BP::itos(min);

    std::string ssec = "";
    int sec = _now->tm_sec;
    if(sec < 10)
      ssec = "0" + BP::itos(sec);
    else
      ssec = BP::itos(sec);


    ss << (_now->tm_year + 1900) << '-'
       << (_now->tm_mon + 1) << '-'
       <<  _now->tm_mday
       << " " << _now->tm_hour << ":" << smin << ":" << ssec;
    return ss.str();

  };

}

#endif // BP_StopWatch_CC
