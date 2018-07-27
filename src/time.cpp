#include "time.h"
#include <sstream>
#include <iomanip>
#include <ctime>
#include <stdexcept>

#include <Rcpp.h>

Time::Time ()
{
    _seconds = 0;
}

Time::Time (int s)
{
    if (s < 0)
    {
        throw std::invalid_argument ("Seconds must be 0 or greater.");
    }
    _seconds = s;
}

Time::Time (int h, int m, int s)
{
    if (h < 0)
    {
        throw std::invalid_argument ("Hours must be 0 or greater.");
    }
    if (m < 0 || m >= 60)
    {
        throw std::invalid_argument ("Minutes must be between 0 and 60.");
    }
    if (s < 0 || s >= 60)
    {
        throw std::invalid_argument ("Seconds must be between 0 and 60.");
    }

    int sec = h * SECONDS_IN_HOUR + m * SECONDS_IN_MIN + s;
    _seconds = sec;
}

Time::Time (std::string& t)
{
    struct std::tm tm;
    std::istringstream input (t);
    input >> std::get_time (&tm, "%H:%M:%S");
    if (input.fail ())
    {
        throw std::runtime_error ("Unable to read the time :(");
    }
    int sec = tm.tm_hour * SECONDS_IN_HOUR + 
        tm.tm_min * SECONDS_IN_MIN + tm.tm_sec;
    _seconds = sec;
}

Time Time::now ()
{
    time_t t = time (0);
    struct tm * now = localtime (&t);
    return Time (now->tm_hour, now->tm_min, now->tm_sec);
}

int Time::hour () const
{
    return _seconds / SECONDS_IN_HOUR;
}

int Time::minute () const
{
    return (_seconds % SECONDS_IN_HOUR) / SECONDS_IN_MIN;
}

int Time::second () const
{
    return _seconds % SECONDS_IN_MIN;
}

int Time::seconds () const
{
    return _seconds;
}

void Time::print () const
{
    Rcpp::Rcout << _seconds;
}