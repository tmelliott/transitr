#include "time.h"
#include <ctime>
#include <stdexcept>

#include <Rcpp.h>

Time::Time ()
{
    _seconds = 0;
}

Time::Time (int s)
{
    _seconds = 0;
    if (s < 0)
    {
        return;
        // throw std::invalid_argument ("Seconds must be 0 or greater.");
    }
    _seconds = s;
}

Time::Time (int h, int m, int s)
{
    _seconds = 0;

    if (h < 0)
    {
        return;
        // throw std::invalid_argument ("Hours must be 0 or greater.");
    }
    if (m < 0 || m >= 60)
    {
        return;
        // throw std::invalid_argument ("Minutes must be between 0 and 60.");
    }
    if (s < 0 || s >= 60)
    {
        return;
        // throw std::invalid_argument ("Seconds must be between 0 and 60.");
    }

    int sec = h * SECONDS_IN_HOUR + m * SECONDS_IN_MIN + s;
    _seconds = sec;
}

Time::Time (std::string& t)
{
    if (t.length () != 8)
    {
        throw std::invalid_argument ("Time should be in the form HH:MM:SS");
    }
// #ifdef __linux__
//     struct tm tm;
//     strptime (t.c_str (), "%H:%M:%S", &tm);
//     int sec = tm.tm_hour * SECONDS_IN_HOUR +
//         tm.tm_min * SECONDS_IN_MIN + tm.tm_sec;
// #else
    // necessary to manually pull apart the string
    std::string delim = ":";
    size_t pos = 0;
    std::string token;
    std::vector<int> tmv;
    while ((pos = t.find (delim)) != std::string::npos)
    {
        token = t.substr (0, pos);
        tmv.push_back (stoi (token));
        t.erase (0, pos + delim.length ());
    }
    if (tmv.size () != 2)
    {
        throw std::invalid_argument ("Time string should be of the form HH:MM:SS");
    }
    // should leave the seconds as the remaining characters in `t`
    int sec = tmv[0] * SECONDS_IN_HOUR +
        tmv[1] * SECONDS_IN_MIN + stoi (t);
// #endif

    _seconds = sec;
}

Time::Time (uint64_t& t)
{
    // convert UNIX timestamp to a time object
    if (t < 0)
    {
        _seconds = 0;
        return;
    }
    std::time_t tx (t);
    struct tm *t0 = localtime(&tx);
    _seconds = t0->tm_hour * SECONDS_IN_HOUR +
        t0->tm_min * SECONDS_IN_MIN + t0->tm_sec;
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

uint64_t Time::asUNIX (uint64_t& day) const
{
    // get midnight of 'day'
    Time today (day);
    uint64_t t0 (day - today.seconds ());
    return t0 + _seconds;
}
