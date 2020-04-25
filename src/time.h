#ifndef TIME_H
#define TIME_H

#include <string>
#include <iostream>

char *strptime(const char * __restrict, const char * __restrict, struct tm * __restrict);

class Time {
private:
    int _seconds;

    static const int SECONDS_IN_DAY = 86400;
    static const int SECONDS_IN_HOUR = 3600;
    static const int SECONDS_IN_MIN = 60;

public:
    Time ();
    Time (int s);
    Time (int h, int m, int s);
    Time (std::string& t);
    Time (uint64_t& t);

    // let you call Time::now () anywhere
    static Time now ();

    int hour () const;
    int minute () const;
    int second () const;
    int seconds () const;

    uint64_t asUNIX (uint64_t& day) const;
};

inline int operator-(const Time& lhs, const Time& rhs)
{
    return lhs.seconds () - rhs.seconds ();
}

inline int operator+(const Time& lhs, const Time& rhs)
{
    return lhs.seconds () + rhs.seconds ();
}

inline int operator-(const Time& lhs, const int& rhs)
{
    return lhs.seconds () - rhs;
}

inline int operator+(const Time& lhs, const int& rhs)
{
    return lhs.seconds () + rhs;
}

inline bool operator==(const Time& lhs, const Time& rhs)
{
    return lhs.seconds () == rhs.seconds ();
}
inline bool operator!=(const Time& lhs, const Time& rhs)
{
    return !(lhs == rhs);
}
inline bool operator<(const Time& lhs, const Time& rhs)
{
    return lhs.seconds () < rhs.seconds ();
}
inline bool operator>(const Time& lhs, const Time& rhs)
{
    return rhs < lhs;
}
inline bool operator<=(const Time& lhs, const Time& rhs)
{
    return !(lhs > rhs);
}
inline bool operator>=(const Time& lhs, const Time& rhs)
{
    return !(lhs < rhs);
}

inline std::ostream& operator<<(std::ostream& os, const Time& t)
{
    os
        << (t.hour () < 10 ? "0" : "")
        << t.hour () << ":"
        << (t.minute () < 10 ? "0" : "")
        << t.minute () << ":"
        << (t.second () < 10 ? "0" : "")
        << t.second ();
    return os;
}

#endif
