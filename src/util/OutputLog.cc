//
// Created by Bernardo Clavijo (TGAC) on 31/12/2016.
//

#include "OutputLog.h"

template <class cT, class traits = std::char_traits<cT> >
class basic_nullbuf: public std::basic_streambuf<cT, traits> {
    typename traits::int_type overflow(typename traits::int_type c)
    {
        return traits::not_eof(c); // indicate success
    }
};

template <class cT, class traits = std::char_traits<cT> >
class basic_onullstream: public std::basic_ostream<cT, traits> {
public:
    basic_onullstream():
            std::basic_ios<cT, traits>(&m_sbuf),
            std::basic_ostream<cT, traits>(&m_sbuf)
    {
        init(&m_sbuf);
    }

private:
    basic_nullbuf<cT, traits> m_sbuf;
};

typedef basic_onullstream<char> onullstream;


onullstream ons;
unsigned OutputLogLevel;

std::ostream & OutputLog(const unsigned level,const bool include_date){
    if (OutputLogLevel>=level){
        if (include_date) {
            std::time_t t = std::time(NULL);
            char date_str[20];
            std::strftime(date_str, 20, "%Y-%m-%d %H:%M:%S", std::localtime(&t));
            std::cout << date_str << ": ";
        }
        return std::cout;
    }
    return ons;
}

