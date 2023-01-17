// =====================================================================================
// 
//       Filename:  myutil.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 04:53:23 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//
// =====================================================================================



#include "myutil.h"

#ifdef __unix

#include <linux/limits.h>

#endif

void split(const std::string &s, char &delim, std::vector <std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (item.length() > 0) {
            elems.push_back(item);
        }
    }
}

std::string songStrReplaceAll(std::string &str, const std::string &pattern, const std::string &pattern2) {
    std::regex vowel_re(pattern);
    str = std::regex_replace(str, vowel_re, pattern2);
    return str;
}

int32_t min(const int32_t &a, const int32_t &b) {
    if (a < b) {
        return a;
    }

    return b;
}

void splitCIGAR(std::string str, std::vector <std::string> &cigarElems) {
    while (str.length() > 0) {
        size_t p_d = str.find_first_of("1234567890");
        size_t p_c = str.find_first_of("MIDNSHPX=");

        if(p_c > p_d) {
            std::string s = str.substr(p_d, p_c + 1);
            cigarElems.push_back(s);
        }
        str = str.substr(p_c + 1);
    }
}
