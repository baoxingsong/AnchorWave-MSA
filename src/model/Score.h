/*
 * =====================================================================================
 *
 *       Filename:  Score.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  04/04/2019 19:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************


Reading the score from configure and provide score query functions

 ************************************************************************/

#pragma once

#include <map>
#include <iostream>
#include <fstream>

class Score {
private:
    std::map<int16_t, int32_t **> m;
    std::map<int16_t, int32_t> openGapPenalty1;
    std::map<int16_t, int32_t> extendGapPenalty1;
    std::map<int16_t, int32_t> openGapPenalty2;
    std::map<int16_t, int32_t> extendGapPenalty2;
    std::map<int16_t, int32_t> zdrop;
public:
    Score(const std::string &folder);

    ~Score();

    int32_t **getM(const int16_t &category);

    int32_t &getOpenGapPenalty1(const int16_t &category);

    int32_t &getExtendGapPenalty1(const int16_t &category);

    int32_t &getOpenGapPenalty2(const int16_t &category);

    int32_t &getExtendGapPenalty2(const int16_t &category);

    int32_t &getZdrop(const int16_t &category);
};


class Scorei {
private:
    int8_t **m;
public:
    Scorei(const int8_t &matchingScore, const int8_t &mismatchingPenalty);

    ~Scorei();

    const int8_t &getScore(const int8_t &a, const int8_t &b) const;

    int8_t **getScore() const;
};
