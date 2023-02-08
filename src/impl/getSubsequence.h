//
// Created by baoxing on 10/10/17.
//

#pragma once

#include "../model/STRAND.h"
#include "GetReverseComplementary.h"
#include <unistd.h>
#include <fcntl.h>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <set>

std::string getSubsequence2(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName, const int &_start, const int &_end);

std::string getSubsequence3(std::map<std::string, std::tuple<std::string, long, long, int> > &map, int &fd, const std::string &seqName, const int &_start, const int &_end);

char getCharByPos(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName, const int &_pos);

std::string getSubsequence2(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName);

std::string getSubsequence3(std::map<std::string, std::tuple<std::string, long, long, int> > &map, int &fd, const std::string &seqName);

std::string getSubsequence3(std::map<std::string, std::tuple<std::string, long, long, int> > &map, int &fd, const std::string &seqName, const int &_start, const int &_end, const STRAND &strand);

std::string getSubsequence2(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName, const int &_start, const int &_end, const STRAND &strand);

size_t getSequenceSizeFromPath2(std::tuple<std::string, long, long, int> &t);

std::string getSubsequence(const std::string &sequence, const int &_start, const int &_end);

std::string getSubsequence(const std::string &sequence, const int &_start, const int &_end, const STRAND &strand);
