//
// Created by baoxing on 10/10/17.
//

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <unistd.h>
#include <sys/time.h>
#include <algorithm>

void readFastaFile(const std::string &filePath, std::map<std::string, std::tuple<std::string, long, long, int> > &map);

void readFastaFileWorkWithIUPACcode(const std::string &filePath, std::map<std::string, std::string> &sequences);
