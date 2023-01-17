//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_READFASTAFILE_H
#define ANNOTATIONLIFTOVER_READFASTAFILE_H

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include "../model/model.h"
#include <sstream>
#include <unistd.h>
#include <sys/time.h>

void readFastaFile(const std::string &filePath, std::map<std::string, std::tuple<std::string, long, long, int> > &map);

void readFastaFileWorkWithIUPACcode(const std::string &filePath, std::map<std::string, std::string> &sequences);

#endif //ANNOTATIONLIFTOVER_READFASTAFILE_H
