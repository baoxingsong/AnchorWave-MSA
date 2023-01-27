//
// Created by Baoxing song on 20.10.18.
//

#ifndef PROALI_DENOVOGENOMEVARIANTCALLING_H
#define PROALI_DENOVOGENOMEVARIANTCALLING_H

#include <ctime>
#include "../model/model.h"
#include "../util/util.h"
#include "./readFastaFile.h"
#include "../myImportandFunction/myImportantFunction.h"
#include <iomanip>

#include <atomic>
#include <mutex>
#include <unistd.h>
#include <thread>
#include <iostream>
#include <chrono>


void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                      const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                                      const int32_t &widownWidth, /*const int32_t &wfaSize, const int32_t &wfaSize2, */const std::string &outPutMafFile,
                                      const std::string &outPutFragedFile, const int32_t &matchingScore,
                                      const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1,
                                      const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                      const int &maxThread);

void genomeAlignment(std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                     const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                     const int32_t &widownWidth, /*const int32_t &wfaSize, const int32_t &wfaSize2,*/
                     const std::string &outPutMafFile, const std::string &outPutFragedFile,
                     const int32_t &matchingScore, const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                     const int32_t &extendGapPenalty1,
                     const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                     const int &maxThread);

#endif //PROALI_DENOVOGENOMEVARIANTCALLING_H
