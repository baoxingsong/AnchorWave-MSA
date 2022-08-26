//
// Created by zhf on 22-8-23.
//

#ifndef ANCHORWAVE_MUTISEQUENCESCOMPARE_H
#define ANCHORWAVE_MUTISEQUENCESCOMPARE_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <stack>
#include "./sequencesToKmerDIstanceMatrix.h"
#include "../myImportandFunction/myImportantFunction.h"
#include "../util/util.h"
void mutiSequencesCompre(std::map<std::string, std::string> &sequences, std::vector<std::string> &seqNames, int8_t kmer_size, int32_t loopCountend,
                         int32_t mismatchingPenalty, int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                         int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2, std::string& outResultPath);
#endif //ANCHORWAVE_MUTISEQUENCESCOMPARE_H

