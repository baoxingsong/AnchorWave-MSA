//
// Created by song on 8/5/18.
//

#ifndef PROALI_ALIGNSLIDINGWINDOW_H
#define PROALI_ALIGNSLIDINGWINDOW_H

#include <string>
#include <stack>
#ifdef __SSE2NEON__
#include "sse2neon.h"
#else
#include <immintrin.h>
#endif // __SSE2NEON__
#include <map>
#include "../impl/impl.h"
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include "../../WFA2-lib/bindings/cpp/WFAligner.hpp"
#include <stdlib.h>
#include "../../minimap2/ksw2.h"


int64_t alignSlidingWindow(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int64_t &slidingWindowSize, const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                           const Scorei &m);



int64_t alignSlidingWindowNW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                             const int64_t &slidingWindowSize, const int32_t &matchingScore,
                             const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                             const Scorei &m);



#endif //PROALI_ALIGNSLIDINGWINDOW_H
