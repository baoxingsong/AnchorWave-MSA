//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
#define ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H

#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include <map>
#include <cstdlib>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "../../minimap2/minimap.h"
#include "../../minimap2/kseq.h"


void setupAnchorsWithSpliceAlignmentResult(const std::string &gffFilePath, const std::string &cdsSequenceFile, const std::string &samFile, std::map<std::string, std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                           double &inversion_PENALTY, double &MIN_ALIGNMENT_SCORE, bool &considerInversion, const int &minExon, const int64_t &windownWidth, const double &minimumSimilarity, const double &minimumSimilarity2,
                                           std::map<std::string, std::tuple<std::string, long, long, int> > &referenceGenome2,
                                           std::map<std::string, std::tuple<std::string, long, long, int> > &queryGenome2,
                                           int &expectedCopies, double &maximumSimilarity, const std::string &referenceSamFilePath, const int32_t &wfaSize3, const bool &searchForNewAnchors, const bool &exonModel);

void setupAnchorsWithSpliceAlignmentResultQuota(const std::string &gffFilePath, const std::string &samFile, const std::string &cdsSequenceFile, std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                                double &INDEL_SCORE, double &GAP_OPEN_PENALTY, double &MIN_ALIGNMENT_SCORE, int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                                                double &calculateIndelDistance, const int &minExon, const int64_t &windownWidth, const double &minimumSimilarity, const double &minimumSimilarity2,
                                                std::map<std::string, std::tuple<std::string, long, long, int> > &referenceGenome2,
                                                std::map<std::string, std::tuple<std::string, long, long, int> > &queryGenome2,
                                                int &expectedCopies, const int32_t &wfaSize3,
                                                double &maximumSimilarity, const std::string &referenceSamFilePath,
                                                bool &searchForNewAnchors, const bool &exonModel);

#endif //ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
