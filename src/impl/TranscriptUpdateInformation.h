//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_TRANSCRIPTUPDATEINFORMATION_H
#define ANNOTATIONLIFTOVER_TRANSCRIPTUPDATEINFORMATION_H

#include "../model/Transcript.h"
#include "getSubsequence.h"

void TranscriptUpdateCdsInformation(Transcript &transcript, std::map<std::string, std::tuple<std::string, long, long, int> > &genome);

#endif //ANNOTATIONLIFTOVER_TRANSCRIPTUPDATEINFORMATION_H
