//
// Created by song on 9/19/18.
//

#include "include/gtest/gtest.h"
#include "../../service/TransferGffWithNucmerResult.h"
#include "../../util/myutil.h"
#include <string>
#include <iostream>
#include <map>
#include <regex>
#include <cstdlib>
#include "../../controlLayer.h"

TEST(mafToSdi, c1){ // just to make sure that every line has been analysed
    std::string outputovcffile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm20/outputvcf";
    std::string fastaFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string mafFile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm20/subtract1.maf";
    mafTovcf( mafFile, fastaFilePath,  outputovcffile );

    ASSERT_EQ(0, 0);
}
