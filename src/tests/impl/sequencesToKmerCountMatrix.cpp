//
// Created by song on 8/8/18.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
TEST(sequencesToKmerCountMatrix, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/zhf/testfa.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );
    for( std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > >::iterator it = kmerCountMatrix.begin(); it!=kmerCountMatrix.end(); ++it ){
        std::cout << it->first << std::endl;
        for( std::map<std::string/*kmer*/, int16_t >::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2 ){
            std::cout << it2->first << "\t" << it2->second << std::endl;
        }
    }
    ASSERT_EQ(0, 0);
}

TEST(sequencesToKmerDistanceMatrix, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/zhf/aaa1.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( fastaFilePath, sequences, seqNames);
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );

    float ** distanceMatrix = new float * [sequences.size()];

    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }
    kmerCountMatrixToDistanceMatrix(  kmerCountMatrix, seqNames, distanceMatrix );
    std::cout << seqNames.size()<<std::endl;
    for(size_t i = 0; i < seqNames.size(); ++i){
        std::cout << seqNames[i] ;
        for(size_t j = 0; j < seqNames.size(); ++j){
            std::cout << "\t" << distanceMatrix[i][j];
        }
        std::cout<< std::endl;
    }

    for(size_t i = 0; i < sequences.size(); ++i){
        delete[] distanceMatrix[i];
    }
    delete[] distanceMatrix;

    ASSERT_EQ(0, 0);
}


TEST(sequencesToKmeruUpgma, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/zhf/aaa1.fa";
//    std::string fastaFilePath = "/home/zhf/normalfa.fa";
    std::string transquences;
    std::string temp;
    int8_t kmer_size = 12;
    int8_t differ_size=1;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( fastaFilePath, sequences, seqNames);
    //delete
    for(std::map<std::string, std::string>::iterator itsv=sequences.begin();itsv!=sequences.end();++itsv){
        std::cout<<"             init sequences        "<<   itsv->first<<std::endl;
    }
    for(int i=0; i<seqNames.size(); ++i){
        std::string seqName = seqNames[i];
        seqName=songStrReplaceAll( seqName, ",", "_"  );
        seqName=songStrReplaceAll( seqName, "\\(", "_"  );
        seqName=songStrReplaceAll( seqName, "\\)", "_"  );
        std::string seq = sequences[seqNames[i]];
        sequences.erase(seqNames[i]);
        sequences.insert(make_pair(seqName,seq));
        std::cout << seqNames[i] << "\t" << seqName << std::endl;
        seqNames[i] = seqName;
    }
    //delete
    for (std::map<std::string,std::string>::iterator it6=sequences.begin();it6!=sequences.end();++it6){
        std::cout <<"traverse    "<< it6->first <<"  "<< it6->second << std::endl;
    }
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;
//    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > differCountMatrix;
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );
    std::vector<std::string> dna_refs;
    std::vector<std::string> dna_queries;
    std::vector<std::string> Rnames;
    std::vector<std::string> Qnames;

    float ** distanceMatrix = new float * [sequences.size()];
    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }

    kmerCountMatrixToDistanceMatrix(  kmerCountMatrix, seqNames, distanceMatrix );
    //delete
    for(size_t i = 0; i < seqNames.size(); ++i){
        std::cout << seqNames[i] ;
        for(size_t j = 0; j < seqNames.size(); ++j){
            std::cout << "\t" << distanceMatrix[i][j];
        }
        std::cout<< std::endl;
    }
    using namespace std;
    // let's start with empty DynMatrix:
    ClusterNode* head = NULL;
    ClusterNode* tail = NULL;

    for (string seqName : seqNames){
        addCluster(head, tail, seqName);
    }
    ClusterNode *node = head;
    for (int i=0; i<seqNames.size(); i++) {
        DistanceNode *newDistance = node->row;
        for (int j=0; j<seqNames.size(); j++) {
            double d = distanceMatrix[i][j];
            newDistance->distance = d;
            newDistance = newDistance->nextInRow;
        }
        node  = node->next;
    }
    // Run UPGMA until there is only one node left
    while (head->next) {
        // Create temporary Cluster and Distance Nodes
        ClusterNode *tempHead = head;

        ClusterNode *clusterOne = NULL;
        ClusterNode *clusterTwo = NULL;
        DistanceNode *distanceOne;
        DistanceNode *distanceTwo;
        // Create vectors to hold distances of DistanceNodes
        std::vector<double> clusterOneDistances;
        std::vector<double> clusterTwoDistances;
        std::vector<double> resultantDistances;
        smatch result;
        regex pattern("[^(,)]+");
        int32_t mismatchingPenalty = -3;
        int32_t _open_gap_penalty1 = -4;
        int32_t _extend_gap_penalty1 = -2;
        int32_t _open_gap_penalty2 = -80;
        int32_t _extend_gap_penalty2 = -1;
        findMinimum(tempHead, clusterOne, clusterTwo);
        std::cout<<"167"<<std::endl;
        for (distanceOne = clusterOne->column; distanceOne != NULL; distanceOne = distanceOne->nextInColumn){
            clusterOneDistances.push_back(distanceOne->distance);
        }
        for (distanceTwo = clusterTwo->column; distanceTwo != NULL; distanceTwo = distanceTwo->nextInColumn){
            clusterTwoDistances.push_back(distanceTwo->distance);
        }
        resultantDistances = useFormula(clusterOne, clusterTwo, clusterOneDistances, clusterTwoDistances);
        std::cout<<"175"<<std::endl;
        string::const_iterator iterStartR = clusterOne->name.begin();
        string::const_iterator iterEndR = clusterOne->name.end();
        string::const_iterator iterStartQ = clusterTwo->name.begin();
        string::const_iterator iterEndQ = clusterTwo->name.end();
        std::cout<<"180"<<std::endl;
        while (regex_search(iterStartR, iterEndR, result, pattern)){
            Rnames.push_back(result[0]);
            iterStartR = result[0].second;
        }
        while (regex_search(iterStartQ, iterEndQ, result, pattern)){
            Qnames.push_back(result[0]);
            iterStartQ = result[0].second;
        }
        std::cout<<"189"<<std::endl;
        combineCluster(head, tail, clusterOne, clusterTwo, resultantDistances);
        std::cout<<"191"<<std::endl;
        for(int j=0;j<Rnames.size();++j){
            for (std::map<std::string,std::string>::iterator it1=sequences.begin();it1!=sequences.end();++it1){
                if(sequences.find(Rnames[j])==it1){
                    dna_refs.push_back(it1->second);
                }
            }
        }
        for(int k=0;k<Qnames.size();++k){
            for (std::map<std::string,std::string>::iterator it2=sequences.begin();it2!=sequences.end();++it2){
                if(sequences.find(Qnames[k])==it2){
                    dna_queries.push_back(it2->second);
                }
            }
        }
        std::cout<<"273"<<std::endl;
        for (int i=0;i<dna_refs.size();++i){
            std::cout<<"cin "<<"dna_refs "<<dna_refs[i]<<std::endl;
        }
        std::cout<<std::endl;
        for (int i = 0; i < dna_queries.size(); ++i) {
            std::cout << "cin  "<<"dna_queries "<<dna_queries[i] << std::endl;
        }
        std::vector<std::stack<char>> SQs(dna_queries.size());
        std::vector<std::stack<char>> SRs(dna_refs.size());
        int32_t score = needleAlignment(dna_refs, dna_queries, SQs, SRs, mismatchingPenalty,_open_gap_penalty1,_extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2);

        std::vector<std::string> _alignment_qs(dna_queries.size());
        std::vector<std::string> _alignment_ds(dna_refs.size());

        while (!SQs[0].empty()) {
            for( int k=0; k<dna_refs.size(); ++k ){
                _alignment_ds[k] += SRs[k].top();
                SRs[k].pop();
            }
            for( int l=0; l<dna_queries.size(); ++l ){
                _alignment_qs[l] += SQs[l].top();
                SQs[l].pop();
            }
        }
        //delete
        for (int i=0;i<_alignment_ds.size();++i){
            std::cout <<"Compare complete "<<"dna_refs      "<<_alignment_ds[i]<< std::endl;
        }
        for (int i = 0; i < _alignment_qs.size(); ++i) {
            std::cout <<"Compare complete "<<"dna_queries   "<<_alignment_qs[i] << std::endl;
        }
        for(int m=0 ,z=0;m<Rnames.size();++m,++z){
            for (std::map<std::string, std::string>::iterator it3 = sequences.begin();it3 != sequences.end(); ++it3) {
                if (it3->first == Rnames[m]){
                    sequences[Rnames[m]] = _alignment_ds[z];
                }
            }
        }

        for(int n=0,z=0;n<Qnames.size();++n,++z){
            for (std::map<std::string, std::string>::iterator it4 = sequences.begin();it4 != sequences.end(); ++it4) {
                if (it4->first == Qnames[n]) {
                    sequences[Qnames[n]] = _alignment_qs[z];
                }
            }
        }
        if ((dna_refs.size() + dna_queries.size() ) == sequences.size()) {
              for (std::map<std::string,std::string>::iterator it5=sequences.begin();it5!=sequences.end();++it5) {
                std::cout << it5->first <<"  "<< it5->second << std::endl;
              }
                std::cout << score << std::endl;
            break;
        }
        for (std::map<std::string,std::string>::iterator it6=sequences.begin();it6!=sequences.end();++it6) {
            std::cout << it6->first <<"  "<< it6->second << std::endl;
        }
        std::cout << score << std::endl;
        dna_refs.clear();
        dna_queries.clear();
        Rnames.clear();
        Qnames.clear();
        std::cout << "end============================================================" << std::endl;
    }
    // Print out name of last node
    if (head) {
        cout  << head->name <<endl;
    }
    // BONUS (optional): print the tree in a nice way
    std::cout<<"273"<<std::endl;
    std::ofstream ofile;
    ofile.open("/home/zhf/completed.fa");
    //for( std::map<std::string, std::string>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
    for( std::string seqName : seqNames ){
        int lineWidth=60;
        writeFasta(ofile, seqName, sequences[seqName],  lineWidth);
    }

    ofile.close();
    ASSERT_EQ(0, 0);
}





TEST(sequencesToKmeruUpgma, c2){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/zhf/aaa1.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    std::vector<std::string> dna_refs;
    std::vector<std::string> dna_queries;
    std::vector<std::string> Rnames;
    std::vector<std::string> Qnames;
    int32_t mismatchingPenalty = -3;
    int32_t _open_gap_penalty1 = -4;
    int32_t _extend_gap_penalty1 = -2;
    int32_t _open_gap_penalty2 = -80;
    int32_t _extend_gap_penalty2 = -1;
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;

    readFastaFile( fastaFilePath, sequences, seqNames);

    for(int i=0; i<seqNames.size(); ++i){
        std::string seqName = seqNames[i];
        seqName=songStrReplaceAll( seqName, ",", "_"  );
        seqName=songStrReplaceAll( seqName, "\\(", "_"  );
        seqName=songStrReplaceAll( seqName, "\\)", "_"  );
        std::string seq = sequences[seqNames[i]];
        sequences.erase(seqNames[i]);
        sequences.insert(make_pair(seqName,seq));
        std::cout << seqNames[i] << "\t" << seqName << std::endl;
        seqNames[i] = seqName;
    }

    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );

    float ** distanceMatrix = new float * [sequences.size()];
    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }

    kmerCountMatrixToDistanceMatrix(  kmerCountMatrix, seqNames, distanceMatrix );

    for(size_t i = 0; i < seqNames.size(); ++i){
        std::cout << seqNames[i] ;
        for(size_t j = 0; j < seqNames.size(); ++j){
            std::cout << "\t" << distanceMatrix[i][j];
        }
        std::cout<< std::endl;
    }

    using namespace std;
    // let's start with empty DynMatrix:
    ClusterNode* head = NULL;
    ClusterNode* tail = NULL;

    for (string seqName : seqNames){
        addCluster(head, tail, seqName);
    }
    ClusterNode *node = head;
    for (int i=0; i<seqNames.size(); i++) {
        DistanceNode *newDistance = node->row;
        for (int j=0; j<seqNames.size(); j++) {
            double d = distanceMatrix[i][j];
            newDistance->distance = d;
            newDistance = newDistance->nextInRow;
        }
        node  = node->next;
    }
//     Run UPGMA until there is only one node left
    while (head->next) {

        printRowByRow(head);
        cout << "--------------------------------------------------------------" << endl;

        // Create temporary Cluster and Distance Nodes
        UPGMA(head, tail,sequences,Rnames,Qnames,dna_refs,dna_queries);
        std::vector<std::stack<char>> SQs(dna_queries.size());
        std::vector<std::stack<char>> SRs(dna_refs.size());

        int32_t score = needleAlignment(dna_refs, dna_queries, SQs, SRs, mismatchingPenalty,_open_gap_penalty1,_extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2);

        std::vector<std::string> _alignment_qs(dna_queries.size());
        std::vector<std::string> _alignment_ds(dna_refs.size());

        updataSequences(sequences,_alignment_ds, _alignment_qs,SRs,SQs,Rnames,Qnames,dna_refs,dna_queries);

        if ((dna_refs.size() + dna_queries.size() ) == sequences.size()) {
            std::cout << score << std::endl;
            break;
        }
        std::cout << score << std::endl;
        dna_refs.clear();
        dna_queries.clear();
        Rnames.clear();
        Qnames.clear();
        std::cout << "end============================================================" << std::endl;
    }
    if (head) {
        cout  << head->name <<endl;
    }
    std::cout<<"273"<<std::endl;

    std::ofstream ofile;
    ofile.open("/home/zhf/completed.fa");
    for( std::string seqName : seqNames ){
        int lineWidth=60;
        writeFasta(ofile, seqName, sequences[seqName],  lineWidth);
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}














TEST(sequencesToKmeruUpgma, c3){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/zhf/aaa1.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    std::vector<std::string> dna_refs;
    std::vector<std::string> dna_queries;
    std::vector<std::string> Rnames;
    std::vector<std::string> Qnames;
    int32_t mismatchingPenalty = -3;
    int32_t _open_gap_penalty1 = -4;
    int32_t _extend_gap_penalty1 = -2;
    int32_t _open_gap_penalty2 = -80;
    int32_t _extend_gap_penalty2 = -1;

    readFastaFile( fastaFilePath, sequences, seqNames);

    float ** distanceMatrix = new float * [sequences.size()];
    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }
    sequencesToDifferDistanceMatrix( sequences, seqNames,distanceMatrix);

    for(size_t i = 0; i < seqNames.size(); ++i){
        std::cout << seqNames[i] ;
        for(size_t j = 0; j < seqNames.size(); ++j){
            std::cout << "\t" << distanceMatrix[i][j];
        }
        std::cout<< std::endl;
    }

    using namespace std;
    // let's start with empty DynMatrix:
    ClusterNode* head = NULL;
    ClusterNode* tail = NULL;

    for (string seqName : seqNames){
        std::cout<<"name    "<<seqName<<std::endl;
        addCluster(head, tail, seqName);
    }
    ClusterNode *node = head;
    for (int i=0; i<seqNames.size(); i++) {
        DistanceNode *newDistance = node->row;
        for (int j=0; j<seqNames.size(); j++) {
            double d = distanceMatrix[i][j];
            newDistance->distance = d;
            newDistance = newDistance->nextInRow;
        }
        node  = node->next;
    }
//     Run UPGMA until there is only one node left
    while (head->next) {

        printRowByRow(head);
        cout << "--------------------------------------------------------------" << endl;

        // Create temporary Cluster and Distance Nodes
        UPGMA(head, tail,sequences,Rnames,Qnames,dna_refs,dna_queries);
        std::vector<std::stack<char>> SQs(dna_queries.size());
        std::vector<std::stack<char>> SRs(dna_refs.size());
        int32_t score = needleAlignment(dna_refs, dna_queries, SQs, SRs, mismatchingPenalty,_open_gap_penalty1,_extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2);
        std::vector<std::string> _alignment_qs(dna_queries.size());
        std::vector<std::string> _alignment_ds(dna_refs.size());
        updataSequences(sequences,_alignment_ds, _alignment_qs,SRs,SQs,Rnames,Qnames,dna_refs,dna_queries);
        if ((dna_refs.size() + dna_queries.size() ) == sequences.size()) {
            std::cout << score << std::endl;
            break;
        }
        std::cout << score << std::endl;
        dna_refs.clear();
        dna_queries.clear();
        Rnames.clear();
        Qnames.clear();
        std::cout << "end============================================================" << std::endl;
    }
    if (head) {
        cout  << head->name <<endl;
    }
    std::cout<<"273"<<std::endl;
    std::ofstream ofile;
    ofile.open("/home/zhf/completed.fa");
    for( std::string seqName : seqNames ){
        int lineWidth=60;
        writeFasta(ofile, seqName, sequences[seqName],  lineWidth);
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}







TEST(sequencesToKmeruUpgma, c99){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/zhf/aaa2.fa";
    int8_t kmer_size = 12;
    int loopCount=0;
    int loopCountend=50;
    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    std::vector<std::string> dna_refs;
    std::vector<std::string> dna_queries;
    std::vector<std::string> Rnames;
    std::vector<std::string> Qnames;
    int32_t mismatchingPenalty = -3;
    int32_t _open_gap_penalty1 = -4;
    int32_t _extend_gap_penalty1 = -2;
    int32_t _open_gap_penalty2 = -80;
    int32_t _extend_gap_penalty2 = -1;
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;

    readFastaFile( fastaFilePath, sequences, seqNames);

    for(int i=0; i<seqNames.size(); ++i){
        std::string seqName = seqNames[i];
        seqName=songStrReplaceAll( seqName, ",", "_"  );
        seqName=songStrReplaceAll( seqName, "\\(", "_"  );
        seqName=songStrReplaceAll( seqName, "\\)", "_"  );
        std::string seq = sequences[seqNames[i]];
        sequences.erase(seqNames[i]);
        sequences.insert(make_pair(seqName,seq));
        std::cout << seqNames[i] << "\t" << seqName << std::endl;
        seqNames[i] = seqName;
    }
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );
    float ** distanceMatrix = new float * [sequences.size()];
    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }
    kmerCountMatrixToDistanceMatrix(  kmerCountMatrix, seqNames, distanceMatrix );
    using namespace std;
    // let's start with empty DynMatrix:
    ClusterNode* head = NULL;
    ClusterNode* tail = NULL;

    for (string seqName : seqNames){
        addCluster(head, tail, seqName);
    }
    ClusterNode *node = head;
    for (int i=0; i<seqNames.size(); i++) {
        DistanceNode *newDistance = node->row;
        for (int j=0; j<seqNames.size(); j++) {
            double d = distanceMatrix[i][j];
            newDistance->distance = d;
            newDistance = newDistance->nextInRow;
        }
        node  = node->next;
    }
//     Run UPGMA until there is only one node left
    while (head->next) {
        // Create temporary Cluster and Distance Nodes
        UPGMA(head, tail,sequences,Rnames,Qnames,dna_refs,dna_queries);
        std::vector<std::stack<char>> SQs(dna_queries.size());
        std::vector<std::stack<char>> SRs(dna_refs.size());
        int32_t score = needleAlignment(dna_refs, dna_queries, SQs, SRs, mismatchingPenalty,_open_gap_penalty1,_extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2);

        std::vector<std::string> _alignment_qs(dna_queries.size());
        std::vector<std::string> _alignment_ds(dna_refs.size());
        updataSequences(sequences,_alignment_ds, _alignment_qs,SRs,SQs,Rnames,Qnames,dna_refs,dna_queries);

        if ((dna_refs.size() + dna_queries.size() ) == sequences.size()) {
            std::cout << score << std::endl;
            break;
        }
        std::cout << score << std::endl;
        dna_refs.clear();
        dna_queries.clear();
        Rnames.clear();
        Qnames.clear();
    }
    if (head) {
        cout  << head->name << endl;
        dna_refs.clear();
        dna_queries.clear();
        Rnames.clear();
        Qnames.clear();
    }
    std::map<std::string,std::string> temporary_sequences;
    do {
        loopCount++;
        float ** distanceMatrix_l = new float * [sequences.size()];
        for(size_t i = 0; i < sequences.size(); ++i){
            distanceMatrix_l[i] = new float [sequences.size()];
        }

        sequencesToDifferDistanceMatrix( sequences, seqNames,distanceMatrix_l);

        for(size_t i = 0; i < seqNames.size(); ++i){
            std::cout << seqNames[i] ;
            for(size_t j = 0; j < seqNames.size(); ++j){
                std::cout << "\t" << distanceMatrix_l[i][j];
            }
            std::cout<< std::endl;
        }
        temporary_sequences=sequences;
        using namespace std;
        // let's start with empty DynMatrix:
        ClusterNode* head_l = NULL;
        ClusterNode* tail_l = NULL;

        for (string seqName : seqNames){
            addCluster(head_l, tail_l, seqName);
        }
        ClusterNode *node_l = head_l;
        for (int i=0; i<seqNames.size(); i++) {
            DistanceNode *newDistance_l = node_l->row;
            for (int j=0; j<seqNames.size(); j++) {
                double d = distanceMatrix_l[i][j];
                newDistance_l->distance = d;
                newDistance_l = newDistance_l->nextInRow;
            }
            node_l  = node_l->next;
        }
//     Run UPGMA until there is only one node left
        while (head_l->next) {
//            printRowByRow(head_l);
//            printColumnByColumn(head_l);
//            cout << "--------------------------------------------------------------" << endl;
            // Create temporary Cluster and Distance Nodes
            UPGMA(head_l, tail_l,sequences,Rnames,Qnames,dna_refs,dna_queries);
            std::vector<std::stack<char>> SQs(dna_queries.size());
            std::vector<std::stack<char>> SRs(dna_refs.size());
            int32_t score = needleAlignment(dna_refs, dna_queries, SQs, SRs, mismatchingPenalty,_open_gap_penalty1,_extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2);
            std::vector<std::string> _alignment_qs(dna_queries.size());
            std::vector<std::string> _alignment_ds(dna_refs.size());
            updataSequences(sequences,_alignment_ds, _alignment_qs,SRs,SQs,Rnames,Qnames,dna_refs,dna_queries);
            if ((dna_refs.size() + dna_queries.size() ) == sequences.size()) {
                std::cout << score << std::endl;
                break;
            }
            std::cout << score << std::endl;
            dna_refs.clear();
            dna_queries.clear();
            Rnames.clear();
            Qnames.clear();
        }
        std::cout<<"641"<<std::endl;
    }while((temporary_sequences!=sequences)&&(loopCount<loopCountend));

    std::cout<<"completed compare "<<std::endl;
    std::ofstream ofile;
    ofile.open("/home/zhf/completed.fa");
    for( std::string seqName : seqNames ){
        int lineWidth=60;
        writeFasta(ofile, seqName, sequences[seqName],  lineWidth);
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}

//Same code can run  fffff

TEST(sequencesToKmeruUpgma, c199) { // just to make sure that every line has been analysed

    std::map<std::string,std::string> temporary_sequences;
    std::map<std::string, std::string> sequences_l;
    do {
        std::string fastaFilePath = "/home/zhf/completed.fa";
        std::vector<std::string> seqNames_l;
        std::vector<std::string> dna_refs_l;
        std::vector<std::string> dna_queries_l;
        std::vector<std::string> Rnames;
        std::vector<std::string> Qnames;
        int32_t mismatchingPenalty = -3;
        int32_t _open_gap_penalty1 = -4;
        int32_t _extend_gap_penalty1 = -2;
        int32_t _open_gap_penalty2 = -80;
        int32_t _extend_gap_penalty2 = -1;
        readFastaFile(fastaFilePath, sequences_l, seqNames_l);
        float **distanceMatrix_l = new float *[sequences_l.size()];
        for (size_t i = 0; i < sequences_l.size(); ++i) {
            distanceMatrix_l[i] = new float[sequences_l.size()];
        }
        sequencesToDifferDistanceMatrix(sequences_l, seqNames_l, distanceMatrix_l);

        for (size_t i = 0; i < seqNames_l.size(); ++i) {
            std::cout << seqNames_l[i];
            for (size_t j = 0; j < seqNames_l.size(); ++j) {
                std::cout << "\t" << distanceMatrix_l[i][j];
            }
            std::cout << std::endl;
        }

        using namespace std;
        // let's start with empty DynMatrix:
        ClusterNode *head_l = NULL;
        ClusterNode *tail_l = NULL;

        for (string seqName: seqNames_l) {
//            std::cout << "name    " << seqName << std::endl;
            addCluster(head_l, tail_l, seqName);
        }
        ClusterNode *node_l = head_l;
        for (int i = 0; i < seqNames_l.size(); i++) {
            DistanceNode *newDistance = node_l->row;
            for (int j = 0; j < seqNames_l.size(); j++) {
                double d = distanceMatrix_l[i][j];
                newDistance->distance = d;
                newDistance = newDistance->nextInRow;
            }
            node_l = node_l->next;
        }
//     Run UPGMA until there is only one node left
        while (head_l->next) {
            // Create temporary Cluster and Distance Nodes
            UPGMA(head_l, tail_l, sequences_l, Rnames, Qnames, dna_refs_l, dna_queries_l);
            std::vector<std::stack<char>> SQs(dna_queries_l.size());
            std::vector<std::stack<char>> SRs(dna_refs_l.size());
            int32_t score = needleAlignment(dna_refs_l, dna_queries_l, SQs, SRs, mismatchingPenalty, _open_gap_penalty1,
                                            _extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2);
            std::vector<std::string> _alignment_qs(dna_queries_l.size());
            std::vector<std::string> _alignment_ds(dna_refs_l.size());
            updataSequences(sequences_l, _alignment_ds, _alignment_qs, SRs, SQs, Rnames, Qnames, dna_refs_l,
                            dna_queries_l);
            if ((dna_refs_l.size() + dna_queries_l.size()) == sequences_l.size()) {
                std::cout << score << std::endl;
                break;
            }
            std::cout << score << std::endl;
            dna_refs_l.clear();
            dna_queries_l.clear();
            Rnames.clear();
            Qnames.clear();
            std::cout << "end============================================================" << std::endl;
            if (head_l) {
                cout  << head_l->name <<endl;
            }
            std::cout<<"273"<<std::endl;

        }
    }while(temporary_sequences==sequences_l);
    ASSERT_EQ(0, 0);
}



