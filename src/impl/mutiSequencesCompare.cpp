//
// Created by zhf on 22-8-23.
//
#include "mutiSequencesCompare.h"


void mutiSequencesCompre(std::map<std::string, std::string> &sequences,std::vector<std::string> &seqNames,int8_t kmer_size,int32_t loopCountend,
                         int32_t mismatchingPenalty, int32_t _open_gap_penalty1,int32_t _extend_gap_penalty1,
                         int32_t _open_gap_penalty2,int32_t _extend_gap_penalty2,std::string& outResultPath ){
    int32_t loopCount=0;
    std::vector<std::string> dna_refs;
    std::vector<std::string> dna_queries;
    std::vector<std::string> Rnames;
    std::vector<std::string> Qnames;
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;

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
        for (std::map<std::string,std::string>::iterator it6=sequences.begin();it6!=sequences.end();++it6) {
            std::cout << it6->first <<"  "<< it6->second << std::endl;
        }
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
    }while((temporary_sequences!=sequences)&&(loopCount<loopCountend));
    std::cout<<"completed compare "<<std::endl;
    std::ofstream ofile;
    ofile.open(outResultPath);
    for( std::string seqName : seqNames ){
        int lineWidth=60;
        writeFasta(ofile, seqName, sequences[seqName],  lineWidth);
    }
    ofile.close();
}