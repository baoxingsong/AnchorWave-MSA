//
// Created by baoxing on 10/10/17.
//

#include "readGffFile.h"

void readGffFile(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const int &minExon) {
    std::string cdsParentRegex = "([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
    readGffFile(filePath, transcriptHashSet, cdsParentRegex, minExon);
}

void get_transcript_to_gene_map_from_gff(const std::string &filePath, std::map <std::string, std::string> &transcript_to_gene_map) {
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.size() < 19 || line[0] == '#') {
            continue;
        }

        std::string id;
        size_t p_id_b = line.find("ID=");
        if (p_id_b != std::string::npos) {
            std::string line_s = line.substr(p_id_b + 3);
            size_t p_id_e = line_s.find(";");
            id = line_s.substr(0, p_id_e);
        } else {
            continue;
        }

        size_t p_p_b = line.find("Parent=");
        if (p_p_b != std::string::npos) {
            std::string line_s = line.substr(p_p_b + 7);
            size_t p_p_e = line_s.find(";");
            std::string parent = line_s.substr(0, p_p_e);
            transcript_to_gene_map[id] = parent;
            continue;
        }

        size_t p_g_b = line.find("geneID=");
        if (p_g_b != std::string::npos) {
            std::string line_s = line.substr(p_g_b + 7);
            size_t p_g_e = line_s.find(";");
            std::string gene_id = line_s.substr(0, p_g_e);
            transcript_to_gene_map[id] = gene_id;
        }
    }

    infile.close();
}

void readGffFile(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
    std::map <std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        size_t pos_cds = line.find("CDS");
        if (pos_cds == std::string::npos || line.size() < 19 || line[0] == '#') {
            continue;
        }

        const char *ss = line.c_str();
        int size = line.size();
        char name[size];
        char type[size];
        int start, end;
        char str6[1];
        char str8[size];

        int ret = sscanf(ss, "%s%*s%s%d%d%*s%s%*s%s", name, type, &start, &end, str6, str8);
        if (ret == 6 && std::string(type) == "CDS") {
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }

            if ((end - start + 1) >= minExon) {
                std::string str8_ = std::string(str8);
                // --1 like :    chr1	NAM	CDS	34722	35318	.	+	0	ID=Zm00001eb000010_P001;Parent=Zm00001eb000010_T001;protein_id=Zm00001eb000010_P001
                size_t p_p_b = str8_.find("Parent=");
                if (p_p_b == std::string::npos)
                    continue;
                str8_ = str8_.substr(p_p_b + 7);
                size_t p_p_e = str8_.find_first_of(";");
                str8_ = str8_.substr(0, p_p_e);

                // --2

                std::string information = str8_;
                if (transcriptHashMap.find(information) != transcriptHashMap.end()) {
                } else {
                    std::string chromosomeName = std::string(name);

                    STRAND strand;
                    if (str6[0] == '-') {
                        strand = NEGATIVE;
                    } else {
                        strand = POSITIVE;
                    }

                    Transcript transcript1(information, chromosomeName, strand);
                    transcriptHashMap[information] = transcript1;
                }

                GenomeBasicFeature cds(start, end);
                transcriptHashMap[information].addCds(cds);
            }
        }
    }

    infile.close();

    for (std::map<std::string, Transcript>::iterator it = transcriptHashMap.begin(); it != transcriptHashMap.end(); ++it) {
        if (transcriptHashSet.find(it->second.getChromeSomeName()) == transcriptHashSet.end()) {
            transcriptHashSet[it->second.getChromeSomeName()] = std::vector<Transcript>();
        }

        it->second.updateInforCds();
        transcriptHashSet[it->second.getChromeSomeName()].push_back(it->second);
    }

    for (std::map < std::string, std::vector < Transcript >> ::iterator it = transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](Transcript a, Transcript b) {
            return a.getPStart() < b.getPStart();
        });
    }
}

void readGffFile_exon(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
    std::map <std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        size_t pos_exon = line.find("exon");
        if (pos_exon == std::string::npos || line.size() < 19 || line[0] == '#') {
            continue;
        }

        const char *ss = line.c_str();
        int size = line.size();
        char name[size];
        char type[size];
        int start, end;
        char str6[1];
        char str8[size];

        int ret = sscanf(ss, "%s%*s%s%d%d%*s%s%*s%s", name, type, &start, &end, str6, str8);
        if (ret == 6 && std::string(type) == "exon") {
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }

            if ((end - start + 1) >= minExon) {
                std::string str8_ = std::string(str8);
                // --1 like :    chr1	NAM	exon	34617	35318	.	+	.	Parent=Zm00001eb000010_T001;Name=Zm00001eb000010_T001.exon.1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001eb000010_T001.exon.1;rank=1
                size_t p_p_b = str8_.find("Parent=");
                if (p_p_b == std::string::npos)
                    continue;
                str8_ = str8_.substr(p_p_b + 7);
                size_t p_p_e = str8_.find_first_of(";");
                str8_ = str8_.substr(0, p_p_e);

                // --2

                std::string information = str8_;
                if (transcriptHashMap.find(information) != transcriptHashMap.end()) {
                } else {
                    std::string chromosomeName = std::string(name);

                    STRAND strand;
                    if (str6[0] == '-') {
                        strand = NEGATIVE;
                    } else {
                        strand = POSITIVE;
                    }

                    Transcript transcript1(information, chromosomeName, strand);
                    transcriptHashMap[information] = transcript1;
                }

                GenomeBasicFeature cds(start, end);
                transcriptHashMap[information].addCds(cds);
            }
        }
    }

    infile.close();

    for (std::map<std::string, Transcript>::iterator it = transcriptHashMap.begin(); it != transcriptHashMap.end(); ++it) {
        if (transcriptHashSet.find(it->second.getChromeSomeName()) == transcriptHashSet.end()) {
            transcriptHashSet[it->second.getChromeSomeName()] = std::vector<Transcript>();
        }

        it->second.updateInforCds();
        transcriptHashSet[it->second.getChromeSomeName()].push_back(it->second);
    }

    for (std::map < std::string, std::vector < Transcript >> ::iterator it = transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](Transcript a, Transcript b) {
            return a.getPStart() < b.getPStart();
        });
    }
}
