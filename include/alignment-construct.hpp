#ifndef construct_hpp
#define construct_hpp
#include <iostream>
#include <vector>
#include "zlib.h"
#include "map"
#include <cmath>
#include <algorithm>
#include "kseq.h"
#include "get_time.hpp"
#include "split.hpp"
// #include "src/RBTree.hpp"

using namespace std;

// kseq.h ���ļ�
KSEQ_INIT(gzFile, gzread)

struct minimizerSrt
{
    long long int _end;
    string _chromosome;
};


namespace construct
{
    long long int hash_fun(string sequence);
    pair<string, string> find_node_up_down_seq(const string genotype, 
                                               const long long int seqLen,
                                               const map<long long int, vector<map<string, string>>>::iterator & iter, 
                                               const map<long long int, vector<map<string, string>>> & startMapTmp);

    string sequence_select(string sequence);

    //////////////////////////////////////////////////////////////////////////

    /**
     * @author zezhen du
     * @date 2022/12/13
     * @version v1.0
	 * @brief build refgenome indexs
     * 
     * @param inputFasta     ��Ҫ����������fasta�ļ�
     * 
     * @return fastaMap      map<chr, sequence>
	**/
    map<string, string> build_fasta_index(string inputFasta)
    {
        map<string, string> fastaMap;

        // open fasta file
        gzFile gzfp = gzopen(inputFasta.c_str(), "rb");

        // ���ļ�
        if(!gzfp)
        {
            cerr << "[" << getTime() << "] " 
                << inputFasta 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }
        else
        {
            kseq_t *ks;
            ks = kseq_init(gzfp);

            while( kseq_read(ks) >= 0 )
            {
                // ks->name.s ��¼��������
                // ks->seq.s ��¼��������
                string chromosome = ks->name.s;
                string sequence = ks->seq.s;

                // ����fasta����
                fastaMap[chromosome] = sequence;
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        return fastaMap;
    }

    // ��vcf�ļ�
    void vcf_construct(string inputVcf, map<string, string> & fastaMap, 
                        map<string, map<long long int, vector<map<string, string>>>> & vcfGraph)
    {
        // open base fastafile
        gzFile gzfp = gzopen(inputVcf.c_str(), "rb");

        // ��¼��һ���ڵ��λ��
        long int tmpRefEnd;
        string tmpChromosome;

        // ���ļ�
        if(!gzfp)
        {
            cerr << "[" << getTime() << "] " 
                << inputVcf 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }
        else
        {
            string information;

            while (!gzeof(gzfp))
            {
                char tmpChar = gzgetc(gzfp);

                if (tmpChar != '\n') // ���û�ж������һ���ַ�����һֱ��ӡ�
                {
                    information += tmpChar;
                }
                else // ������й�ͼ
                {
                    // ���ж��ǲ���ע���У��ǵĻ�����
                    if (information.find("#") != string::npos)
                    {
                        // ����ַ���
                        information.clear();
                        string().swap(information);
                    }
                    else // ��ͼ
                    {
                        vector<string> informationVec = split(information, "\t");

                        if (informationVec.size() < 9) // �ȼ���ļ��Բ��ԣ���������������
                        {
                            cerr << "[" << getTime() << "] " 
                                << inputVcf 
                                << ": File format error." 
                                << " [" << information << "]"
                                << endl;
                            exit(1);
                        }

                        // ����ĸ�����Ϣ
                        string chromosome = informationVec[0];
                        long long int refStart = stol(informationVec[1]);
                        string refSeq = informationVec[3];
                        long long int refEnd = refStart + refSeq.size() - 1;
                        string qrySeq = informationVec[4];

                        // ���vcf��qry_seq�ͻ������Ӧ����Ӧ
                        if (fastaMap[chromosome].substr(refStart-1, refSeq.length()) != refSeq)
                        {
                            cerr << "[" << getTime() << "] " 
                                << "Sequence difference between refgenome and vcf: "
                                << fastaMap[chromosome].substr(refStart-1, refSeq.length()) 
                                << " --> " 
                                << refSeq 
                                << endl;
                            exit(1);
                        }

                        // ͼ�еĽڵ���
                        long long int nodeNum = vcfGraph.size();

                        // ����vcfͼ��
                        map<string, string> tmpMap; // ���кͻ�������Ϣ
                        
                        // ��� genotype=0 �Ľڵ�
                        // ���Ⱦɫ���һ��vcfλ�ò���1����ѻ�����ǰ�ߵ����й�����һ��node
                        string genotype = "0";
                        if (chromosome != tmpChromosome && refStart > 0)
                        {
                            // �����Ⱦɫ���ǰ����
                            string preRefSeq;
                            long long int preRefStart;
                            long long int preRefEnd;
                            long long int preRefLen;
                            
                            // ��Ӿ�Ⱦɫ��ĺ�벿��
                            if (nodeNum > 1 && tmpRefEnd < fastaMap[tmpChromosome].length()) // �ж����һ�������ڲ���Ⱦɫ�����
                            {
                                preRefStart = tmpRefEnd + 1;
                                preRefEnd = fastaMap[tmpChromosome].length();
                                preRefLen = preRefEnd - preRefStart + 1;

                                preRefSeq = fastaMap[tmpChromosome].substr(tmpRefEnd, preRefLen);
                                
                                tmpMap["sequence"] = preRefSeq;
                                tmpMap["genotype"] = genotype;
                                vcfGraph[tmpChromosome][preRefStart].push_back(tmpMap);
                            }

                            // �����Ⱦɫ���ǰ����
                            preRefStart = 1;
                            preRefEnd = refStart-1;
                            preRefLen = preRefEnd - preRefStart + 1;
                            preRefSeq = fastaMap[chromosome].substr(0, preRefLen);
                            
                            tmpMap["sequence"] = preRefSeq;
                            tmpMap["genotype"] = genotype;
                            vcfGraph[chromosome][preRefStart].push_back(tmpMap);

                            tmpRefEnd = refEnd;
                            tmpChromosome = chromosome;
                        }
                        else // �����vcf�м�����й���node
                        {
                            string preRefSeq;
                            long long int preRefStart;
                            long long int preRefEnd;
                            long long int preRefLen;
                            
                            preRefStart = tmpRefEnd + 1;
                            preRefEnd = refStart - 1;
                            preRefLen = preRefEnd - preRefStart + 1;

                            // �ж��ǲ������ڣ����ڵĻ�preRefLenС�ڵ���0������
                            if (preRefLen > 0)
                            {
                                preRefSeq = fastaMap[chromosome].substr(tmpRefEnd, preRefLen);

                                tmpMap["sequence"] = preRefSeq;
                                tmpMap["genotype"] = genotype;
                                vcfGraph[chromosome][preRefStart].push_back(tmpMap);
                                
                                tmpRefEnd = refEnd;
                                tmpChromosome = chromosome;
                            }
                        }
                        
                        // ���vcf�ڵ�
                        // �����ref�Ľڵ�
                        tmpMap["sequence"] = refSeq;
                        tmpMap["genotype"] = genotype;
                        vcfGraph[chromosome][refStart].push_back(tmpMap);

                        // qry���ܺ��ж����λ���򣬰������
                        vector<string> qrySeqVec = split(qrySeq, ",");
                        for (auto it : qrySeqVec)
                        {
                            genotype = to_string(stoi(genotype) + 1);
                            tmpMap["sequence"] = it;
                            tmpMap["genotype"] = genotype;
                            vcfGraph[chromosome][refStart].push_back(tmpMap);
                        }
                    }

                    // ����ַ���
                    information.clear();
                    string().swap(information);
                }
            }

            // ���һ��Ⱦɫ��ĺ�벿��
            if (tmpRefEnd < fastaMap[tmpChromosome].length()) // �ж����һ�������ڲ���Ⱦɫ�����
            {
                long long int preRefStart = tmpRefEnd + 1;
                long long int preRefEnd = fastaMap[tmpChromosome].length();
                long long int preRefLen = preRefEnd - preRefStart + 1;

                string preRefSeq = fastaMap[tmpChromosome].substr(tmpRefEnd, preRefLen);

                map<string, string> tmpMap; // �洢��ʼ����ֹ�����е���Ϣ
                
                tmpMap["sequence"] = preRefSeq;
                tmpMap["genotype"] = "0";
                vcfGraph[tmpChromosome][preRefStart].push_back(tmpMap);
            }

            // �ر��ļ����ͷ��ڴ�
            gzclose(gzfp);
        }
    }

    // ΪvcfGraph��������
    void graphIndex(const map<string, map<long long int, vector<map<string, string>>>> & vcfGraph,
               int32_t kmerLen, 
               int minimizerK, 
               int minimizerW,
               int threads_num)
    {
        // ��ӦȾɫ���map
        map<long long int, vector<map<string, string>>> startMapTmp;
        
        for (auto it1 : vcfGraph)
        {
            startMapTmp = it1.second;

            for (map<long long int, vector<map<string, string>>>::iterator iter = startMapTmp.begin();
                 iter != startMapTmp.end(); iter++)
            {
                // �ڵ����ʼλ��
                long long int nodeStart = iter->first;

                // ��ʱ�����������ڼ���һ��node�Ľ�β����
                auto iterTmp = iter;

                for (auto it2 : iter->second)
                {
                    string sequence = it2["sequence"];
                    string genotype = it2["genotype"];

                    // ��������һ��read������
                    pair<string, string> outSeq = find_node_up_down_seq(genotype, minimizerK-1, iter, startMapTmp);
                    sequence = outSeq.first + sequence + outSeq.second;

                    cerr << "chromosome: "
                        << it1.first 
                        << "\tstart:" 
                        << nodeStart 
                        << "\tend:" 
                        << sequence.size()+nodeStart-1 
                        << "\tGT:" 
                        << genotype
                        << "\tfirst:" 
                        << outSeq.first
                        << "\tsecond:"
                        << outSeq.second
                        << endl;
                }
            }
        }
    }

    // �ҽڵ�ǰ�� n bp�����У����ڹ���minimizer
    // find_node_up_down_seq(this node seq, this node GT, this node iterator)
    pair<string, string> find_node_up_down_seq(const string genotype, 
                                               const long long int seqLen,
                                               const map<long long int, vector<map<string, string>>>::iterator & iter, 
                                               const map<long long int, vector<map<string, string>>> & startMapTmp)
    {
        // ��¼�ڵ�����
        string upSeq = "";
        string downSeq = "";

        // ��ʱ������
        auto iterTmp = iter;

        // ��������һ���ڵ��ĩβ����
        while (upSeq.size() < seqLen && iterTmp != startMapTmp.begin())
        {
            // ������ǰ��
            iterTmp--;

            string seqTmp = "";
            string refSeqTmp = "";

            // �Ըýڵ���������н��б���
            for (auto it : iterTmp->second)
            {
                // ����͸ýڵ�Ļ�����һ�£����sequenceTmp��ֵ
                if (it["genotype"] == genotype)
                {
                    seqTmp = it["sequence"];
                }
                if (it["genotype"] == "0")
                {
                    refSeqTmp = it["sequence"];
                }
            }

            // �ж���û�к͸ýڵ�Ļ�����һ�µ�map��û�еĻ���sequenceTmp��ֵref������
            if (seqTmp.size() == 0)
            {
                seqTmp = refSeqTmp;
            }
            
            if (seqTmp.size() >= seqLen-upSeq.size())
            {
                upSeq = seqTmp.substr(seqTmp.size()-(seqLen-upSeq.size()), seqLen-upSeq.size()) + upSeq;
            }
            else
            {
                upSeq = seqTmp + upSeq;
            }
        }

        // ���õ�����
        iterTmp = iter;

        // ��������һ���ڵ�Ŀ�ʼ���У����������Ƽӵ��ж������
        while ((downSeq.size() < seqLen) && (++iterTmp != startMapTmp.end()))
        {
            string seqTmp = "";
            string refSeqTmp = "";

            // �Ըýڵ���������н��б���
            for (auto it : iterTmp->second)
            {
                // ����͸ýڵ�Ļ�����һ�£����sequenceTmp��ֵ
                if (it["genotype"] == genotype)
                {
                    seqTmp = it["sequence"];
                }
                if (it["genotype"] == "0")
                {
                    refSeqTmp = it["sequence"];
                }
            }

            // �ж���û�к͸ýڵ�Ļ�����һ�µ�map��û�еĻ���sequenceTmp��ֵref������
            if (seqTmp.size() == 0)
            {
                seqTmp = refSeqTmp;
            }
            
            if (seqTmp.size() >= seqLen-downSeq.size())
            {
                downSeq = downSeq + seqTmp.substr(0, seqLen-downSeq.size());
            }
            else
            {
                downSeq = downSeq + seqTmp;
            }
        }

        return make_pair(upSeq, downSeq);
    }

    // build kmer index
    map<long long int,string> build_kmer_index(long long int seqStart, string sequence, int32_t kmerLen)
    {
        // ת��д
        transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

        // �洢�����map<start, sequence>
        map<long long int,string> kmerMap;

        int seqLen = sequence.size();
        int kmerNumber = seqLen - kmerLen + 1;

        int seqStartTmp = 0;

        for (int i = 0; i < kmerNumber; i++)
        {
            // ����hash��
            string kmerSeq = sequence.substr(seqStartTmp, kmerLen);
            kmerMap[seqStart + seqStartTmp] = kmerSeq;

            seqStartTmp += 1;
        }

        return kmerMap;
    }

    // build_minimizer_index
    // kmerVector, minimizer����, minimizer����
    int build_minimizer_index(const map<long long int,string> & kmerMap, 
                              const long long int sequenceLen, 
                              const string chromosome, 
                              map<string,vector<minimizerSrt>> & minimizerMap, 
                              const int minimizerK, 
                              const int minimizerW,
                              const int indexLeft,
                              const int indexRight)
    {
        // ����kmerMap������minimizer����
        for (auto iter = kmerMap.begin(); iter != kmerMap.end(); iter++)
        {
            // begin��minimizerW - 1��ʼ������kmer�������Ĵӵ�iter->second.length() - minimizerK����
            // ������minimizerK���ַ���ֹͣ
            int indexTmp;
            if (iter == kmerMap.begin())
            {
                indexTmp = minimizerW - 1;
            }
            else
            {
                indexTmp = iter->second.length() - minimizerK;
            }
            
            for (long long int i = indexTmp; i < iter->second.length() - minimizerK + 1; i++)
            {
                // ��ʱ������
                auto iterTmp = iter;

                // ��ȡkmer��minimizer�Ӵ�
                string possibleMinimizer = iter->second.substr(i, minimizerK);
                
                // ѭ�������minimizerW������������kmer��û�и�minimizer�Ӵ�
                for (int j = 1; j < minimizerW; j++)
                {
                    // ���жϵ�������û��ָ��������ָ��end�Ļ�������ѭ��
                    if (iterTmp == kmerMap.end())
                    {
                        possibleMinimizer.clear(); // ���������possibleMinimizer���
                        break;
                    }
                    
                    // �ڸ�kmer��֮��minimizerW��kmer�в����Ƿ�������ִ�����������������
                    if (iterTmp->second.find(possibleMinimizer) == string::npos)
                    {
                        possibleMinimizer.clear(); // ���������possibleMinimizer���
                        break;
                    }

                    // ������еĻ�������������һλ
                    iterTmp++;
                }
                // ���possibleMinimizer��Ϊ�գ�������minimizer
                if (possibleMinimizer.length() > 0)
                {
                    // ѡ���С������
                    string minimizer = sequence_select(possibleMinimizer);

                    if (minimizer.length() == 0) // ��������˿��ַ�����������ATGC����ַ�����������minimizer
                    {
                        continue;
                    }

                    // ��ʱ�ṹ��
                    minimizerSrt minimizerSrtTmp;

                    // ��minimizer����ֹλ�ú�Ⱦɫ���
                    minimizerSrtTmp._end = iter->first + i + minimizerK - 1;
                    minimizerSrtTmp._chromosome = chromosome;

                    minimizerMap[minimizer].push_back(minimizerSrtTmp);
                }
            }
        }

        // ��ʱ�ṹ�壬�洢minimizer��Ϣ
        minimizerSrt minimizerSrtTmp;
        auto iterTmp = kmerMap.begin();

        // ��ͷ
        string minimizer = iterTmp->second.substr(0, minimizerK);
        minimizer = sequence_select(minimizer);
        if (minimizer.length() > 0)
        {
            minimizerSrtTmp._end = minimizerK;
            minimizerSrtTmp._chromosome = "chr1";
            minimizerMap[minimizer].push_back(minimizerSrtTmp);
        }

        // ��β
        iterTmp = kmerMap.end();
        iterTmp--;
        minimizer = iterTmp->second.substr(iterTmp->second.length()-minimizerK, minimizerK);
        minimizer = sequence_select(minimizer);

        if (minimizer.length() > 0)
        {
            minimizerSrtTmp._end = sequenceLen;
            minimizerSrtTmp._chromosome = "chr1";
            minimizerMap[minimizer].push_back(minimizerSrtTmp);
        }

        return 0;
    }

    // ���з�ת�ͱȽ�
    string sequence_select(string sequence)
    {
        string sequenceRev;
        
        // ���з�ת�ͻ���
        string tmpSeq = sequence;
        reverse(tmpSeq.begin(), tmpSeq.end());

        for (int i = 0; i < tmpSeq.size(); i++)
        {
            if (tmpSeq[i] == 'a')
            {
                sequenceRev += "t";
            }
            else if (tmpSeq[i] == 't')
            {
                sequenceRev += "a";
            }
            else if (tmpSeq[i] == 'g')
            {
                sequenceRev += "c";
            }
            else if (tmpSeq[i] == 'c')
            {
                sequenceRev += "g";
            }
            else if (tmpSeq[i] == 'A')
            {
                sequenceRev += "T";
            }
            else if (tmpSeq[i] == 'T')
            {
                sequenceRev += "A";
            }
            else if (tmpSeq[i] == 'G')
            {
                sequenceRev += "C";
            }
            else if (tmpSeq[i] == 'C')
            {
                sequenceRev += "G";
            }
            else // �����atgcATGC����ַ���������ѭ�������ؿ��ַ���
            {
                sequenceRev = "";
                return sequenceRev;
            }
        }

        return min(sequence, sequenceRev);
    }

    // ?(A)=0, ?(C)=1, ?(G)=2, ?(T)=3
    // hash���� (E.g. ?(ACATAC) = ?(A)*4^5 + ?(C)*4^4 + ?(A)*4^3 + ?(T)*4^2 + ?(A)*4^1 + ?(C)*4^0 = 305)
    long long int hash_fun(string sequence)
    {
        long long int outIndex = 0;
        map<char, int> base2num = 
        {
            {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
        };

        int tmpIndex = sequence.size() - 1;
        
        for (auto it : sequence)
        {
            outIndex += base2num[it]*pow(4,tmpIndex);
            tmpIndex--;
        }

        return outIndex;
    }
}

#endif