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

// kseq.h 打开文件
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
     * @param inputFasta     需要构建索引的fasta文件
     * 
     * @return fastaMap      map<chr, sequence>
	**/
    map<string, string> build_fasta_index(string inputFasta)
    {
        map<string, string> fastaMap;

        // open fasta file
        gzFile gzfp = gzopen(inputFasta.c_str(), "rb");

        // 打开文件
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
                // ks->name.s 记录的是名字
                // ks->seq.s 记录的是序列
                string chromosome = ks->name.s;
                string sequence = ks->seq.s;

                // 构建fasta索引
                fastaMap[chromosome] = sequence;
            }

            // 释放内存，关闭文件
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        return fastaMap;
    }

    // 打开vcf文件
    void vcf_construct(string inputVcf, map<string, string> & fastaMap, 
                        map<string, map<long long int, vector<map<string, string>>>> & vcfGraph)
    {
        // open base fastafile
        gzFile gzfp = gzopen(inputVcf.c_str(), "rb");

        // 记录上一个节点的位置
        long int tmpRefEnd;
        string tmpChromosome;

        // 打开文件
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

                if (tmpChar != '\n') // 如果没有读到最后一个字符，则一直相加。
                {
                    information += tmpChar;
                }
                else // 否则进行构图
                {
                    // 先判断是不是注释行，是的话跳过
                    if (information.find("#") != string::npos)
                    {
                        // 清空字符串
                        information.clear();
                        string().swap(information);
                    }
                    else // 构图
                    {
                        vector<string> informationVec = split(information, "\t");

                        if (informationVec.size() < 9) // 先检查文件对不对，不对则跳出代码
                        {
                            cerr << "[" << getTime() << "] " 
                                << inputVcf 
                                << ": File format error." 
                                << " [" << information << "]"
                                << endl;
                            exit(1);
                        }

                        // 变异的各种信息
                        string chromosome = informationVec[0];
                        long long int refStart = stol(informationVec[1]);
                        string refSeq = informationVec[3];
                        long long int refEnd = refStart + refSeq.size() - 1;
                        string qrySeq = informationVec[4];

                        // 检查vcf的qry_seq和基因组对应不对应
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

                        // 图中的节点数
                        long long int nodeNum = vcfGraph.size();

                        // 构建vcf图形
                        map<string, string> tmpMap; // 序列和基因型信息
                        
                        // 添加 genotype=0 的节点
                        // 如果染色体第一个vcf位置不是1，则把基因组前边的序列构建第一个node
                        string genotype = "0";
                        if (chromosome != tmpChromosome && refStart > 0)
                        {
                            // 添加新染色体的前部分
                            string preRefSeq;
                            long long int preRefStart;
                            long long int preRefEnd;
                            long long int preRefLen;
                            
                            // 添加旧染色体的后半部分
                            if (nodeNum > 1 && tmpRefEnd < fastaMap[tmpChromosome].length()) // 判断最后一个变异在不在染色体最后
                            {
                                preRefStart = tmpRefEnd + 1;
                                preRefEnd = fastaMap[tmpChromosome].length();
                                preRefLen = preRefEnd - preRefStart + 1;

                                preRefSeq = fastaMap[tmpChromosome].substr(tmpRefEnd, preRefLen);
                                
                                tmpMap["sequence"] = preRefSeq;
                                tmpMap["genotype"] = genotype;
                                vcfGraph[tmpChromosome][preRefStart].push_back(tmpMap);
                            }

                            // 添加新染色体的前部分
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
                        else // 否则把vcf中间的序列构建node
                        {
                            string preRefSeq;
                            long long int preRefStart;
                            long long int preRefEnd;
                            long long int preRefLen;
                            
                            preRefStart = tmpRefEnd + 1;
                            preRefEnd = refStart - 1;
                            preRefLen = preRefEnd - preRefStart + 1;

                            // 判断是不是相邻，相邻的话preRefLen小于等于0，跳过
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
                        
                        // 添加vcf节点
                        // 先添加ref的节点
                        tmpMap["sequence"] = refSeq;
                        tmpMap["genotype"] = genotype;
                        vcfGraph[chromosome][refStart].push_back(tmpMap);

                        // qry可能含有多个等位基因，按，拆分
                        vector<string> qrySeqVec = split(qrySeq, ",");
                        for (auto it : qrySeqVec)
                        {
                            genotype = to_string(stoi(genotype) + 1);
                            tmpMap["sequence"] = it;
                            tmpMap["genotype"] = genotype;
                            vcfGraph[chromosome][refStart].push_back(tmpMap);
                        }
                    }

                    // 清空字符串
                    information.clear();
                    string().swap(information);
                }
            }

            // 最后一个染色体的后半部分
            if (tmpRefEnd < fastaMap[tmpChromosome].length()) // 判断最后一个变异在不在染色体最后
            {
                long long int preRefStart = tmpRefEnd + 1;
                long long int preRefEnd = fastaMap[tmpChromosome].length();
                long long int preRefLen = preRefEnd - preRefStart + 1;

                string preRefSeq = fastaMap[tmpChromosome].substr(tmpRefEnd, preRefLen);

                map<string, string> tmpMap; // 存储起始，终止，序列等信息
                
                tmpMap["sequence"] = preRefSeq;
                tmpMap["genotype"] = "0";
                vcfGraph[tmpChromosome][preRefStart].push_back(tmpMap);
            }

            // 关闭文件，释放内存
            gzclose(gzfp);
        }
    }

    // 为vcfGraph构建索引
    void graphIndex(const map<string, map<long long int, vector<map<string, string>>>> & vcfGraph,
               int32_t kmerLen, 
               int minimizerK, 
               int minimizerW,
               int threads_num)
    {
        // 对应染色体的map
        map<long long int, vector<map<string, string>>> startMapTmp;
        
        for (auto it1 : vcfGraph)
        {
            startMapTmp = it1.second;

            for (map<long long int, vector<map<string, string>>>::iterator iter = startMapTmp.begin();
                 iter != startMapTmp.end(); iter++)
            {
                // 节点的起始位置
                long long int nodeStart = iter->first;

                // 临时迭代器，用于加上一个node的结尾序列
                auto iterTmp = iter;

                for (auto it2 : iter->second)
                {
                    string sequence = it2["sequence"];
                    string genotype = it2["genotype"];

                    // 找上下游一个read的序列
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

    // 找节点前后 n bp的序列，用于构建minimizer
    // find_node_up_down_seq(this node seq, this node GT, this node iterator)
    pair<string, string> find_node_up_down_seq(const string genotype, 
                                               const long long int seqLen,
                                               const map<long long int, vector<map<string, string>>>::iterator & iter, 
                                               const map<long long int, vector<map<string, string>>> & startMapTmp)
    {
        // 记录节点序列
        string upSeq = "";
        string downSeq = "";

        // 临时迭代器
        auto iterTmp = iter;

        // 迭代加上一个节点的末尾序列
        while (upSeq.size() < seqLen && iterTmp != startMapTmp.begin())
        {
            // 迭代器前移
            iterTmp--;

            string seqTmp = "";
            string refSeqTmp = "";

            // 对该节点的所有序列进行遍历
            for (auto it : iterTmp->second)
            {
                // 如果和该节点的基因型一致，则给sequenceTmp赋值
                if (it["genotype"] == genotype)
                {
                    seqTmp = it["sequence"];
                }
                if (it["genotype"] == "0")
                {
                    refSeqTmp = it["sequence"];
                }
            }

            // 判断有没有和该节点的基因型一致的map，没有的话给sequenceTmp赋值ref的序列
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

        // 重置迭代器
        iterTmp = iter;

        // 迭代加下一个节点的开始序列，迭代器后移加到判断条件里，
        while ((downSeq.size() < seqLen) && (++iterTmp != startMapTmp.end()))
        {
            string seqTmp = "";
            string refSeqTmp = "";

            // 对该节点的所有序列进行遍历
            for (auto it : iterTmp->second)
            {
                // 如果和该节点的基因型一致，则给sequenceTmp赋值
                if (it["genotype"] == genotype)
                {
                    seqTmp = it["sequence"];
                }
                if (it["genotype"] == "0")
                {
                    refSeqTmp = it["sequence"];
                }
            }

            // 判断有没有和该节点的基因型一致的map，没有的话给sequenceTmp赋值ref的序列
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
        // 转大写
        transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

        // 存储结果，map<start, sequence>
        map<long long int,string> kmerMap;

        int seqLen = sequence.size();
        int kmerNumber = seqLen - kmerLen + 1;

        int seqStartTmp = 0;

        for (int i = 0; i < kmerNumber; i++)
        {
            // 构建hash表
            string kmerSeq = sequence.substr(seqStartTmp, kmerLen);
            kmerMap[seqStart + seqStartTmp] = kmerSeq;

            seqStartTmp += 1;
        }

        return kmerMap;
    }

    // build_minimizer_index
    // kmerVector, minimizer长度, minimizer数量
    int build_minimizer_index(const map<long long int,string> & kmerMap, 
                              const long long int sequenceLen, 
                              const string chromosome, 
                              map<string,vector<minimizerSrt>> & minimizerMap, 
                              const int minimizerK, 
                              const int minimizerW,
                              const int indexLeft,
                              const int indexRight)
    {
        // 遍历kmerMap并构建minimizer索引
        for (auto iter = kmerMap.begin(); iter != kmerMap.end(); iter++)
        {
            // begin从minimizerW - 1开始索引该kmer，其它的从到iter->second.length() - minimizerK索引
            // 倒数第minimizerK个字符后停止
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
                // 临时迭代器
                auto iterTmp = iter;

                // 提取kmer的minimizer子串
                string possibleMinimizer = iter->second.substr(i, minimizerK);
                
                // 循环看后边minimizerW（包括自身）个kmer有没有该minimizer子串
                for (int j = 1; j < minimizerW; j++)
                {
                    // 先判断迭代器有没有指向最后，如果指向end的话，跳出循环
                    if (iterTmp == kmerMap.end())
                    {
                        possibleMinimizer.clear(); // 如果不是则将possibleMinimizer清空
                        break;
                    }
                    
                    // 在该kmer和之后minimizerW个kmer中查找是否包含该字串，不包含表明不是
                    if (iterTmp->second.find(possibleMinimizer) == string::npos)
                    {
                        possibleMinimizer.clear(); // 如果不是则将possibleMinimizer清空
                        break;
                    }

                    // 如果含有的话，迭代器后移一位
                    iterTmp++;
                }
                // 如果possibleMinimizer不为空，代表是minimizer
                if (possibleMinimizer.length() > 0)
                {
                    // 选择较小的序列
                    string minimizer = sequence_select(possibleMinimizer);

                    if (minimizer.length() == 0) // 如果返回了空字符串，代表有ATGC外的字符，则跳过该minimizer
                    {
                        continue;
                    }

                    // 临时结构体
                    minimizerSrt minimizerSrtTmp;

                    // 该minimizer的终止位置和染色体号
                    minimizerSrtTmp._end = iter->first + i + minimizerK - 1;
                    minimizerSrtTmp._chromosome = chromosome;

                    minimizerMap[minimizer].push_back(minimizerSrtTmp);
                }
            }
        }

        // 临时结构体，存储minimizer信息
        minimizerSrt minimizerSrtTmp;
        auto iterTmp = kmerMap.begin();

        // 开头
        string minimizer = iterTmp->second.substr(0, minimizerK);
        minimizer = sequence_select(minimizer);
        if (minimizer.length() > 0)
        {
            minimizerSrtTmp._end = minimizerK;
            minimizerSrtTmp._chromosome = "chr1";
            minimizerMap[minimizer].push_back(minimizerSrtTmp);
        }

        // 结尾
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

    // 序列反转和比较
    string sequence_select(string sequence)
    {
        string sequenceRev;
        
        // 序列反转和互补
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
            else // 如果有atgcATGC外的字符，则跳出循环并返回空字符串
            {
                sequenceRev = "";
                return sequenceRev;
            }
        }

        return min(sequence, sequenceRev);
    }

    // ?(A)=0, ?(C)=1, ?(G)=2, ?(T)=3
    // hash函数 (E.g. ?(ACATAC) = ?(A)*4^5 + ?(C)*4^4 + ?(A)*4^3 + ?(T)*4^2 + ?(A)*4^1 + ?(C)*4^0 = 305)
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