#ifndef alignment_hpp
#define alignment_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <regex>
#include <map>
#include <list>
#include <unordered_map>
#include <cmath>
#include <iomanip>
#include <getopt.h>
#include <limits.h>
#include <malloc.h>
#include "NeedlemanWunsch.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "zlib.h"
#include "ThreadPool.hpp"
#include "split.hpp"
#include "minimizer.hpp"
#include "sequence_reverse.hpp"

using namespace std;

// kseq.h 打开文件
KSEQ_INIT(gzFile, gzread)


bool sort_tup_by_sec(tuple<string, long long int, int> & a, 
                     tuple<string, long long int, int> & b)
{
    return (get<1>(a) < get<1>(b));
}

namespace alignment
{
    //**********************************************************************//

    // fasta多线程
    void read_fasta_run(const uint64_t chromosomeId, 
                        string sequence, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW,
                        unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap);

    // p_read_fastq多线程
    void p_read_fastq_run(const string chromosome, 
                        const string sequence1, 
                        const string sequence2, 
                        const string quality1, 
                        const string quality2, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW, 
                        const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                        const unordered_map<uint64_t, string> & refFastaMap, 
                        const unordered_map<uint64_t, string> & chromosomeIdMap);

    // s_read_fastq多线程
    void s_read_fastq_run(const string chromosome, 
                        const string sequence, 
                        const string  quality, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW,
                        const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                        const unordered_map<uint64_t, string> & refFastaMap, 
                        const unordered_map<uint64_t, string> & chromosomeIdMap);

    // 寻找公共子串 map<chromosome, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
    map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seeding(const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                                                                                                            const unordered_map<uint64_t, string> & chromosomeIdMap, 
                                                                                                            const unordered_map<uint64_t, list<uint64_t>> & minimizerQryMap, 
                                                                                                            const unsigned int & minimizerK);

    // 计算chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<unsigned long long int, unsigned long long int, unsigned int, uint64_t> & seedingA,
                                          const tuple<unsigned long long int, unsigned long long int, unsigned int, uint64_t> & seedingB, 
                                          const unsigned int & minimizerK);
    
    // chaining -> 迭代法 tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chaining(
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> & seedingChrMap, 
        const unsigned int & minimizerK
        );

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair);

    // 序列比对
    vector<string> run_alignment(const tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair, 
                                 const unordered_map<uint64_t, string> & refFastaMap, 
                                 const unordered_map<uint64_t, string> & chromosomeIdMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality);

    // 将比对的结果转为cigra
    string alignment_cigar(const string & refAliSeq, 
                            const string & midAliSeq,
                            const string & qryAliSeq,
                            const unsigned long long int & softClipNumLeft,
                            const unsigned long long int & softClipNumRight,
                            const unsigned long long int & hardClipNumLeft,
                            const unsigned long long int & hardClipNumRight);

    //**********************************************************************//





    /** 读取fasta文件 
     * @param fastaFile  fasta文件
     * @param minimizerK  minimizer长度
     * @param minimizerW  minimizer窗口大小
     * @param threads  线程数
     * @param minimizerRefMap  minimizer输出哈希表
     *          unordered_map<uint64_t, list<uint64_t>>
     *                        (uint64_t) -> kMer<<8 | kmerSpan
     *                        (uint64_t) -> chromosomeId<<32 | lastPos<<1 | strand
     * @output get<1>(tuple) -> chromosomeIdMap
     * @output get<2>(tuple) -> fastaMap
     */
    tuple<unordered_map<uint64_t, string>, unordered_map<uint64_t, string>> read_fasta
        (const string & fastaFile, 
         const unsigned int & minimizerK, 
         const unsigned int & minimizerW, 
         const uint32_t & threads, 
         unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap)
    {
        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Building index." 
             << endl;

        // 进程池
        ThreadPool pool(threads);

        // 初始化线程池
        pool.init();

        // 输出序列哈希表
        unordered_map<uint64_t, string> chromosomeIdMap;
        unordered_map<uint64_t, string> fastaMap;

        // 初始化染色体id
        uint64_t chromosomeId = 0;

        // open fasta file
        gzFile gzfp = gzopen(fastaFile.c_str(), "rb");

        // 打开文件
        if(!gzfp)
        {
            cerr << "[" << getTime() << "] [alignment] " 
                << fastaFile 
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
                // 转大写
                transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

                pool.submit(read_fasta_run, 
                            chromosomeId, 
                            sequence, 
                            minimizerK,
                            minimizerW,
                            ref(minimizerRefMap));

                // 构建fasta索引
                fastaMap[chromosomeId] = sequence;
                chromosomeIdMap[chromosomeId] = chromosome;

                // 更新染色体ID
                chromosomeId++;

                // 清空字符串，释放内存
                chromosome.clear();
                sequence.clear();
                string().swap(chromosome);
                string().swap(sequence);

                // 检查任务队列是否超过阈值，超过了等待，以防把数据一次性加载到内存中
                while (pool.get_queue() >= threads*10)
                {
                    // 每隔0.5秒检查一次
                    sleep(0.5);
                }
            }

            // 释放内存，关闭文件
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        // 检查任务队列是否执行完，执行完则关闭线程池，否则每隔0.5s检查一次
        while (pool.get_queue() > 0)
        {
            // 每隔0.5秒检查一次
            sleep(0.5);
        }

        // 关闭线程池
        pool.shutdown();

        // log        
        cerr << "[" << getTime() << "] [alignment] " 
                 << "Index build succeeded." 
                 << endl;

        return make_tuple(fastaMap, chromosomeIdMap);
    }

    // fasta多线程
    void read_fasta_run(const uint64_t chromosomeId, 
                        string sequence, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW,
                        unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap)
    {
        // 构建minimizer索引
        mm_sketch(chromosomeId, sequence, minimizerK, minimizerW, minimizerRefMap);
    }

    // 读取fastq文件 双端测序
    void p_read_fastq(const string & fastqFile1, 
                    const string & fastqFile2, 
                    const unsigned int & minimizerK, 
                    const unsigned int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                    const unordered_map<uint64_t, string> & refFastaMap, 
                    const unordered_map<uint64_t, string> & chromosomeIdMap)
    {
        // 进程池
        ThreadPool pool(threads);

        // 初始化线程池
        pool.init();

        // open fastq1 file
        gzFile gzfp1 = gzopen(fastqFile1.c_str(), "rb");

        // open fastq2 file
        gzFile gzfp2 = gzopen(fastqFile2.c_str(), "rb");

        // 打开文件
        if(!gzfp1 || !gzfp2)
        {
            cerr << "[" << getTime() << "] [alignment] " 
                 << fastqFile1 << "/" << fastqFile2
                << ": No such file or directory." 
                << endl;
            exit(1);
        }
        else
        {
            kseq_t *ks1;
            ks1 = kseq_init(gzfp1);
            kseq_t *ks2;
            ks2 = kseq_init(gzfp2);

            while( kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0 )
            {
                // ks->name.s 记录的是名字
                // ks->seq.s 记录的是序列
                // ks->qual.s 记录的是测序质量
                string chromosome1 = ks1->name.s;
                string sequence1 = ks1->seq.s;
                string quality1 = ks1->qual.s;

                string chromosome2 = ks2->name.s;
                string sequence2 = ks2->seq.s;
                string quality2 = ks2->qual.s;

                // 转大写
                transform(sequence1.begin(),sequence1.end(),sequence1.begin(),::toupper);
                transform(sequence2.begin(),sequence2.end(),sequence2.begin(),::toupper);

                string chromosome = split(chromosome1, "/")[0];

                pool.submit(p_read_fastq_run, 
                            chromosome, 
                            sequence1, 
                            sequence2, 
                            quality1, 
                            quality2, 
                            minimizerK, 
                            minimizerW, 
                            ref(minimizerRefMap), 
                            ref(refFastaMap), 
                            ref(chromosomeIdMap));

                // 清空字符串，释放内存
                chromosome1.clear();
                sequence1.clear();
                quality1.clear();
                chromosome2.clear();
                sequence2.clear();
                quality2.clear();
                chromosome.clear();
                string().swap(chromosome1);
                string().swap(sequence1);
                string().swap(quality1);
                string().swap(chromosome2);
                string().swap(sequence2);
                string().swap(quality2);
                string().swap(chromosome);

                // 检查任务队列是否超过阈值，超过了等待，以防把数据一次性加载到内存中
                while (pool.get_queue() >= threads*10)
                {
                    // 每隔0.5秒检查一次
                    sleep(0.5);
                }
            }

            // 释放内存，关闭文件
            kseq_destroy(ks1);
            gzclose(gzfp1);
            kseq_destroy(ks2);
            gzclose(gzfp2);
        }

        // 检查任务队列是否执行完，执行完则关闭线程池，否则每隔0.5s检查一次
        while (pool.get_queue() > 0)
        {
            // 每隔0.5秒检查一次
            sleep(0.5);
        }

        // 关闭线程池
        pool.shutdown();
    }

    // p_read_fastq多线程
    void p_read_fastq_run(const string chromosome, 
                          const string sequence1, 
                          const string sequence2, 
                          const string quality1, 
                          const string quality2, 
                          const unsigned int & minimizerK, 
                          const unsigned int & minimizerW, 
                          const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                          const unordered_map<uint64_t, string> & refFastaMap, 
                          const unordered_map<uint64_t, string> & chromosomeIdMap)
    {
        // 记录minimizer索引
        unordered_map<uint64_t, list<uint64_t>> minimizerQryMap1;
        unordered_map<uint64_t, list<uint64_t>> minimizerQryMap2;
                          
        // 构建minimizer索引
        mm_sketch(0, sequence1, minimizerK, minimizerW, minimizerQryMap1);
        mm_sketch(0, sequence2, minimizerK, minimizerW, minimizerQryMap2);

        // seeding map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap1 = seeding(minimizerRefMap, chromosomeIdMap, minimizerQryMap1, minimizerK);
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap2 = seeding(minimizerRefMap, chromosomeIdMap, minimizerQryMap2, minimizerK);

        // 释放内存
        minimizerQryMap1.clear();
        unordered_map<uint64_t, list<uint64_t>>().swap(minimizerQryMap1);
        minimizerQryMap2.clear();
        unordered_map<uint64_t, list<uint64_t>>().swap(minimizerQryMap2);
   
        // chaining and alignment
        if (seedingMap1.size() > 0)
        {
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chainingOutPair1 = chaining(seedingMap1, minimizerK);
            
            // 释放内存
            seedingMap1.clear();
            map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>>().swap(seedingMap1);

            if (get<0>(chainingOutPair1).size() > 0)
            {
                anchor_merge(chainingOutPair1);
                vector<string> alignmentVec = run_alignment(chainingOutPair1, 
                                                            refFastaMap, 
                                                            chromosomeIdMap, 
                                                            chromosome, 
                                                            sequence1,
                                                            quality1);
            }

            // 释放内存
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int>().swap(chainingOutPair1);
        }

        if (seedingMap2.size() > 0)
        {
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chainingOutPair2 = chaining(seedingMap2, minimizerK);
            
            // 释放内存
            seedingMap2.clear();
            map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>>().swap(seedingMap2);

            if (get<0>(chainingOutPair2).size() > 0)
            {
                anchor_merge(chainingOutPair2);
                vector<string> alignmentVec = run_alignment(chainingOutPair2, 
                                                            refFastaMap, 
                                                            chromosomeIdMap, 
                                                            chromosome, 
                                                            sequence1,
                                                            quality1);
            }

            // 释放内存
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int>().swap(chainingOutPair2);
        }
    }

    // 读取fastq文件  单端测序
    void s_read_fastq(const string & fastqFile, 
                    const unsigned int & minimizerK, 
                    const unsigned int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                    const unordered_map<uint64_t, string> & refFastaMap, 
                    const unordered_map<uint64_t, string> & chromosomeIdMap)
    {
        // 进程池
        ThreadPool pool(threads);

        // 初始化线程池
        pool.init();

        // open fastq1 file
        gzFile gzfp = gzopen(fastqFile.c_str(), "rb");

        // 打开文件
        if(!gzfp)
        {
            cerr << "[" << getTime() << "] [alignment] " 
                 << fastqFile
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
                // map<"sequence1"/"quality1", sequence/quality>
                map<string, string> tmpMap;

                // ks->name.s 记录的是名字
                // ks->seq.s 记录的是序列
                // ks->qual.s 记录的是测序质量
                string chromosome = ks->name.s;
                string sequence = ks->seq.s;
                string quality = ks->qual.s;

                // 转大写
                transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

                pool.submit(s_read_fastq_run, 
                            chromosome, 
                            sequence,  
                            quality, 
                            minimizerK, 
                            minimizerW, 
                            ref(minimizerRefMap), 
                            ref(refFastaMap), 
                            ref(chromosomeIdMap));

                // 清空字符串，释放内存
                chromosome.clear();
                sequence.clear();
                quality.clear();
                string().swap(chromosome);
                string().swap(sequence);
                string().swap(quality);

                // 检查任务队列是否超过阈值，超过了等待，以防把数据一次性加载到内存中
                while (pool.get_queue() >= threads*10)
                {
                    // 每隔0.5秒检查一次
                    sleep(0.5);
                }
            }

            // 释放内存，关闭文件
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        // 检查任务队列是否执行完，执行完则关闭线程池，否则每隔0.5s检查一次
        while (pool.get_queue() > 0)
        {
            // 每隔0.5秒检查一次
            sleep(0.5);
        }

        // 关闭线程池
        pool.shutdown();
    }

    // s_read_fastq多线程
    void s_read_fastq_run(const string chromosome, 
                          const string sequence, 
                          const string quality, 
                          const unsigned int & minimizerK, 
                          const unsigned int & minimizerW,
                          const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                          const unordered_map<uint64_t, string> & refFastaMap, 
                          const unordered_map<uint64_t, string> & chromosomeIdMap)
    {
        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Building index: " << chromosome
             << endl;
        
        // 记录minimizer索引
        unordered_map<uint64_t, list<uint64_t>> minimizerQryMap;

        // 构建minimizer索引
        mm_sketch(0, sequence, minimizerK, minimizerW, minimizerQryMap);

        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Seeding: " << chromosome
             << endl;

        // seeding map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap = seeding(minimizerRefMap, chromosomeIdMap, minimizerQryMap, minimizerK);

        // 释放内存
        minimizerQryMap.clear();
        unordered_map<uint64_t, list<uint64_t>>().swap(minimizerQryMap);
   
        // chaining and alignment
        if (seedingMap.size() > 0)
        {
            // log
            cerr << "[" << getTime() << "] [alignment] " 
                << "Chaining: " << chromosome
                << endl;

            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chainingPair = chaining(seedingMap, minimizerK);
           
            // 释放内存
            seedingMap.clear();
            map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>>().swap(seedingMap);
            
            if (get<0>(chainingPair).size() > 0)
            {
                anchor_merge(chainingPair);

                // log
                cerr << "[" << getTime() << "] [alignment] " 
                    << "Run_alignment: " << chromosome
                    << endl;

                vector<string> alignmentVec = run_alignment(chainingPair, 
                                                            refFastaMap, 
                                                            chromosomeIdMap,  
                                                            chromosome, 
                                                            sequence,
                                                            quality);
                alignmentVec.clear();
                vector<string>().swap(alignmentVec);
            }

            // log
                cerr << "[" << getTime() << "] [alignment] " 
                    << "Done: " << chromosome
                    << endl;

            // 释放内存
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int>().swap(chainingPair);
        }
        malloc_trim(0);	// 0 is for heap memory
    }


    // 寻找公共子串 map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
    map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seeding (
        const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
        const unordered_map<uint64_t, string> & chromosomeIdMap, 
        const unordered_map<uint64_t, list<uint64_t>> & minimizerQryMap, 
        const unsigned int & minimizerK
        )
    {
        // map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap;

        // 遍历qry的minimizer哈希表，寻找hit
        for (auto iter1 = minimizerQryMap.begin(); iter1 != minimizerQryMap.end(); iter1++)
        {
            // 在ref的哈希表中寻找对应的seeding
            auto iter2 = minimizerRefMap.find(iter1->first);
            if (iter2 != minimizerRefMap.end())
            {
                // 遍历qry的list
                for (auto iter3 = iter1->second.begin(); iter3 != iter1->second.end(); iter3++)
                {
                    uint64_t qryChromosomeId = *iter3>>32;
                    unsigned long long int qryEnd = (*iter3 - (qryChromosomeId<<32))>>1;
                    int qryStrand = *iter3 - (qryChromosomeId<<32) - (qryEnd<<1);

                    // 先判断是不是重复序列，是的话跳过该list
                    if (iter2->second.size() > 100)
                    {
                        continue;
                    }
                    
                    // 遍历ref的list
                    for (auto iter4 = iter2->second.begin(); iter4 != iter2->second.end(); iter4++)
                    {
                        uint64_t refChromosomeId = *iter4>>32;
                        unsigned long long int refEnd = (*iter4 - (refChromosomeId<<32))>>1;
                        int refStrand = *iter4 - (refChromosomeId<<32) - (refEnd<<1);

                        // 判断两条链方向是否一致
                        int strand = 1;
                        if (qryStrand == refStrand)
                        {
                            strand = 0;
                        }
                        // seedingMap赋值
                        seedingMap[refChromosomeId][qryEnd].push_back(make_tuple(refEnd, minimizerK, strand));
                    }
                }
            }
        }

        return seedingMap;
    }


    // 计算chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> & seedingA,
                                          const tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> & seedingB, 
                                          const unsigned int & minimizerK)
    {
        // 计算gap数量
        unsigned long long int refEndA = get<0>(seedingA);
        unsigned long long int qryEndA = get<1>(seedingA);

        unsigned long long int refEndB = get<0>(seedingB);
        unsigned long long int qryEndB = get<1>(seedingB);

        // seeding中序列的长度
        long long int seqScore = 0;
        if (qryEndB - get<2>(seedingB) + 1 <= qryEndA)
        {
            seqScore = qryEndB - qryEndA;
        }
        else
        {
            seqScore = get<2>(seedingB);
        }

        long long int gapNumber = refEndB - refEndA - (qryEndB - qryEndA);
        gapNumber = abs(gapNumber);

        // gap得分
        float gapScore = 0;

        if (gapNumber != 0)
        {
            // 计算gap罚分：γ(s) = 0.01*|w|*s + 0.5*log2s
            gapScore = -(0.01*minimizerK*gapNumber + 0.5*(log(gapNumber)/log(2)));
        }
        
        int outScore = 0;

        outScore = seqScore + gapScore;

        return make_tuple(outScore, gapNumber);
    }

    // chaining -> 迭代法 tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chaining(
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> & seedingChrMap, // map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        const unsigned int & minimizerK
        )
    {
        // 记录chaining信息
        // map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>
        map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>> chainingTupListMap;

        // 对总的map进行循环，其中iter->first是ref染色体号，iter->second是对应seeding的哈希表
        for (auto iter = seedingChrMap.begin(); iter != seedingChrMap.end(); iter++)
        {
            uint64_t refChromosomeId = iter->first;
            map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>> seedingMap = iter->second;
        
            for (auto iter1 = seedingMap.begin(); iter1 != seedingMap.end(); iter1++)
            {
                // seeding list1   更新list，记录被选择过的节点
                list<tuple<unsigned long long int, unsigned int, int>> seedingList1Tmp = iter1->second;

                // 先检查是不是空的list，是的话跳过该节点
                if (seedingList1Tmp.size() == 0)
                {
                    continue;
                }

                for (auto iter2 = seedingList1Tmp.begin(); iter2 != seedingList1Tmp.end(); iter2++)
                {
                    // 记录临时chaining信息
                    // list<tuple(refEnd, qryEnd, length, strand, refChromosomeId)>
                    list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>> chainingTupListTmp;

                    // this 节点的信息
                    unsigned long long int refEnd = get<0>(*iter2);
                    unsigned long long int qryEnd = iter1->first;
                    int length = get<1>(*iter2);

                    // 记录chaining score
                    float score = length;

                    // 第一个seeding，添加到chainingTupListTmp中
                    chainingTupListTmp.push_back(make_tuple(refEnd, qryEnd, length, get<2>(*iter2), refChromosomeId));

                    // 向后找节点，计算score
                    auto iterTmp1 = iter1;
                    // 迭代器后移
                    iterTmp1++;
                    
                    while (iterTmp1 != seedingMap.end())
                    {
                        tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> seedingA;

                        auto iterTmp2 = chainingTupListTmp.end();
                        iterTmp2--;
                        seedingA = *iterTmp2;

                        // seeding list2   更新list，记录被选择过的节点
                        list<tuple<unsigned long long int, unsigned int, int>> seedingList2Tmp = iterTmp1->second;
                        // 先检查是不是空的list，是的话跳过该节点
                        if (seedingList2Tmp.size() == 0)
                        {
                            // 迭代器后移
                            iterTmp1++;
                            continue;
                        }

                        // 记录得分最高的索引和分值，以及位置信息
                        float bestScore = 0;
                        unsigned long long int bestRefEnd = 0;
                        unsigned long long int bestQryEnd = 0;
                        long long bestLength = 0;
                        int bestStrand = 0;
                        std::list<tuple<unsigned long long int, unsigned int, int>>::iterator bestIter;

                        // 对下一个节点的List进行循环
                        tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> seedingB;

                        for (auto iter2 = seedingList2Tmp.begin(); iter2 != seedingList2Tmp.end(); iter2++)
                        {
                            // 1. 检查refEnd在不在上一个的后边，不在的话跳过该seeding
                            // 2. 检查strand是不是一个方向，不是的话跳过该seeding
                            // 3. 检查seeding之间的ref长度
                            // long long int betweenSeedLenRef = get<1>(*iter2) - get<1>(seedingA);
                            // betweenSeedLenRef = abs(betweenSeedLenRef);
                            if (get<0>(*iter2) < refEnd || get<2>(*iter2) != get<3>(seedingA))
                            {
                                continue;
                            }

                            // 计算两个seeding之间的score
                            seedingB = make_tuple(get<0>(*iter2), 
                                                  iterTmp1->first, 
                                                  get<1>(*iter2), 
                                                  get<2>(*iter2), 
                                                  refChromosomeId);
                            tuple<float, long long int> outTuple= cal_score(seedingA,
                                                                            seedingB, 
                                                                            minimizerK);

                            // 判断该seeding加上后的score得分是否大于阈值，gap数量要小于1000
                            if ((score + get<0>(outTuple) < 0 && get<1>(outTuple) > 1000))
                            {
                                continue;
                            }
                            else
                            {
                                // seeding 得分
                                int outScore = get<0>(outTuple);

                                if (outScore > bestScore)
                                {
                                    bestScore = outScore;
                                    bestRefEnd = get<0>(*iter2);
                                    bestQryEnd = iterTmp1->first;
                                    bestLength = get<1>(*iter2);
                                    bestStrand = get<2>(*iter2);
                                    bestIter = iter2;
                                }
                            }
                        }

                        // 先判断bestRenEnd是不是0，是的话代表没有匹配，跳过该seeding
                        if (bestRefEnd == 0)
                        {
                            // 迭代器后移
                            iterTmp1++;
                            continue;
                        }
                        
                        score += bestScore;
                        chainingTupListTmp.push_back(make_tuple(bestRefEnd, 
                                                                bestQryEnd, 
                                                                bestLength, 
                                                                bestStrand, 
                                                                refChromosomeId));
                        // 删除被用过的节点
                        seedingList2Tmp.erase(bestIter);
                        // 更新list，删除被选择过的节点
                        seedingMap[iterTmp1->first] = seedingList2Tmp;
                        
                        // 迭代器后移
                        iterTmp1++;
                    }

                    // 检查chaining中的seeding数量，以及score，少于阈值的话跳过
                    if (chainingTupListTmp.size() < 3 || score < 40)
                    {
                        continue;
                    }

                    // 输出chaining
                    // 先检查有没有一样的得分，有的话score加1
                    while (chainingTupListMap.find(score) != chainingTupListMap.end())
                    {
                        score += 0.01;
                    }
                    chainingTupListMap[score] = chainingTupListTmp;
                }
            }
        }

        // 先检查哈希表是不是空的
        int mapQ = 0;
        if (chainingTupListMap.size() >= 1)
        {
            auto iter = chainingTupListMap.end();
            iter--;
            int anchorsNum = iter->second.size();
            float scoreA = iter->first;
            float scoreB;

            // 先判断是不是只有一个chaining，是的话scoreB为0
            if (iter == chainingTupListMap.begin())
            {
                scoreB = 0;
            }
            else
            {
                iter--;
                scoreB = iter->first;
            }

            mapQ = 40*(1-scoreB/scoreA)*min(1, anchorsNum/10)*log(scoreA);
        }
        
        return make_tuple(chainingTupListMap, min(60, mapQ));
    }

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosomeId)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair)
    {
        // chainingOutPair的score不用变，将map中相交的anchor合并
        // 遍历chainging哈希表
        for (auto iter = get<0>(chainingPair).begin(); iter != get<0>(chainingPair).end(); iter++)
        {
            for (auto iter1 = iter->second.begin(); iter1 != iter->second.end(); iter1++)
            {
                auto iter1Tmp = iter1;
                iter1Tmp++;
                
                while ((get<0>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<0>(*iter1) && // 第二个的ref起始位置在第一个ref终止位置的前边
                    (get<1>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<1>(*iter1) && // 第二个的qry起始位置在第一个qry终止位置的前边
                    get<1>(*iter1Tmp) - get<1>(*iter1) - (get<0>(*iter1Tmp) - get<0>(*iter1)) == 0 && // 两个anchor之间没有gap
                    iter1Tmp != iter->second.end()) // iter1Tmp不指向哈希表的end
                {
                    // 符合条件的话，更新第一个anchor的坐标
                    get<2>(*iter1) += get<0>(*iter1Tmp) - get<0>(*iter1);
                    get<0>(*iter1) = get<0>(*iter1Tmp);
                    get<1>(*iter1) = get<1>(*iter1Tmp);
                    iter->second.erase(iter1Tmp);
                    iter1Tmp = iter1;
                    iter1Tmp++;
                }
            }
        }
    }

    // 序列比对
    vector<string> run_alignment(const tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair, 
                                 const unordered_map<uint64_t, string> & refFastaMap, 
                                 const unordered_map<uint64_t, string> & chromosomeIdMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality)
    {
        // 记录比对条数
        int aliNum = 0;

        for (auto iter1 = get<0>(chainingPair).rbegin(); iter1 != get<0>(chainingPair).rend(); iter1++)
        {
            // 只输出score前两个
            if (aliNum > 1)
            {
                break;
            }
            
            float score = iter1->first;
            list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>> chainingTupList = iter1->second;

            // 染色体号
            uint64_t refChromosomeId = get<4>(*chainingTupList.begin());

            // chainging方向
            int strand = get<3>(*chainingTupList.begin());
            string qrySeqConvert = qrySeq;

            // 根据strand的方向选择序列
            if (strand == 1)
            {
                qrySeqConvert = sequence_reverse(qrySeq);
            }
            

            // 参考基因组的序列
            string refSeq = "";
            string refChromosome = "";
            auto iterFind1 = refFastaMap.find(refChromosomeId);
            auto iterFind2 = chromosomeIdMap.find(refChromosomeId);
            if (iterFind1 == refFastaMap.end() && iterFind2 == chromosomeIdMap.end())
            {
                cerr << "[" << getTime() << "] [alignment] " 
                    << "Chromosome Id not present in reference genome: " << refChromosomeId 
                    << endl;
                exit(1);
            }

            // 比对结果
            string refAliSeqOut = "";
            string qryAliSeqOut = "";
            string midAliSeqOut = "";

            // 记录比对的坐标
            long long int aliRefEndTmp = 0;
            long long int aliQryEndTmp = 0;

            // 遍历chaining中seeding，进行比对
            for (auto iter = chainingTupList.begin(); iter != chainingTupList.end(); iter++)
            {
                long long int length = get<2>(*iter);
                long long int refEnd = get<0>(*iter);
                long long int refStart = refEnd - length + 1;

                long long int qryEnd = get<1>(*iter);
                long long int qryStart = qryEnd - length + 1;

                // 第一个anchor
                if (iter == chainingTupList.begin())
                {
                    for (size_t i = 0; i < length; i++)
                    {
                        // 加上seeding的序列，不提序列，匹配的地方都用N替代
                        refAliSeqOut += "N";
                        qryAliSeqOut += "N";
                        midAliSeqOut += "|";
                    }

                    // 更新比对的坐标
                    aliRefEndTmp = refEnd;
                    aliQryEndTmp = qryEnd;
                }
                // 不是第一个anchor时
                else
                {
                    long long int aliRefStart = aliRefEndTmp + 1;
                    long long int aliRefEnd = refStart - 1;

                    long long int aliQryStart = aliQryEndTmp + 1;
                    long long int aliQryEnd = qryStart - 1;

                    // 判断变异类型
                    // snp
                    if ((aliRefEnd - aliRefStart == 1) && (aliQryEnd - aliQryStart == 1))
                    {
                        string refAliSeq = iterFind1->second.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeqConvert.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

                        // 加上序列
                        refAliSeqOut += refAliSeq;
                        qryAliSeqOut += qryAliSeq;
                        midAliSeqOut += " ";
                        
                        for (size_t i = 0; i < length; i++)
                        {
                            // 加上seeding的序列，不提序列，匹配的地方都用N替代
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // 更新比对的坐标
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;
                    }
                    // 替换
                    else if((aliRefEnd - aliRefStart > 1) && (aliQryEnd - aliQryStart > 1))
                    {
                        string refAliSeq = iterFind1->second.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeqConvert.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

                        auto alignmentOutTup = alignment::alignment(refAliSeq, qryAliSeq, 9, -3, -2);

                        // 加上序列
                        refAliSeqOut += get<0>(alignmentOutTup);
                        qryAliSeqOut += get<1>(alignmentOutTup);
                        midAliSeqOut += get<2>(alignmentOutTup);
                        
                        for (size_t i = 0; i < length; i++)
                        {
                            // 加上seeding的序列，不提序列，匹配的地方都用N替代
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // 更新比对的坐标
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;
                    }
                    // 判断起始位置是不是大于终止位置，ref上出现
                    // 1_RefStart:2728 1_RefEnd:4304 1_QryStart:1647 1_QryEnd:3223  2_RefStart:4305 2_RefEnd:4328 2_QryStart:3256 2_QryEnd:3279
                    // INS
                    else if ((aliRefStart - aliRefEnd == 1) && aliQryStart <= aliQryEnd)
                    {
                        // qry加上序列
                        qryAliSeqOut += qrySeqConvert.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

                        // ref直接加gap
                       for (size_t i = 0; i < (aliQryEnd - aliQryStart + 1); i++)
                        {
                            refAliSeqOut += "-";
                            midAliSeqOut += " ";
                        }
                        
                        for (size_t i = 0; i < length; i++)
                        {
                            // 加上seeding的序列，不提序列，匹配的地方都用N替代
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // 更新比对的坐标
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;

                        continue;
                    }
                    // qry上的起始位置大于终止位置
                    // 1_RefStart:2674 1_RefEnd:2723 1_QryStart:1597 1_QryEnd:1646  2_RefStart:2728 2_RefEnd:4304 2_QryStart:1647 2_QryEnd:3223
                    // DEL
                    else if ((aliQryStart - aliQryEnd == 1) && aliRefStart <= aliRefEnd)
                    {
                        // qry直接加gap
                       for (size_t i = 0; i < (aliRefEnd - aliRefStart + 1); i++)
                        {
                            qryAliSeqOut += "-";
                            midAliSeqOut += " ";
                        }
                        // ref加上序列
                        refAliSeqOut += iterFind1->second.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);

                        for (size_t i = 0; i < length; i++)
                        {
                            // 加上seeding的序列，不提序列，匹配的地方都用N替代
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // 更新比对的坐标
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;

                        continue;
                    }
                    // Translocation
                    else if (aliRefStart - aliRefEnd > 1)
                    {
                        cerr << "[" << getTime() << "] [alignment] Incorrect coordinates, delete the seeding.\n";
                        continue;
                        // exit(1);
                    }
                    // Translocation
                    else if(aliQryStart - aliQryEnd > 1)
                    {
                        cerr << "[" << getTime() << "] [alignment] Incorrect coordinates, delete the seeding.\n";
                        continue;
                        // exit(1);
                    }
                    else
                    {
                        cerr << "[" << getTime() << "] [alignment] Unknow type.\n";
                        continue;
                    }
                }
            }

            // qry两端序列比对
            auto iter = chainingTupList.begin();
            long long int length = get<2>(*iter);
            long long int refMinStart = get<0>(*iter) - length + 1;
            long long int qryMinStart = get<1>(*iter) - length + 1;
            if (qryMinStart > 0)
            {
                string refAliSeq = iterFind1->second.substr(refMinStart - qryMinStart, qryMinStart - 1);
                string qryAliSeq = qrySeqConvert.substr(0, qryMinStart - 1);
                refMinStart -= qryMinStart - 1;
                qryMinStart = 1;

                auto alignmentOutTup = alignment::alignment(refAliSeq, qryAliSeq, 9, -3, -2);

                // 加上序列
                refAliSeqOut = get<0>(alignmentOutTup) + refAliSeqOut;
                qryAliSeqOut = get<1>(alignmentOutTup) + qryAliSeqOut;
                midAliSeqOut = get<2>(alignmentOutTup) + midAliSeqOut;
            }
            
            iter = chainingTupList.end();
            iter--;
            long long int refMaxEnd = get<0>(*iter);
            long long int qryMaxEnd = get<1>(*iter);
            if (qryMaxEnd < qrySeqConvert.size())
            {
                string refAliSeq = iterFind1->second.substr(refMaxEnd, qrySeqConvert.size() - qryMaxEnd);
                string qryAliSeq = qrySeqConvert.substr(qryMaxEnd, qrySeqConvert.size() - qryMaxEnd);
                refMaxEnd += qrySeqConvert.size() - qryMaxEnd;
                qryMaxEnd = qrySeqConvert.size();

                auto alignmentOutTup = alignment::alignment(refAliSeq, qryAliSeq, 9, -3, -2);

                // 加上序列
                refAliSeqOut += get<0>(alignmentOutTup);
                qryAliSeqOut += get<1>(alignmentOutTup);
                midAliSeqOut += get<2>(alignmentOutTup);
            }
            

            string cigra = alignment_cigar(refAliSeqOut, midAliSeqOut, qryAliSeqOut, 
                                            0, 0, 0, 0);

            string samOut = qryChromosome + "\t" + 
                            "2048\t" + 
                            iterFind2->second + "\t" + 
                            to_string(refMinStart) + "\t" + 
                            to_string(get<1>(chainingPair)) + "\t" + 
                            cigra + "\t" + 
                            "*\t*\t" + 
                            qrySeq + "\t" + 
                            qryQuality;
                            // qrySeq.substr(qryMinStart - 1, qryMaxEnd + 1 - qryMinStart) + "\t" + 
                            // qryQuality.substr(qryMinStart - 1, qryMaxEnd + 1 - qryMinStart);
            cout << samOut << endl;

            // 1. 单端->0   双端->1
            // 2. 正常比对，如果是PE测序，还代表PE的两条read之间的比对距离没有明显的偏离插入片段长度
            // 4. 序列没有mapping到参考序列上
            // 8. 序列的另一端序列没有比对到参考序列上，比如这条序列是R1,它对应的R2端序列没有比对到参考序列上
            // 16. 序列比对到参考序列的负链上
            // 32. 序列对应的另一端序列比对到参考序列的负链上
            // 64. 序列是R1端序列，read1
            // 128. 序列是R2端序列，read2
            // 256. 序列不是主要的比对
            // 512. 该read没有通过质量控制
            // 1024. 序列是PCR重复序列
            // 2048. read可能存在嵌合
            // 将比对结果转为FLAG
            int flagScore = 0;
            if (midAliSeqOut.find("|") != string::npos)
            {
                flagScore += 2048;
            }

            refAliSeqOut.clear();
            midAliSeqOut.clear();
            qryAliSeqOut.clear();

            // 更新比对条数
            aliNum++;
        }
        

        vector<string> alignmentVec;
        malloc_trim(0);	// 0 is for heap memory

        return alignmentVec;
    }

    // 将比对的结果转为cigra
    string alignment_cigar(const string & refAliSeq, 
                           const string & midAliSeq,
                           const string & qryAliSeq,
                           const unsigned long long int & softClipNumLeft,
                           const unsigned long long int & softClipNumRight,
                           const unsigned long long int & hardClipNumLeft,
                           const unsigned long long int & hardClipNumRight)
    {
        // 记录转化的结果
        string cigar;

        // // 先判断是不是完全匹配，是的话直接用=
        // if (midAliSeq.find(" ") == string::npos)
        // {
        //     cigar += to_string(midAliSeq.length()) + "=";
        //     return cigar;
        // }
        
        // 跳过的序列，左
        if (softClipNumLeft != 0)
        {
            cigar += to_string(softClipNumLeft) + "S";
        }
        else if (hardClipNumLeft != 0)
        {
            cigar += to_string(hardClipNumLeft) + "H";
        }
        
        // 记录总的数量
        unsigned long long int matchNum = 0;
        unsigned long long int delNum = 0;
        unsigned long long int insNum = 0;

        // 记录单个的数量
        unsigned long long int matchNumTmp = 0;
        unsigned long long int delNumTmp = 0;
        unsigned long long int insNumTmp = 0;

        // 记录上一个比对的类型
        string typeTmp = "";

        // 循环midAliSeq转换
        for (size_t i = 0; i < midAliSeq.size(); i++)
        {
            // 记录当前的比对类型
            string type;
        
            if (midAliSeq[i] == '|')
            {
                type = "match";
            }
            else
            {
                if (refAliSeq[i] != '-' && qryAliSeq[i] != '-')
                {
                    type = "match";
                }
                else if (refAliSeq[i] != '-' && qryAliSeq[i] == '-')
                {
                    type = "del";
                }
                else if (refAliSeq[i] == '-' && qryAliSeq[i] != '-')
                {
                    type = "ins";
                }
                else
                {
                    cerr << "Alignment type unknown: " << type << endl
                         << refAliSeq << endl 
                         << midAliSeq << endl 
                         << qryAliSeq << endl;
                    exit(1);
                }
            }

            // 比对的类型和上一个比对类型不一样
            if (type != typeTmp)
            {
                if (typeTmp.length() == 0)
                {
                    typeTmp = type;
                }
                else
                {
                    if (typeTmp == "match")
                    {
                        cigar += to_string(matchNumTmp) + "M";
                        matchNumTmp = 0;
                        typeTmp = type;
                    }
                    else if (typeTmp == "del")
                    {
                        cigar += to_string(delNumTmp) + "D";
                        delNumTmp = 0;
                        typeTmp = type;
                    }
                    else if (typeTmp == "ins")
                    {
                        cigar += to_string(insNumTmp) + "I";
                        insNumTmp = 0;
                        typeTmp = type;
                    }
                    else
                    {
                        cerr << "Alignment type unknown: " << typeTmp << endl
                            << refAliSeq << endl 
                            << midAliSeq << endl 
                            << qryAliSeq << endl;
                        exit(1);
                    }
                }
            }

            // 数量增加1
            if (type == "match")
            {
                matchNum++;
                matchNumTmp++;
            }
            else if (type == "del")
            {
                delNum++;
                delNumTmp++;
            }
            else if (type == "ins")
            {
                insNum++;
                insNumTmp++;
            }
            else
            {
                cerr << "Alignment type unknown: " << type << endl
                    << refAliSeq << endl 
                    << midAliSeq << endl 
                    << qryAliSeq << endl;
                exit(1);
            }

            // 最后一个比对情况
            if (i == midAliSeq.size() - 1)
            {
                if (typeTmp == "match")
                {
                    cigar += to_string(matchNumTmp) + "M";
                    matchNumTmp = 0;
                }
                else if (typeTmp == "del")
                {
                    cigar += to_string(delNumTmp) + "D";
                    delNumTmp = 0;
                }
                else if (typeTmp == "ins")
                {
                    cigar += to_string(insNumTmp) + "I";
                    insNumTmp = 0;
                }
                else
                {
                    cerr << "Alignment type unknown: " << typeTmp << endl
                        << refAliSeq << endl 
                        << midAliSeq << endl 
                        << qryAliSeq << endl;
                    exit(1);
                }
            }
        }

        // 跳过的序列，右
        if (softClipNumRight != 0)
        {
            cigar += to_string(softClipNumRight) + "S";
        }
        else if (hardClipNumRight != 0)
        {
            cigar += to_string(hardClipNumRight) + "H";
        }

        return cigar;
    }
}

#endif