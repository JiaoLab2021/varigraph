#ifndef alignment_hpp
#define alignment_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <regex>
#include <map>
#include <cmath>
#include <iomanip>
#include <getopt.h>
#include <limits.h>
#include <mutex>
#include "NeedlemanWunsch.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "zlib.h"
#include "ThreadPool.hpp"
#include "split.hpp"
#include <list>
#include <unordered_map>

using namespace std;

// kseq.h 打开文件
KSEQ_INIT(gzFile, gzread)

std::mutex mtx;

// 对tup按qryEnd排序，正序
bool sort_tup_by_sec(tuple<long long int, long long int, int, string> & a, 
                     tuple<long long int, long long int, int, string> & b)
{
    return (get<1>(a) < get<1>(b));
}

namespace alignment
{
    //**********************************************************************//

    // fasta多线程
    void read_fasta_run(string sequence, 
                        const string & chromosome,
                        const uint32_t& kmerLen, 
                        const int & minimizerK, 
                        const int & minimizerW,
                        unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap);

    // p_read_fastq多线程
    void p_read_fastq_run(string sequence1, 
                        string sequence2, 
                        const string & quality1, 
                        const string & quality2, 
                        const string & chromosome,
                        const uint32_t& kmerLen, 
                        const int & minimizerK, 
                        const int & minimizerW, 
                        const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                        const map<string, string> & refFastaMap);

    // s_read_fastq多线程
    void s_read_fastq_run(string sequence, 
                        const string & sequenceQuality, 
                        const string & chromosome,
                        const uint32_t& kmerLen, 
                        const int & minimizerK, 
                        const int & minimizerW,
                        const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                        const map<string, string> & refFastaMap);

    // build kmer index
    map<long long int,string> build_kmer_index(const string & sequence, const uint32_t& kmerLen);

    // build_minimizer_index
    int build_minimizer_index(const map<long long int, string> & kmerMap, 
                            const long long int & sequenceLen, 
                            const string & chromosome, 
                            unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerMap, 
                            const int & minimizerK, 
                            const int & minimizerW);

    // build_minimizer_index
    int build_minimizer_index(const map<long long int, string> & kmerMap, 
                            const long long int & sequenceLen, 
                            const string & chromosome, 
                            map<string, map<long long int, vector<long long int>>> & minimizerMap, 
                            const int & minimizerK, 
                            const int & minimizerW);

    // 寻找公共子串
    map<string, map<long long int, list<tuple<long long int, int>>>> seeding(const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                                                                             const map<string, map<long long int, vector<long long int>>> & minimizerQryMap, 
                                                                             const int & minimizerK);

    // 计算chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<long long int, long long int, int, string> & seedingA,
                                          const tuple<long long int, long long int, int, string> & seedingB, 
                                          const int & minimizerK);
    
    // chaining -> 迭代法 tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chaining(map<string, map<long long int, list<tuple<long long int, int>>>> & seedingChrMap,                                                                           
                                                                                              const string & qryChromosome);

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> & chainingOutPair);

    // 序列比对
    vector<string> run_alignment(const tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> & chainingOutPair, 
                                 const map<string, string> & refFastaMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality);

    // 将比对的结果转为cigra
    string alignment_cigar(const string & refAliSeq, 
                        const string & midAliSeq,
                        const string & qryAliSeq,
                        const long long int & softClipNumLeft,
                        const long long int & softClipNumRight,
                        const long long int & hardClipNumLeft,
                        const long long int & hardClipNumRight);

    // 反向互补
    string sequence_reverse(const string & sequence);

    // hash函数 (E.g. α(ACATAC) = α(A)*4^5 + α(C)*4^4 + α(A)*4^3 + α(T)*4^2 + α(A)*4^1 + α(C)*4^0 = 305)
    long long int hash_fun(const string & sequence);

    //**********************************************************************//





    // 读取fasta文件
    map<string, string> read_fasta(const string & fastaFile, 
                                   const uint32_t& kmerLen, 
                                   const int & minimizerK, 
                                   const int & minimizerW, 
                                   const uint32_t & threads, 
                                   unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap)
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
        map<string, string> fastaMap;

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

                pool.submit(read_fasta_run, 
                            sequence, 
                            chromosome, 
                            kmerLen, 
                            minimizerK,
                            minimizerW,
                            ref(minimizerRefMap));

                // 构建fasta索引
                fastaMap[chromosome] = sequence;

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

        return fastaMap;
    }

    // fasta多线程
    void read_fasta_run(string sequence, 
                        const string & chromosome,
                        const uint32_t& kmerLen, 
                        const int & minimizerK, 
                        const int & minimizerW,
                        unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap)
    {
        // 转大写
        transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

        // 构建kmer索引
        map<long long int,string> kmerMap = build_kmer_index(sequence, kmerLen);
        
        // 构建minimizer索引
        build_minimizer_index(kmerMap, 
                              sequence.length(), 
                              chromosome, 
                              minimizerRefMap, 
                              minimizerK, 
                              minimizerW);
        // 释放内存
        kmerMap.clear();
        map<long long int,string>().swap(kmerMap);
    }

    // 读取fastq文件 双端测序
    void p_read_fastq(const string & fastqFile1, 
                    const string & fastqFile2, 
                    const uint32_t& kmerLen, 
                    const int & minimizerK, 
                    const int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                    const map<string, string> & refFastaMap)
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
                // map<"sequence1/2"/"quality1/2", sequence/quality>
                map<string, string> tmpMap;

                // ks->name.s 记录的是名字
                // ks->seq.s 记录的是序列
                // ks->qual.s 记录的是测序质量
                string chromosome1 = ks1->name.s;
                string sequence1 = ks1->seq.s;
                string quality1 = ks1->qual.s;

                string chromosome2 = ks2->name.s;
                string sequence2 = ks2->seq.s;
                string quality2 = ks2->qual.s;

                string chromosome = split(chromosome1, "/")[0];

                pool.submit(p_read_fastq_run, 
                            sequence1, 
                            sequence2, 
                            quality1, 
                            quality2, 
                            chromosome, 
                            kmerLen, 
                            minimizerK, 
                            minimizerW, 
                            ref(minimizerRefMap), 
                            ref(refFastaMap));

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
    void p_read_fastq_run(string sequence1, 
                          string sequence2, 
                          const string & quality1, 
                          const string & quality2, 
                          const string & chromosome,
                          const uint32_t& kmerLen, 
                          const int & minimizerK, 
                          const int & minimizerW, 
                          const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                          const map<string, string> & refFastaMap)
    {
        // 转大写
        transform(sequence1.begin(),sequence1.end(),sequence1.begin(),::toupper);
        transform(sequence2.begin(),sequence2.end(),sequence2.begin(),::toupper);

        // 构建kmer索引
        map<long long int,string> kmerMap1 = build_kmer_index(sequence1, kmerLen);
        map<long long int,string> kmerMap2 = build_kmer_index(sequence2, kmerLen);
        
        // 记录minimizer索引
        map<string, map<long long int, vector<long long int>>> minimizerQryMap1;
        map<string, map<long long int, vector<long long int>>> minimizerQryMap2;
                          
        // 构建minimizer索引
        build_minimizer_index(kmerMap1, 
                              sequence1.length(), 
                              chromosome, 
                              minimizerQryMap1, 
                              minimizerK, 
                              minimizerW);

        build_minimizer_index(kmerMap2, 
                              sequence2.length(), 
                              chromosome, 
                              minimizerQryMap2, 
                              minimizerK, 
                              minimizerW);
        // 释放内存
        kmerMap1.clear();
        map<long long int,string>().swap(kmerMap1);
        kmerMap2.clear();
        map<long long int,string>().swap(kmerMap2);

        // seeding
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap1 = seeding(minimizerRefMap, minimizerQryMap1, minimizerK);
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap2 = seeding(minimizerRefMap, minimizerQryMap2, minimizerK);

        // 释放内存
        minimizerQryMap1.clear();
        map<string, map<long long int, vector<long long int>>>().swap(minimizerQryMap1);
        minimizerQryMap2.clear();
        map<string, map<long long int, vector<long long int>>>().swap(minimizerQryMap2);
   
        // chaining and alignment
        if (seedingMap1.size() > 0)
        {
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chainingOutPair1 = chaining(seedingMap1, chromosome);
            
            // 释放内存
            seedingMap1.clear();
            map<string, map<long long int, list<tuple<long long int, int>>>>().swap(seedingMap1);

            if (get<0>(chainingOutPair1).size() > 0)
            {
                anchor_merge(chainingOutPair1);
                vector<string> alignmentVec = run_alignment(chainingOutPair1, 
                                                            refFastaMap, 
                                                            chromosome, 
                                                            sequence1,
                                                            quality1);
            }

            // 释放内存
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int>().swap(chainingOutPair1);
        }

        if (seedingMap2.size() > 0)
        {
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chainingOutPair2 = chaining(seedingMap2, chromosome);
            
            // 释放内存
            seedingMap2.clear();
            map<string, map<long long int, list<tuple<long long int, int>>>>().swap(seedingMap2);

            if (get<0>(chainingOutPair2).size() > 0)
            {
                anchor_merge(chainingOutPair2);
                vector<string> alignmentVec = run_alignment(chainingOutPair2, 
                                                            refFastaMap, 
                                                            chromosome, 
                                                            sequence1,
                                                            quality1);
            }

            // 释放内存
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int>().swap(chainingOutPair2);
        }
    }

    // 读取fastq文件  单端测序
    void s_read_fastq(const string & fastqFile, 
                    const uint32_t& kmerLen, 
                    const int & minimizerK, 
                    const int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                    const map<string, string> & refFastaMap)
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

                pool.submit(s_read_fastq_run, 
                            sequence,  
                            quality, 
                            chromosome, 
                            kmerLen, 
                            minimizerK, 
                            minimizerW, 
                            ref(minimizerRefMap), 
                            ref(refFastaMap));

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
    void s_read_fastq_run(string sequence, 
                          const string & quality, 
                          const string & chromosome,
                          const uint32_t& kmerLen, 
                          const int & minimizerK, 
                          const int & minimizerW,
                          const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                          const map<string, string> & refFastaMap)
    {
        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Building index: " << chromosome
             << endl;

        // 转大写
        transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

        // 构建kmer索引
        map<long long int,string> kmerMap = build_kmer_index(sequence, kmerLen);
        
        // 记录minimizer索引
        map<string, map<long long int, vector<long long int>>> minimizerQryMap;

        // 构建minimizer索引
        build_minimizer_index(kmerMap, 
                              sequence.length(), 
                              chromosome, 
                              minimizerQryMap, 
                              minimizerK, 
                              minimizerW);

        // 释放内存
        kmerMap.clear();
        map<long long int,string>().swap(kmerMap);

        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Seeding: " << chromosome
             << endl;

        // seeding
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap = seeding(minimizerRefMap, minimizerQryMap, minimizerK);

        // 释放内存
        minimizerQryMap.clear();
        map<string, map<long long int, vector<long long int>>>().swap(minimizerQryMap);
   
        // chaining and alignment
        if (seedingMap.size() > 0)
        {
            // log
            cerr << "[" << getTime() << "] [alignment] " 
                << "Chaining: " << chromosome
                << endl;

            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chainingOutPair = chaining(seedingMap, chromosome);
           
            // 释放内存
            seedingMap.clear();
            map<string, map<long long int, list<tuple<long long int, int>>>>().swap(seedingMap);

            if (get<0>(chainingOutPair).size() > 0)
            {
                anchor_merge(chainingOutPair);

                // log
                cerr << "[" << getTime() << "] [alignment] " 
                    << "Run_alignment: " << chromosome
                    << endl;

                vector<string> alignmentVec = run_alignment(chainingOutPair, 
                                                            refFastaMap, 
                                                            chromosome, 
                                                            sequence,
                                                            quality);
            }

            // log
                cerr << "[" << getTime() << "] [alignment] " 
                    << "Done: " << chromosome
                    << endl;

            // 释放内存
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int>().swap(chainingOutPair);
        }
    }

    // build kmer index
    map<long long int,string> build_kmer_index(const string & sequence, const uint32_t& kmerLen)
    {
        // kmerMap<start, sequence> 输出结果
        map<long long int,string> kmerMap;

        long long int seqLen = sequence.size();
        long long int kmerNumber = seqLen - kmerLen + 1;

        for (long long int i = 0; i < kmerNumber; i++)
        {
            // 存储kmer到哈希表中
            string kmerSeq = sequence.substr(i, kmerLen);

            // 必须加1，索引从0开始，但是序列是从1开始
            kmerMap[i + 1] = kmerSeq;
        }

        return kmerMap;
    }

    // build_minimizer_index
    int build_minimizer_index(const map<long long int, string> & kmerMap, 
                              const long long int & sequenceLen, 
                              const string & chromosome, 
                              map<string, map<long long int, vector<long long int>>> & minimizerMap, 
                              const int & minimizerK, 
                              const int & minimizerW)
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
                iterTmp++;
                for (int j = 1; j < minimizerW; j++)
                {
                    // 先判断迭代器有没有指向最后，如果指向end的话，跳出循环
                    if (iterTmp == kmerMap.end())
                    {
                        possibleMinimizer.clear(); // possibleMinimizer清空
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
                    // 多线程数据锁
                    std::lock_guard<std::mutex> mtx_locker(mtx);

                    // 计算序列对应的hash值
                    long long int sequenceHashNum = hash_fun(possibleMinimizer);
                    // 染色体 minimizer序列 终止位置
                    minimizerMap[chromosome][sequenceHashNum].push_back(iter->first + i + minimizerK - 1);
                }
            }
        }

        // 开头
        auto iterTmp = kmerMap.begin();
        string minimizer = iterTmp->second.substr(0, minimizerK);
        if (minimizer.length() > 0)
        {
            // 多线程数据锁
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // 计算序列对应的hash值
            long long int sequenceHashNum = hash_fun(minimizer);
            // 染色体 minimizer序列 终止位置
            minimizerMap[chromosome][sequenceHashNum].push_back(minimizerK);
        }

        // 结尾
        iterTmp = kmerMap.end();
        iterTmp--;
        minimizer.clear();
        minimizer = iterTmp->second.substr(iterTmp->second.length()-minimizerK, minimizerK);
        if (minimizer.length() > 0)
        {
            // 多线程数据锁
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // 计算序列对应的hash值
            long long int sequenceHashNum = hash_fun(minimizer);
            // 染色体 minimizer序列 终止位置
            minimizerMap[chromosome][sequenceHashNum].push_back(sequenceLen);
        }
        
        return 0;
    }


    // build_minimizer_index unordered_map
    int build_minimizer_index(const map<long long int, string> & kmerMap, 
                              const long long int & sequenceLen, 
                              const string & chromosome, 
                              unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerMap, 
                              const int & minimizerK, 
                              const int & minimizerW)
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
                iterTmp++;
                for (int j = 1; j < minimizerW; j++)
                {
                    // 先判断迭代器有没有指向最后，如果指向end的话，跳出循环
                    if (iterTmp == kmerMap.end())
                    {
                        possibleMinimizer.clear(); // possibleMinimizer清空
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
                    // 多线程数据锁
                    std::lock_guard<std::mutex> mtx_locker(mtx);

                    // 计算序列对应的hash值
                    long long int sequenceHashNum = hash_fun(possibleMinimizer);
                    // 染色体 minimizer序列 终止位置
                    minimizerMap[chromosome][sequenceHashNum].push_back(iter->first + i + minimizerK - 1);
                }
            }
        }

        // 开头
        auto iterTmp = kmerMap.begin();
        string minimizer = iterTmp->second.substr(0, minimizerK);
        if (minimizer.length() > 0)
        {
            // 多线程数据锁
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // 计算序列对应的hash值
            long long int sequenceHashNum = hash_fun(minimizer);
            // 染色体 minimizer序列 终止位置
            minimizerMap[chromosome][sequenceHashNum].push_back(minimizerK);
        }

        // 结尾
        iterTmp = kmerMap.end();
        iterTmp--;
        minimizer.clear();
        minimizer = iterTmp->second.substr(iterTmp->second.length()-minimizerK, minimizerK);
        if (minimizer.length() > 0)
        {
            // 多线程数据锁
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // 计算序列对应的hash值
            long long int sequenceHashNum = hash_fun(minimizer);
            // 染色体 minimizer序列 终止位置
            minimizerMap[chromosome][sequenceHashNum].push_back(sequenceLen);
        }
        
        return 0;
    }

    // 寻找公共子串 map<chromosome, map<qryEnd, list<tuple<refEnd, length>>>>
    map<string, map<long long int, list<tuple<long long int, int>>>> seeding(const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                                                                             const map<string, map<long long int, vector<long long int>>> & minimizerQryMap, 
                                                                             const int & minimizerK)
    {
        // map<chromosome, map<qryEnd, list<tuple<refEnd, length>>>>
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap;

        for (auto it1 : minimizerQryMap)
        {
            // it1.first为染色体号
            string qryChromosome = it1.first;
            
            // 对ref的minimizerMap进行循环，找与每一条染色体的seedNum
            // it2.first为ref的染色体号， it2.second为ref的 map<sequence, vector<end>>
            for (auto it2 : minimizerRefMap)
            {
                string refChromosome = it2.first;

                // it1.second 为 map<sequence, vector<end>>
                for (auto it3 : it1.second)
                {
                    // it3.first 为minimizer序列
                    // it3.second 为minimizer的终止位置Vector

                    // 在ref的minimizerMap中找索引
                    auto iter1 = it2.second.find(it3.first);
                    
                    if (iter1 != it2.second.end())
                    {
                        // 对应到了ref中多个位置，则分别构建结构体
                        for (auto it4 : iter1->second)
                        {
                            // 对qryEndVec进行循环
                            for (auto it5 : it3.second)
                            {
                                // seedingMap赋值
                                seedingMap[refChromosome][it5].push_back(make_tuple(it4, minimizerK));
                            }
                        }
                    }
                }
            }
        }

        return seedingMap;
    }


    // 计算chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<long long int, long long int, int, string> & seedingA,
                                          const tuple<long long int, long long int, int, string> & seedingB, 
                                          const int & minimizerK)
    {
        // 计算gap数量
        long long int refEndA = get<0>(seedingA);
        long long int qryEndA = get<1>(seedingA);

        long long int refEndB = get<0>(seedingB);
        long long int qryEndB = get<1>(seedingB);

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

        long long int gapNumber = abs(refEndB - refEndA - (qryEndB - qryEndA));

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

    // chaining -> 迭代法 tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chaining(map<string, map<long long int, list<tuple<long long int, int>>>> & seedingChrMap,                                                                           
                                                                                              const string & qryChromosome)
    {
        // 记录chaining信息
        // map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>
        map<float, list<tuple<long long int, long long int, int, string>>> chainingTupListMap;

        // 对总的map进行循环，其中it.first是ref染色体号，it.second是对应的seeding的哈希表
        for (auto it : seedingChrMap)
        {
            string refChromosome = it.first;
            map<long long int, list<tuple<long long int, int>>> seedingMap = it.second;
        
            for (auto iter = seedingMap.begin(); iter != seedingMap.end(); iter++)
            {
                // seeding list1   更新list，记录被选择过的节点
                list<tuple<long long int, int>> seedingList1Tmp = iter->second;

                // 先检查是不是空的list，是的话跳过该节点
                if (seedingList1Tmp.size() == 0)
                {
                    continue;
                }

                for (auto iter1 = seedingList1Tmp.begin(); iter1 != seedingList1Tmp.end(); iter1++)
                {
                    // 记录临时chaining信息
                    // list<tuple(refEnd, qryEnd, length, refChromosome)>
                    list<tuple<long long int, long long int, int, string>> chainingTupListTmp;

                    // this 节点的信息
                    long long int refEnd = get<0>(*iter1);
                    long long int qryEnd = iter->first;
                    int length = get<1>(*iter1);

                    // 记录chaining score
                    float score = length;

                    // 第一个seeding，添加到chainingTupListTmp中
                    chainingTupListTmp.push_back(make_tuple(refEnd, qryEnd, length, refChromosome));

                    // 向前找节点，计算score
                    auto iterTmp1 = iter;
                    // 迭代器后移
                    iterTmp1++;
                    
                    while (iterTmp1 != seedingMap.end())
                    {
                        tuple<long long int, long long int, int, string> seedingA;

                        auto iterTmp = chainingTupListTmp.end();
                        iterTmp--;
                        seedingA = *iterTmp;

                        // seeding list2   更新list，记录被选择过的节点
                        list<tuple<long long int, int>> seedingList2Tmp = iterTmp1->second;
                        // 先检查是不是空的list，是的话跳过该节点
                        if (seedingList2Tmp.size() == 0)
                        {
                            // 迭代器后移
                            iterTmp1++;
                            continue;
                        }

                        // 记录得分最高的索引和分值，以及位置信息
                        float bestScore = 0;
                        long long int bestRefEnd = 0;
                        long long int bestQryEnd = 0;
                        long long bestLength = 0;
                        std::list<tuple<long long int, int>>::iterator bestIter;

                        // 对下一个节点的List进行循环
                        tuple<long long int, long long int, int, string> seedingB;

                        for (auto iter2 = seedingList2Tmp.begin(); iter2 != seedingList2Tmp.end(); iter2++)
                        {
                            // 先检查refEnd在不在上一个的后边，不在的话跳过该seeding
                            if (get<0>(*iter2) < refEnd)
                            {
                                continue;
                            }

                            // 计算两个seeding之间的score
                            seedingB = make_tuple(get<0>(*iter2), 
                                                  iterTmp1->first, 
                                                  get<1>(*iter2),
                                                  qryChromosome);
                            tuple<float, long long int> outTuple= cal_score(seedingA,
                                                                            seedingB, 
                                                                            15);

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
                                                                refChromosome));
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

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> & chainingOutPair)
    {
        // chainingOutPair的score不用变，将map中相交的anchor合并
        // 遍历chainging哈希表
        for (auto iter = get<0>(chainingOutPair).begin(); iter != get<0>(chainingOutPair).end(); iter++)
        {
            for (auto iter1 = iter->second.begin(); iter1 != iter->second.end(); iter1++)
            {
                auto iter1Tmp = iter1;
                iter1Tmp++;
                
                while ((get<0>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<0>(*iter1) && // 第二个的ref起始位置在第一个ref终止位置的后边
                    (get<1>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<1>(*iter1) && // 第二个的qry起始位置在第一个qry终止位置的后边
                    abs(get<1>(*iter1Tmp) - get<1>(*iter1) - (get<0>(*iter1Tmp) - get<0>(*iter1))) == 0 && // 两个anchor之间没有gap
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
    vector<string> run_alignment(const tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> & chainingOutPair, 
                                 const map<string, string> & refFastaMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality)
    {
        for (auto it1 : get<0>(chainingOutPair))
        {
            float score = it1.first;
            list<tuple<long long int, long long int, int, string>> chainingTupList = it1.second;

            // 染色体号
            string refChromosome = get<3>(*chainingTupList.begin());

            // 参考基因组的序列
            string refSeq = "";
            auto iterFind = refFastaMap.find(refChromosome);
            if (iterFind != refFastaMap.end())
            {
                refSeq = iterFind->second;
            }
            else
            {
                cerr << "[" << getTime() << "] [alignment] " 
                    << "Chromosome not present in reference genome: " << refChromosome 
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
                int length = get<2>(*iter);
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
                        string refAliSeq = refSeq.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeq.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

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
                        string refAliSeq = refSeq.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeq.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

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
                        qryAliSeqOut += qrySeq.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

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
                        refAliSeqOut += refSeq.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);

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
            int length = get<2>(*iter);
            int refMinStart = get<0>(*iter) - length + 1;
            int qryMinStart = get<1>(*iter) - length + 1;
            if (qryMinStart > 0)
            {
                string refAliSeq = refSeq.substr(refMinStart - qryMinStart, qryMinStart - 1);
                string qryAliSeq = qrySeq.substr(0, qryMinStart - 1);

                auto alignmentOutTup = alignment::alignment(refAliSeq, qryAliSeq, 9, -3, -2);

                // 加上序列
                refAliSeqOut = get<0>(alignmentOutTup) + refAliSeqOut;
                qryAliSeqOut = get<1>(alignmentOutTup) + qryAliSeqOut;
                midAliSeqOut = get<2>(alignmentOutTup) + midAliSeqOut;
            }
            
            iter = chainingTupList.end();
            iter--;
            int refMaxEnd = get<0>(*iter);
            int qryMaxEnd = get<1>(*iter);
            if (qryMaxEnd < qrySeq.size())
            {
                string refAliSeq = refSeq.substr(refMaxEnd, qrySeq.size() - qryMaxEnd);
                string qryAliSeq = qrySeq.substr(qryMaxEnd, qrySeq.size() - qryMaxEnd);

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
                            refChromosome + "\t" + 
                            to_string(refMinStart) + "\t" + 
                            to_string(get<1>(chainingOutPair)) + "\t" + 
                            cigra + "\t" + 
                            "*\t*\t" + 
                            qrySeq.substr(qryMinStart - 1, qryMaxEnd + 1 - qryMinStart) + "\t" + 
                            qryQuality.substr(qryMinStart - 1, qryMaxEnd + 1 - qryMinStart);
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
        }
        

        vector<string> alignmentVec;

        return alignmentVec;
    }

    // 将比对的结果转为cigra
    string alignment_cigar(const string & refAliSeq, 
                           const string & midAliSeq,
                           const string & qryAliSeq,
                           const long long int & softClipNumLeft,
                           const long long int & softClipNumRight,
                           const long long int & hardClipNumLeft,
                           const long long int & hardClipNumRight)
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
        long long int matchNum = 0;
        long long int delNum = 0;
        long long int insNum = 0;

        // 记录单个的数量
        long long int matchNumTmp = 0;
        long long int delNumTmp = 0;
        long long int insNumTmp = 0;

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


    // 反向互补
    string sequence_reverse(const string & sequence)
    {
        string sequenceRev;
        
        // 反向互补
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
            else if (tmpSeq[i] == 'n')
            {
                sequenceRev += "n";
            }
            else if (tmpSeq[i] == 'N')
            {
                sequenceRev += "N";
            }
            else // 如果序列含有atgcnATGCN之外的序列，则返回空序列
            {
                sequenceRev = "";
                return sequenceRev;
            }
        }

        return sequenceRev;
    }

    // α(A)=0, α(C)=1, α(G)=2, α(T)=3
    // hash函数 (E.g. α(ACATAC) = α(A)*4^5 + α(C)*4^4 + α(A)*4^3 + α(T)*4^2 + α(A)*4^1 + α(C)*4^0 = 305)
    long long int hash_fun(const string & sequence)
    {
        long long int score = 0;
        map<char, int> base2num = {
            {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
        };

        // 计算minimizer值
        int indexTmp = sequence.size() - 1;
        for (auto it : sequence)
        {
            score += base2num[it]*pow(4,indexTmp);
            indexTmp--;
        }

        return score;
    }
}

#endif