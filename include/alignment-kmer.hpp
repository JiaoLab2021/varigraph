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

// kseq.h ���ļ�
KSEQ_INIT(gzFile, gzread)

std::mutex mtx;

// ��tup��qryEnd��������
bool sort_tup_by_sec(tuple<long long int, long long int, int, string> & a, 
                     tuple<long long int, long long int, int, string> & b)
{
    return (get<1>(a) < get<1>(b));
}

namespace alignment
{
    //**********************************************************************//

    // fasta���߳�
    void read_fasta_run(string sequence, 
                        const string & chromosome,
                        const uint32_t& kmerLen, 
                        const int & minimizerK, 
                        const int & minimizerW,
                        unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap);

    // p_read_fastq���߳�
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

    // s_read_fastq���߳�
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

    // Ѱ�ҹ����Ӵ�
    map<string, map<long long int, list<tuple<long long int, int>>>> seeding(const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                                                                             const map<string, map<long long int, vector<long long int>>> & minimizerQryMap, 
                                                                             const int & minimizerK);

    // ����chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<long long int, long long int, int, string> & seedingA,
                                          const tuple<long long int, long long int, int, string> & seedingB, 
                                          const int & minimizerK);
    
    // chaining -> ������ tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chaining(map<string, map<long long int, list<tuple<long long int, int>>>> & seedingChrMap,                                                                           
                                                                                              const string & qryChromosome);

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> & chainingOutPair);

    // ���бȶ�
    vector<string> run_alignment(const tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> & chainingOutPair, 
                                 const map<string, string> & refFastaMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality);

    // ���ȶԵĽ��תΪcigra
    string alignment_cigar(const string & refAliSeq, 
                        const string & midAliSeq,
                        const string & qryAliSeq,
                        const long long int & softClipNumLeft,
                        const long long int & softClipNumRight,
                        const long long int & hardClipNumLeft,
                        const long long int & hardClipNumRight);

    // ���򻥲�
    string sequence_reverse(const string & sequence);

    // hash���� (E.g. ��(ACATAC) = ��(A)*4^5 + ��(C)*4^4 + ��(A)*4^3 + ��(T)*4^2 + ��(A)*4^1 + ��(C)*4^0 = 305)
    long long int hash_fun(const string & sequence);

    //**********************************************************************//





    // ��ȡfasta�ļ�
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

        // ���̳�
        ThreadPool pool(threads);

        // ��ʼ���̳߳�
        pool.init();

        // ������й�ϣ��
        map<string, string> fastaMap;

        // open fasta file
        gzFile gzfp = gzopen(fastaFile.c_str(), "rb");

        // ���ļ�
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
                // ks->name.s ��¼��������
                // ks->seq.s ��¼��������
                string chromosome = ks->name.s;
                string sequence = ks->seq.s;

                pool.submit(read_fasta_run, 
                            sequence, 
                            chromosome, 
                            kmerLen, 
                            minimizerK,
                            minimizerW,
                            ref(minimizerRefMap));

                // ����fasta����
                fastaMap[chromosome] = sequence;

                // ����ַ������ͷ��ڴ�
                chromosome.clear();
                sequence.clear();
                string().swap(chromosome);
                string().swap(sequence);

                // �����������Ƿ񳬹���ֵ�������˵ȴ����Է�������һ���Լ��ص��ڴ���
                while (pool.get_queue() >= threads*10)
                {
                    // ÿ��0.5����һ��
                    sleep(0.5);
                }
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        // �����������Ƿ�ִ���ִ꣬������ر��̳߳أ�����ÿ��0.5s���һ��
        while (pool.get_queue() > 0)
        {
            // ÿ��0.5����һ��
            sleep(0.5);
        }

        // �ر��̳߳�
        pool.shutdown();

        // log        
        cerr << "[" << getTime() << "] [alignment] " 
                 << "Index build succeeded." 
                 << endl;

        return fastaMap;
    }

    // fasta���߳�
    void read_fasta_run(string sequence, 
                        const string & chromosome,
                        const uint32_t& kmerLen, 
                        const int & minimizerK, 
                        const int & minimizerW,
                        unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap)
    {
        // ת��д
        transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

        // ����kmer����
        map<long long int,string> kmerMap = build_kmer_index(sequence, kmerLen);
        
        // ����minimizer����
        build_minimizer_index(kmerMap, 
                              sequence.length(), 
                              chromosome, 
                              minimizerRefMap, 
                              minimizerK, 
                              minimizerW);
        // �ͷ��ڴ�
        kmerMap.clear();
        map<long long int,string>().swap(kmerMap);
    }

    // ��ȡfastq�ļ� ˫�˲���
    void p_read_fastq(const string & fastqFile1, 
                    const string & fastqFile2, 
                    const uint32_t& kmerLen, 
                    const int & minimizerK, 
                    const int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                    const map<string, string> & refFastaMap)
    {
        // ���̳�
        ThreadPool pool(threads);

        // ��ʼ���̳߳�
        pool.init();

        // open fastq1 file
        gzFile gzfp1 = gzopen(fastqFile1.c_str(), "rb");

        // open fastq2 file
        gzFile gzfp2 = gzopen(fastqFile2.c_str(), "rb");

        // ���ļ�
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

                // ks->name.s ��¼��������
                // ks->seq.s ��¼��������
                // ks->qual.s ��¼���ǲ�������
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

                // ����ַ������ͷ��ڴ�
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

                // �����������Ƿ񳬹���ֵ�������˵ȴ����Է�������һ���Լ��ص��ڴ���
                while (pool.get_queue() >= threads*10)
                {
                    // ÿ��0.5����һ��
                    sleep(0.5);
                }
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks1);
            gzclose(gzfp1);
            kseq_destroy(ks2);
            gzclose(gzfp2);
        }

        // �����������Ƿ�ִ���ִ꣬������ر��̳߳أ�����ÿ��0.5s���һ��
        while (pool.get_queue() > 0)
        {
            // ÿ��0.5����һ��
            sleep(0.5);
        }

        // �ر��̳߳�
        pool.shutdown();
    }

    // p_read_fastq���߳�
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
        // ת��д
        transform(sequence1.begin(),sequence1.end(),sequence1.begin(),::toupper);
        transform(sequence2.begin(),sequence2.end(),sequence2.begin(),::toupper);

        // ����kmer����
        map<long long int,string> kmerMap1 = build_kmer_index(sequence1, kmerLen);
        map<long long int,string> kmerMap2 = build_kmer_index(sequence2, kmerLen);
        
        // ��¼minimizer����
        map<string, map<long long int, vector<long long int>>> minimizerQryMap1;
        map<string, map<long long int, vector<long long int>>> minimizerQryMap2;
                          
        // ����minimizer����
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
        // �ͷ��ڴ�
        kmerMap1.clear();
        map<long long int,string>().swap(kmerMap1);
        kmerMap2.clear();
        map<long long int,string>().swap(kmerMap2);

        // seeding
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap1 = seeding(minimizerRefMap, minimizerQryMap1, minimizerK);
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap2 = seeding(minimizerRefMap, minimizerQryMap2, minimizerK);

        // �ͷ��ڴ�
        minimizerQryMap1.clear();
        map<string, map<long long int, vector<long long int>>>().swap(minimizerQryMap1);
        minimizerQryMap2.clear();
        map<string, map<long long int, vector<long long int>>>().swap(minimizerQryMap2);
   
        // chaining and alignment
        if (seedingMap1.size() > 0)
        {
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chainingOutPair1 = chaining(seedingMap1, chromosome);
            
            // �ͷ��ڴ�
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

            // �ͷ��ڴ�
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int>().swap(chainingOutPair1);
        }

        if (seedingMap2.size() > 0)
        {
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chainingOutPair2 = chaining(seedingMap2, chromosome);
            
            // �ͷ��ڴ�
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

            // �ͷ��ڴ�
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int>().swap(chainingOutPair2);
        }
    }

    // ��ȡfastq�ļ�  ���˲���
    void s_read_fastq(const string & fastqFile, 
                    const uint32_t& kmerLen, 
                    const int & minimizerK, 
                    const int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                    const map<string, string> & refFastaMap)
    {
        // ���̳�
        ThreadPool pool(threads);

        // ��ʼ���̳߳�
        pool.init();

        // open fastq1 file
        gzFile gzfp = gzopen(fastqFile.c_str(), "rb");

        // ���ļ�
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

                // ks->name.s ��¼��������
                // ks->seq.s ��¼��������
                // ks->qual.s ��¼���ǲ�������
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

                // ����ַ������ͷ��ڴ�
                chromosome.clear();
                sequence.clear();
                quality.clear();
                string().swap(chromosome);
                string().swap(sequence);
                string().swap(quality);

                // �����������Ƿ񳬹���ֵ�������˵ȴ����Է�������һ���Լ��ص��ڴ���
                while (pool.get_queue() >= threads*10)
                {
                    // ÿ��0.5����һ��
                    sleep(0.5);
                }
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        // �����������Ƿ�ִ���ִ꣬������ر��̳߳أ�����ÿ��0.5s���һ��
        while (pool.get_queue() > 0)
        {
            // ÿ��0.5����һ��
            sleep(0.5);
        }

        // �ر��̳߳�
        pool.shutdown();
    }

    // s_read_fastq���߳�
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

        // ת��д
        transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

        // ����kmer����
        map<long long int,string> kmerMap = build_kmer_index(sequence, kmerLen);
        
        // ��¼minimizer����
        map<string, map<long long int, vector<long long int>>> minimizerQryMap;

        // ����minimizer����
        build_minimizer_index(kmerMap, 
                              sequence.length(), 
                              chromosome, 
                              minimizerQryMap, 
                              minimizerK, 
                              minimizerW);

        // �ͷ��ڴ�
        kmerMap.clear();
        map<long long int,string>().swap(kmerMap);

        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Seeding: " << chromosome
             << endl;

        // seeding
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap = seeding(minimizerRefMap, minimizerQryMap, minimizerK);

        // �ͷ��ڴ�
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
           
            // �ͷ��ڴ�
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

            // �ͷ��ڴ�
            tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int>().swap(chainingOutPair);
        }
    }

    // build kmer index
    map<long long int,string> build_kmer_index(const string & sequence, const uint32_t& kmerLen)
    {
        // kmerMap<start, sequence> ������
        map<long long int,string> kmerMap;

        long long int seqLen = sequence.size();
        long long int kmerNumber = seqLen - kmerLen + 1;

        for (long long int i = 0; i < kmerNumber; i++)
        {
            // �洢kmer����ϣ����
            string kmerSeq = sequence.substr(i, kmerLen);

            // �����1��������0��ʼ�����������Ǵ�1��ʼ
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
                iterTmp++;
                for (int j = 1; j < minimizerW; j++)
                {
                    // ���жϵ�������û��ָ��������ָ��end�Ļ�������ѭ��
                    if (iterTmp == kmerMap.end())
                    {
                        possibleMinimizer.clear(); // possibleMinimizer���
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
                    // ���߳�������
                    std::lock_guard<std::mutex> mtx_locker(mtx);

                    // �������ж�Ӧ��hashֵ
                    long long int sequenceHashNum = hash_fun(possibleMinimizer);
                    // Ⱦɫ�� minimizer���� ��ֹλ��
                    minimizerMap[chromosome][sequenceHashNum].push_back(iter->first + i + minimizerK - 1);
                }
            }
        }

        // ��ͷ
        auto iterTmp = kmerMap.begin();
        string minimizer = iterTmp->second.substr(0, minimizerK);
        if (minimizer.length() > 0)
        {
            // ���߳�������
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // �������ж�Ӧ��hashֵ
            long long int sequenceHashNum = hash_fun(minimizer);
            // Ⱦɫ�� minimizer���� ��ֹλ��
            minimizerMap[chromosome][sequenceHashNum].push_back(minimizerK);
        }

        // ��β
        iterTmp = kmerMap.end();
        iterTmp--;
        minimizer.clear();
        minimizer = iterTmp->second.substr(iterTmp->second.length()-minimizerK, minimizerK);
        if (minimizer.length() > 0)
        {
            // ���߳�������
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // �������ж�Ӧ��hashֵ
            long long int sequenceHashNum = hash_fun(minimizer);
            // Ⱦɫ�� minimizer���� ��ֹλ��
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
                iterTmp++;
                for (int j = 1; j < minimizerW; j++)
                {
                    // ���жϵ�������û��ָ��������ָ��end�Ļ�������ѭ��
                    if (iterTmp == kmerMap.end())
                    {
                        possibleMinimizer.clear(); // possibleMinimizer���
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
                    // ���߳�������
                    std::lock_guard<std::mutex> mtx_locker(mtx);

                    // �������ж�Ӧ��hashֵ
                    long long int sequenceHashNum = hash_fun(possibleMinimizer);
                    // Ⱦɫ�� minimizer���� ��ֹλ��
                    minimizerMap[chromosome][sequenceHashNum].push_back(iter->first + i + minimizerK - 1);
                }
            }
        }

        // ��ͷ
        auto iterTmp = kmerMap.begin();
        string minimizer = iterTmp->second.substr(0, minimizerK);
        if (minimizer.length() > 0)
        {
            // ���߳�������
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // �������ж�Ӧ��hashֵ
            long long int sequenceHashNum = hash_fun(minimizer);
            // Ⱦɫ�� minimizer���� ��ֹλ��
            minimizerMap[chromosome][sequenceHashNum].push_back(minimizerK);
        }

        // ��β
        iterTmp = kmerMap.end();
        iterTmp--;
        minimizer.clear();
        minimizer = iterTmp->second.substr(iterTmp->second.length()-minimizerK, minimizerK);
        if (minimizer.length() > 0)
        {
            // ���߳�������
            std::lock_guard<std::mutex> mtx_locker(mtx);

            // �������ж�Ӧ��hashֵ
            long long int sequenceHashNum = hash_fun(minimizer);
            // Ⱦɫ�� minimizer���� ��ֹλ��
            minimizerMap[chromosome][sequenceHashNum].push_back(sequenceLen);
        }
        
        return 0;
    }

    // Ѱ�ҹ����Ӵ� map<chromosome, map<qryEnd, list<tuple<refEnd, length>>>>
    map<string, map<long long int, list<tuple<long long int, int>>>> seeding(const unordered_map<string, unordered_map<long long int, vector<long long int>>> & minimizerRefMap, 
                                                                             const map<string, map<long long int, vector<long long int>>> & minimizerQryMap, 
                                                                             const int & minimizerK)
    {
        // map<chromosome, map<qryEnd, list<tuple<refEnd, length>>>>
        map<string, map<long long int, list<tuple<long long int, int>>>> seedingMap;

        for (auto it1 : minimizerQryMap)
        {
            // it1.firstΪȾɫ���
            string qryChromosome = it1.first;
            
            // ��ref��minimizerMap����ѭ��������ÿһ��Ⱦɫ���seedNum
            // it2.firstΪref��Ⱦɫ��ţ� it2.secondΪref�� map<sequence, vector<end>>
            for (auto it2 : minimizerRefMap)
            {
                string refChromosome = it2.first;

                // it1.second Ϊ map<sequence, vector<end>>
                for (auto it3 : it1.second)
                {
                    // it3.first Ϊminimizer����
                    // it3.second Ϊminimizer����ֹλ��Vector

                    // ��ref��minimizerMap��������
                    auto iter1 = it2.second.find(it3.first);
                    
                    if (iter1 != it2.second.end())
                    {
                        // ��Ӧ����ref�ж��λ�ã���ֱ𹹽��ṹ��
                        for (auto it4 : iter1->second)
                        {
                            // ��qryEndVec����ѭ��
                            for (auto it5 : it3.second)
                            {
                                // seedingMap��ֵ
                                seedingMap[refChromosome][it5].push_back(make_tuple(it4, minimizerK));
                            }
                        }
                    }
                }
            }
        }

        return seedingMap;
    }


    // ����chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<long long int, long long int, int, string> & seedingA,
                                          const tuple<long long int, long long int, int, string> & seedingB, 
                                          const int & minimizerK)
    {
        // ����gap����
        long long int refEndA = get<0>(seedingA);
        long long int qryEndA = get<1>(seedingA);

        long long int refEndB = get<0>(seedingB);
        long long int qryEndB = get<1>(seedingB);

        // seeding�����еĳ���
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

        // gap�÷�
        float gapScore = 0;

        if (gapNumber != 0)
        {
            // ����gap���֣���(s) = 0.01*|w|*s + 0.5*log2s
            gapScore = -(0.01*minimizerK*gapNumber + 0.5*(log(gapNumber)/log(2)));
        }
        
        int outScore = 0;

        outScore = seqScore + gapScore;

        return make_tuple(outScore, gapNumber);
    }

    // chaining -> ������ tuple<map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<long long int, long long int, int, string>>>, int> chaining(map<string, map<long long int, list<tuple<long long int, int>>>> & seedingChrMap,                                                                           
                                                                                              const string & qryChromosome)
    {
        // ��¼chaining��Ϣ
        // map<score, list<tuple(refEnd, qryEnd, length, refChromosome)>>
        map<float, list<tuple<long long int, long long int, int, string>>> chainingTupListMap;

        // ���ܵ�map����ѭ��������it.first��refȾɫ��ţ�it.second�Ƕ�Ӧ��seeding�Ĺ�ϣ��
        for (auto it : seedingChrMap)
        {
            string refChromosome = it.first;
            map<long long int, list<tuple<long long int, int>>> seedingMap = it.second;
        
            for (auto iter = seedingMap.begin(); iter != seedingMap.end(); iter++)
            {
                // seeding list1   ����list����¼��ѡ����Ľڵ�
                list<tuple<long long int, int>> seedingList1Tmp = iter->second;

                // �ȼ���ǲ��ǿյ�list���ǵĻ������ýڵ�
                if (seedingList1Tmp.size() == 0)
                {
                    continue;
                }

                for (auto iter1 = seedingList1Tmp.begin(); iter1 != seedingList1Tmp.end(); iter1++)
                {
                    // ��¼��ʱchaining��Ϣ
                    // list<tuple(refEnd, qryEnd, length, refChromosome)>
                    list<tuple<long long int, long long int, int, string>> chainingTupListTmp;

                    // this �ڵ����Ϣ
                    long long int refEnd = get<0>(*iter1);
                    long long int qryEnd = iter->first;
                    int length = get<1>(*iter1);

                    // ��¼chaining score
                    float score = length;

                    // ��һ��seeding����ӵ�chainingTupListTmp��
                    chainingTupListTmp.push_back(make_tuple(refEnd, qryEnd, length, refChromosome));

                    // ��ǰ�ҽڵ㣬����score
                    auto iterTmp1 = iter;
                    // ����������
                    iterTmp1++;
                    
                    while (iterTmp1 != seedingMap.end())
                    {
                        tuple<long long int, long long int, int, string> seedingA;

                        auto iterTmp = chainingTupListTmp.end();
                        iterTmp--;
                        seedingA = *iterTmp;

                        // seeding list2   ����list����¼��ѡ����Ľڵ�
                        list<tuple<long long int, int>> seedingList2Tmp = iterTmp1->second;
                        // �ȼ���ǲ��ǿյ�list���ǵĻ������ýڵ�
                        if (seedingList2Tmp.size() == 0)
                        {
                            // ����������
                            iterTmp1++;
                            continue;
                        }

                        // ��¼�÷���ߵ������ͷ�ֵ���Լ�λ����Ϣ
                        float bestScore = 0;
                        long long int bestRefEnd = 0;
                        long long int bestQryEnd = 0;
                        long long bestLength = 0;
                        std::list<tuple<long long int, int>>::iterator bestIter;

                        // ����һ���ڵ��List����ѭ��
                        tuple<long long int, long long int, int, string> seedingB;

                        for (auto iter2 = seedingList2Tmp.begin(); iter2 != seedingList2Tmp.end(); iter2++)
                        {
                            // �ȼ��refEnd�ڲ�����һ���ĺ�ߣ����ڵĻ�������seeding
                            if (get<0>(*iter2) < refEnd)
                            {
                                continue;
                            }

                            // ��������seeding֮���score
                            seedingB = make_tuple(get<0>(*iter2), 
                                                  iterTmp1->first, 
                                                  get<1>(*iter2),
                                                  qryChromosome);
                            tuple<float, long long int> outTuple= cal_score(seedingA,
                                                                            seedingB, 
                                                                            15);

                            // �жϸ�seeding���Ϻ��score�÷��Ƿ������ֵ��gap����ҪС��1000
                            if ((score + get<0>(outTuple) < 0 && get<1>(outTuple) > 1000))
                            {
                                continue;
                            }
                            else
                            {
                                // seeding �÷�
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

                        // ���ж�bestRenEnd�ǲ���0���ǵĻ�����û��ƥ�䣬������seeding
                        if (bestRefEnd == 0)
                        {
                            // ����������
                            iterTmp1++;
                            continue;
                        }
                        
                        score += bestScore;
                        chainingTupListTmp.push_back(make_tuple(bestRefEnd, 
                                                                bestQryEnd, 
                                                                bestLength, 
                                                                refChromosome));
                        // ɾ�����ù��Ľڵ�
                        seedingList2Tmp.erase(bestIter);
                        // ����list��ɾ����ѡ����Ľڵ�
                        seedingMap[iterTmp1->first] = seedingList2Tmp;
                        
                        // ����������
                        iterTmp1++;
                    }

                    // ���chaining�е�seeding�������Լ�score��������ֵ�Ļ�����
                    if (chainingTupListTmp.size() < 3 || score < 40)
                    {
                        continue;
                    }

                    // ���chaining
                    // �ȼ����û��һ���ĵ÷֣��еĻ�score��1
                    while (chainingTupListMap.find(score) != chainingTupListMap.end())
                    {
                        score += 0.01;
                    }
                    chainingTupListMap[score] = chainingTupListTmp;
                }
            }
        }

        // �ȼ���ϣ���ǲ��ǿյ�
        int mapQ = 0;
        if (chainingTupListMap.size() >= 1)
        {
            auto iter = chainingTupListMap.end();
            iter--;
            int anchorsNum = iter->second.size();
            float scoreA = iter->first;
            float scoreB;

            // ���ж��ǲ���ֻ��һ��chaining���ǵĻ�scoreBΪ0
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
        // chainingOutPair��score���ñ䣬��map���ཻ��anchor�ϲ�
        // ����chainging��ϣ��
        for (auto iter = get<0>(chainingOutPair).begin(); iter != get<0>(chainingOutPair).end(); iter++)
        {
            for (auto iter1 = iter->second.begin(); iter1 != iter->second.end(); iter1++)
            {
                auto iter1Tmp = iter1;
                iter1Tmp++;
                
                while ((get<0>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<0>(*iter1) && // �ڶ�����ref��ʼλ���ڵ�һ��ref��ֹλ�õĺ��
                    (get<1>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<1>(*iter1) && // �ڶ�����qry��ʼλ���ڵ�һ��qry��ֹλ�õĺ��
                    abs(get<1>(*iter1Tmp) - get<1>(*iter1) - (get<0>(*iter1Tmp) - get<0>(*iter1))) == 0 && // ����anchor֮��û��gap
                    iter1Tmp != iter->second.end()) // iter1Tmp��ָ���ϣ���end
                {
                    // ���������Ļ������µ�һ��anchor������
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

    // ���бȶ�
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

            // Ⱦɫ���
            string refChromosome = get<3>(*chainingTupList.begin());

            // �ο������������
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

            // �ȶԽ��
            string refAliSeqOut = "";
            string qryAliSeqOut = "";
            string midAliSeqOut = "";

            // ��¼�ȶԵ�����
            long long int aliRefEndTmp = 0;
            long long int aliQryEndTmp = 0;

            // ����chaining��seeding�����бȶ�
            for (auto iter = chainingTupList.begin(); iter != chainingTupList.end(); iter++)
            {
                int length = get<2>(*iter);
                long long int refEnd = get<0>(*iter);
                long long int refStart = refEnd - length + 1;

                long long int qryEnd = get<1>(*iter);
                long long int qryStart = qryEnd - length + 1;

                // ��һ��anchor
                if (iter == chainingTupList.begin())
                {
                    for (size_t i = 0; i < length; i++)
                    {
                        // ����seeding�����У��������У�ƥ��ĵط�����N���
                        refAliSeqOut += "N";
                        qryAliSeqOut += "N";
                        midAliSeqOut += "|";
                    }

                    // ���±ȶԵ�����
                    aliRefEndTmp = refEnd;
                    aliQryEndTmp = qryEnd;
                }
                // ���ǵ�һ��anchorʱ
                else
                {
                    long long int aliRefStart = aliRefEndTmp + 1;
                    long long int aliRefEnd = refStart - 1;

                    long long int aliQryStart = aliQryEndTmp + 1;
                    long long int aliQryEnd = qryStart - 1;

                    // �жϱ�������
                    // snp
                    if ((aliRefEnd - aliRefStart == 1) && (aliQryEnd - aliQryStart == 1))
                    {
                        string refAliSeq = refSeq.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeq.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

                        // ��������
                        refAliSeqOut += refAliSeq;
                        qryAliSeqOut += qryAliSeq;
                        midAliSeqOut += " ";
                        
                        for (size_t i = 0; i < length; i++)
                        {
                            // ����seeding�����У��������У�ƥ��ĵط�����N���
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // ���±ȶԵ�����
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;
                    }
                    // �滻
                    else if((aliRefEnd - aliRefStart > 1) && (aliQryEnd - aliQryStart > 1))
                    {
                        string refAliSeq = refSeq.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeq.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

                        auto alignmentOutTup = alignment::alignment(refAliSeq, qryAliSeq, 9, -3, -2);

                        // ��������
                        refAliSeqOut += get<0>(alignmentOutTup);
                        qryAliSeqOut += get<1>(alignmentOutTup);
                        midAliSeqOut += get<2>(alignmentOutTup);
                        
                        for (size_t i = 0; i < length; i++)
                        {
                            // ����seeding�����У��������У�ƥ��ĵط�����N���
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // ���±ȶԵ�����
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;
                    }
                    // �ж���ʼλ���ǲ��Ǵ�����ֹλ�ã�ref�ϳ���
                    // 1_RefStart:2728 1_RefEnd:4304 1_QryStart:1647 1_QryEnd:3223  2_RefStart:4305 2_RefEnd:4328 2_QryStart:3256 2_QryEnd:3279
                    // INS
                    else if ((aliRefStart - aliRefEnd == 1) && aliQryStart <= aliQryEnd)
                    {
                        // qry��������
                        qryAliSeqOut += qrySeq.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

                        // refֱ�Ӽ�gap
                       for (size_t i = 0; i < (aliQryEnd - aliQryStart + 1); i++)
                        {
                            refAliSeqOut += "-";
                            midAliSeqOut += " ";
                        }
                        
                        for (size_t i = 0; i < length; i++)
                        {
                            // ����seeding�����У��������У�ƥ��ĵط�����N���
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // ���±ȶԵ�����
                        aliRefEndTmp = refEnd;
                        aliQryEndTmp = qryEnd;

                        continue;
                    }
                    // qry�ϵ���ʼλ�ô�����ֹλ��
                    // 1_RefStart:2674 1_RefEnd:2723 1_QryStart:1597 1_QryEnd:1646  2_RefStart:2728 2_RefEnd:4304 2_QryStart:1647 2_QryEnd:3223
                    // DEL
                    else if ((aliQryStart - aliQryEnd == 1) && aliRefStart <= aliRefEnd)
                    {
                        // qryֱ�Ӽ�gap
                       for (size_t i = 0; i < (aliRefEnd - aliRefStart + 1); i++)
                        {
                            qryAliSeqOut += "-";
                            midAliSeqOut += " ";
                        }
                        // ref��������
                        refAliSeqOut += refSeq.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);

                        for (size_t i = 0; i < length; i++)
                        {
                            // ����seeding�����У��������У�ƥ��ĵط�����N���
                            refAliSeqOut += "N";
                            qryAliSeqOut += "N";
                            midAliSeqOut += "|";
                        }

                        // ���±ȶԵ�����
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

            // qry�������бȶ�
            auto iter = chainingTupList.begin();
            int length = get<2>(*iter);
            int refMinStart = get<0>(*iter) - length + 1;
            int qryMinStart = get<1>(*iter) - length + 1;
            if (qryMinStart > 0)
            {
                string refAliSeq = refSeq.substr(refMinStart - qryMinStart, qryMinStart - 1);
                string qryAliSeq = qrySeq.substr(0, qryMinStart - 1);

                auto alignmentOutTup = alignment::alignment(refAliSeq, qryAliSeq, 9, -3, -2);

                // ��������
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

                // ��������
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

            // 1. ����->0   ˫��->1
            // 2. �����ȶԣ������PE���򣬻�����PE������read֮��ıȶԾ���û�����Ե�ƫ�����Ƭ�γ���
            // 4. ����û��mapping���ο�������
            // 8. ���е���һ������û�бȶԵ��ο������ϣ���������������R1,����Ӧ��R2������û�бȶԵ��ο�������
            // 16. ���бȶԵ��ο����еĸ�����
            // 32. ���ж�Ӧ����һ�����бȶԵ��ο����еĸ�����
            // 64. ������R1�����У�read1
            // 128. ������R2�����У�read2
            // 256. ���в�����Ҫ�ıȶ�
            // 512. ��readû��ͨ����������
            // 1024. ������PCR�ظ�����
            // 2048. read���ܴ���Ƕ��
            // ���ȶԽ��תΪFLAG
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

    // ���ȶԵĽ��תΪcigra
    string alignment_cigar(const string & refAliSeq, 
                           const string & midAliSeq,
                           const string & qryAliSeq,
                           const long long int & softClipNumLeft,
                           const long long int & softClipNumRight,
                           const long long int & hardClipNumLeft,
                           const long long int & hardClipNumRight)
    {
        // ��¼ת���Ľ��
        string cigar;

        // // ���ж��ǲ�����ȫƥ�䣬�ǵĻ�ֱ����=
        // if (midAliSeq.find(" ") == string::npos)
        // {
        //     cigar += to_string(midAliSeq.length()) + "=";
        //     return cigar;
        // }
        
        // ���������У���
        if (softClipNumLeft != 0)
        {
            cigar += to_string(softClipNumLeft) + "S";
        }
        else if (hardClipNumLeft != 0)
        {
            cigar += to_string(hardClipNumLeft) + "H";
        }
        
        // ��¼�ܵ�����
        long long int matchNum = 0;
        long long int delNum = 0;
        long long int insNum = 0;

        // ��¼����������
        long long int matchNumTmp = 0;
        long long int delNumTmp = 0;
        long long int insNumTmp = 0;

        // ��¼��һ���ȶԵ�����
        string typeTmp = "";

        // ѭ��midAliSeqת��
        for (size_t i = 0; i < midAliSeq.size(); i++)
        {
            // ��¼��ǰ�ıȶ�����
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

            // �ȶԵ����ͺ���һ���ȶ����Ͳ�һ��
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

            // ��������1
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

            // ���һ���ȶ����
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

        // ���������У���
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


    // ���򻥲�
    string sequence_reverse(const string & sequence)
    {
        string sequenceRev;
        
        // ���򻥲�
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
            else // ������к���atgcnATGCN֮������У��򷵻ؿ�����
            {
                sequenceRev = "";
                return sequenceRev;
            }
        }

        return sequenceRev;
    }

    // ��(A)=0, ��(C)=1, ��(G)=2, ��(T)=3
    // hash���� (E.g. ��(ACATAC) = ��(A)*4^5 + ��(C)*4^4 + ��(A)*4^3 + ��(T)*4^2 + ��(A)*4^1 + ��(C)*4^0 = 305)
    long long int hash_fun(const string & sequence)
    {
        long long int score = 0;
        map<char, int> base2num = {
            {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
        };

        // ����minimizerֵ
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