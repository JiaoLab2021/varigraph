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

// kseq.h ���ļ�
KSEQ_INIT(gzFile, gzread)


bool sort_tup_by_sec(tuple<string, long long int, int> & a, 
                     tuple<string, long long int, int> & b)
{
    return (get<1>(a) < get<1>(b));
}

namespace alignment
{
    //**********************************************************************//

    // fasta���߳�
    void read_fasta_run(const uint64_t chromosomeId, 
                        string sequence, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW,
                        unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap);

    // p_read_fastq���߳�
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

    // s_read_fastq���߳�
    void s_read_fastq_run(const string chromosome, 
                        const string sequence, 
                        const string  quality, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW,
                        const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                        const unordered_map<uint64_t, string> & refFastaMap, 
                        const unordered_map<uint64_t, string> & chromosomeIdMap);

    // Ѱ�ҹ����Ӵ� map<chromosome, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
    map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seeding(const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                                                                                                            const unordered_map<uint64_t, string> & chromosomeIdMap, 
                                                                                                            const unordered_map<uint64_t, list<uint64_t>> & minimizerQryMap, 
                                                                                                            const unsigned int & minimizerK);

    // ����chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<unsigned long long int, unsigned long long int, unsigned int, uint64_t> & seedingA,
                                          const tuple<unsigned long long int, unsigned long long int, unsigned int, uint64_t> & seedingB, 
                                          const unsigned int & minimizerK);
    
    // chaining -> ������ tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chaining(
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> & seedingChrMap, 
        const unsigned int & minimizerK
        );

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair);

    // ���бȶ�
    vector<string> run_alignment(const tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair, 
                                 const unordered_map<uint64_t, string> & refFastaMap, 
                                 const unordered_map<uint64_t, string> & chromosomeIdMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality);

    // ���ȶԵĽ��תΪcigra
    string alignment_cigar(const string & refAliSeq, 
                            const string & midAliSeq,
                            const string & qryAliSeq,
                            const unsigned long long int & softClipNumLeft,
                            const unsigned long long int & softClipNumRight,
                            const unsigned long long int & hardClipNumLeft,
                            const unsigned long long int & hardClipNumRight);

    //**********************************************************************//





    /** ��ȡfasta�ļ� 
     * @param fastaFile  fasta�ļ�
     * @param minimizerK  minimizer����
     * @param minimizerW  minimizer���ڴ�С
     * @param threads  �߳���
     * @param minimizerRefMap  minimizer�����ϣ��
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

        // ���̳�
        ThreadPool pool(threads);

        // ��ʼ���̳߳�
        pool.init();

        // ������й�ϣ��
        unordered_map<uint64_t, string> chromosomeIdMap;
        unordered_map<uint64_t, string> fastaMap;

        // ��ʼ��Ⱦɫ��id
        uint64_t chromosomeId = 0;

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
                // ת��д
                transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

                pool.submit(read_fasta_run, 
                            chromosomeId, 
                            sequence, 
                            minimizerK,
                            minimizerW,
                            ref(minimizerRefMap));

                // ����fasta����
                fastaMap[chromosomeId] = sequence;
                chromosomeIdMap[chromosomeId] = chromosome;

                // ����Ⱦɫ��ID
                chromosomeId++;

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

        return make_tuple(fastaMap, chromosomeIdMap);
    }

    // fasta���߳�
    void read_fasta_run(const uint64_t chromosomeId, 
                        string sequence, 
                        const unsigned int & minimizerK, 
                        const unsigned int & minimizerW,
                        unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap)
    {
        // ����minimizer����
        mm_sketch(chromosomeId, sequence, minimizerK, minimizerW, minimizerRefMap);
    }

    // ��ȡfastq�ļ� ˫�˲���
    void p_read_fastq(const string & fastqFile1, 
                    const string & fastqFile2, 
                    const unsigned int & minimizerK, 
                    const unsigned int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                    const unordered_map<uint64_t, string> & refFastaMap, 
                    const unordered_map<uint64_t, string> & chromosomeIdMap)
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
                // ks->name.s ��¼��������
                // ks->seq.s ��¼��������
                // ks->qual.s ��¼���ǲ�������
                string chromosome1 = ks1->name.s;
                string sequence1 = ks1->seq.s;
                string quality1 = ks1->qual.s;

                string chromosome2 = ks2->name.s;
                string sequence2 = ks2->seq.s;
                string quality2 = ks2->qual.s;

                // ת��д
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
        // ��¼minimizer����
        unordered_map<uint64_t, list<uint64_t>> minimizerQryMap1;
        unordered_map<uint64_t, list<uint64_t>> minimizerQryMap2;
                          
        // ����minimizer����
        mm_sketch(0, sequence1, minimizerK, minimizerW, minimizerQryMap1);
        mm_sketch(0, sequence2, minimizerK, minimizerW, minimizerQryMap2);

        // seeding map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap1 = seeding(minimizerRefMap, chromosomeIdMap, minimizerQryMap1, minimizerK);
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap2 = seeding(minimizerRefMap, chromosomeIdMap, minimizerQryMap2, minimizerK);

        // �ͷ��ڴ�
        minimizerQryMap1.clear();
        unordered_map<uint64_t, list<uint64_t>>().swap(minimizerQryMap1);
        minimizerQryMap2.clear();
        unordered_map<uint64_t, list<uint64_t>>().swap(minimizerQryMap2);
   
        // chaining and alignment
        if (seedingMap1.size() > 0)
        {
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chainingOutPair1 = chaining(seedingMap1, minimizerK);
            
            // �ͷ��ڴ�
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

            // �ͷ��ڴ�
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int>().swap(chainingOutPair1);
        }

        if (seedingMap2.size() > 0)
        {
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chainingOutPair2 = chaining(seedingMap2, minimizerK);
            
            // �ͷ��ڴ�
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

            // �ͷ��ڴ�
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int>().swap(chainingOutPair2);
        }
    }

    // ��ȡfastq�ļ�  ���˲���
    void s_read_fastq(const string & fastqFile, 
                    const unsigned int & minimizerK, 
                    const unsigned int & minimizerW, 
                    const uint32_t & threads,
                    const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
                    const unordered_map<uint64_t, string> & refFastaMap, 
                    const unordered_map<uint64_t, string> & chromosomeIdMap)
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

                // ת��д
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
        
        // ��¼minimizer����
        unordered_map<uint64_t, list<uint64_t>> minimizerQryMap;

        // ����minimizer����
        mm_sketch(0, sequence, minimizerK, minimizerW, minimizerQryMap);

        // log
        cerr << "[" << getTime() << "] [alignment] " 
             << "Seeding: " << chromosome
             << endl;

        // seeding map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap = seeding(minimizerRefMap, chromosomeIdMap, minimizerQryMap, minimizerK);

        // �ͷ��ڴ�
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
           
            // �ͷ��ڴ�
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

            // �ͷ��ڴ�
            tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int>().swap(chainingPair);
        }
        malloc_trim(0);	// 0 is for heap memory
    }


    // Ѱ�ҹ����Ӵ� map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
    map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seeding (
        const unordered_map<uint64_t, list<uint64_t>> & minimizerRefMap, 
        const unordered_map<uint64_t, string> & chromosomeIdMap, 
        const unordered_map<uint64_t, list<uint64_t>> & minimizerQryMap, 
        const unsigned int & minimizerK
        )
    {
        // map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> seedingMap;

        // ����qry��minimizer��ϣ��Ѱ��hit
        for (auto iter1 = minimizerQryMap.begin(); iter1 != minimizerQryMap.end(); iter1++)
        {
            // ��ref�Ĺ�ϣ����Ѱ�Ҷ�Ӧ��seeding
            auto iter2 = minimizerRefMap.find(iter1->first);
            if (iter2 != minimizerRefMap.end())
            {
                // ����qry��list
                for (auto iter3 = iter1->second.begin(); iter3 != iter1->second.end(); iter3++)
                {
                    uint64_t qryChromosomeId = *iter3>>32;
                    unsigned long long int qryEnd = (*iter3 - (qryChromosomeId<<32))>>1;
                    int qryStrand = *iter3 - (qryChromosomeId<<32) - (qryEnd<<1);

                    // ���ж��ǲ����ظ����У��ǵĻ�������list
                    if (iter2->second.size() > 100)
                    {
                        continue;
                    }
                    
                    // ����ref��list
                    for (auto iter4 = iter2->second.begin(); iter4 != iter2->second.end(); iter4++)
                    {
                        uint64_t refChromosomeId = *iter4>>32;
                        unsigned long long int refEnd = (*iter4 - (refChromosomeId<<32))>>1;
                        int refStrand = *iter4 - (refChromosomeId<<32) - (refEnd<<1);

                        // �ж������������Ƿ�һ��
                        int strand = 1;
                        if (qryStrand == refStrand)
                        {
                            strand = 0;
                        }
                        // seedingMap��ֵ
                        seedingMap[refChromosomeId][qryEnd].push_back(make_tuple(refEnd, minimizerK, strand));
                    }
                }
            }
        }

        return seedingMap;
    }


    // ����chaining score tuple<outScore, gapNumber>
    tuple<float, long long int> cal_score(const tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> & seedingA,
                                          const tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> & seedingB, 
                                          const unsigned int & minimizerK)
    {
        // ����gap����
        unsigned long long int refEndA = get<0>(seedingA);
        unsigned long long int qryEndA = get<1>(seedingA);

        unsigned long long int refEndB = get<0>(seedingB);
        unsigned long long int qryEndB = get<1>(seedingB);

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

        long long int gapNumber = refEndB - refEndA - (qryEndB - qryEndA);
        gapNumber = abs(gapNumber);

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

    // chaining -> ������ tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>, mapQ>
    tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> chaining(
        map<uint64_t, map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>>> & seedingChrMap, // map<chromosomeId, map<qryEnd, list<tuple<refEnd, length, strand>>>> 0-F 1-R
        const unsigned int & minimizerK
        )
    {
        // ��¼chaining��Ϣ
        // map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosome)>>
        map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>> chainingTupListMap;

        // ���ܵ�map����ѭ��������iter->first��refȾɫ��ţ�iter->second�Ƕ�Ӧseeding�Ĺ�ϣ��
        for (auto iter = seedingChrMap.begin(); iter != seedingChrMap.end(); iter++)
        {
            uint64_t refChromosomeId = iter->first;
            map<unsigned long long int, list<tuple<unsigned long long int, unsigned int, int>>> seedingMap = iter->second;
        
            for (auto iter1 = seedingMap.begin(); iter1 != seedingMap.end(); iter1++)
            {
                // seeding list1   ����list����¼��ѡ����Ľڵ�
                list<tuple<unsigned long long int, unsigned int, int>> seedingList1Tmp = iter1->second;

                // �ȼ���ǲ��ǿյ�list���ǵĻ������ýڵ�
                if (seedingList1Tmp.size() == 0)
                {
                    continue;
                }

                for (auto iter2 = seedingList1Tmp.begin(); iter2 != seedingList1Tmp.end(); iter2++)
                {
                    // ��¼��ʱchaining��Ϣ
                    // list<tuple(refEnd, qryEnd, length, strand, refChromosomeId)>
                    list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>> chainingTupListTmp;

                    // this �ڵ����Ϣ
                    unsigned long long int refEnd = get<0>(*iter2);
                    unsigned long long int qryEnd = iter1->first;
                    int length = get<1>(*iter2);

                    // ��¼chaining score
                    float score = length;

                    // ��һ��seeding����ӵ�chainingTupListTmp��
                    chainingTupListTmp.push_back(make_tuple(refEnd, qryEnd, length, get<2>(*iter2), refChromosomeId));

                    // ����ҽڵ㣬����score
                    auto iterTmp1 = iter1;
                    // ����������
                    iterTmp1++;
                    
                    while (iterTmp1 != seedingMap.end())
                    {
                        tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> seedingA;

                        auto iterTmp2 = chainingTupListTmp.end();
                        iterTmp2--;
                        seedingA = *iterTmp2;

                        // seeding list2   ����list����¼��ѡ����Ľڵ�
                        list<tuple<unsigned long long int, unsigned int, int>> seedingList2Tmp = iterTmp1->second;
                        // �ȼ���ǲ��ǿյ�list���ǵĻ������ýڵ�
                        if (seedingList2Tmp.size() == 0)
                        {
                            // ����������
                            iterTmp1++;
                            continue;
                        }

                        // ��¼�÷���ߵ������ͷ�ֵ���Լ�λ����Ϣ
                        float bestScore = 0;
                        unsigned long long int bestRefEnd = 0;
                        unsigned long long int bestQryEnd = 0;
                        long long bestLength = 0;
                        int bestStrand = 0;
                        std::list<tuple<unsigned long long int, unsigned int, int>>::iterator bestIter;

                        // ����һ���ڵ��List����ѭ��
                        tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t> seedingB;

                        for (auto iter2 = seedingList2Tmp.begin(); iter2 != seedingList2Tmp.end(); iter2++)
                        {
                            // 1. ���refEnd�ڲ�����һ���ĺ�ߣ����ڵĻ�������seeding
                            // 2. ���strand�ǲ���һ�����򣬲��ǵĻ�������seeding
                            // 3. ���seeding֮���ref����
                            // long long int betweenSeedLenRef = get<1>(*iter2) - get<1>(seedingA);
                            // betweenSeedLenRef = abs(betweenSeedLenRef);
                            if (get<0>(*iter2) < refEnd || get<2>(*iter2) != get<3>(seedingA))
                            {
                                continue;
                            }

                            // ��������seeding֮���score
                            seedingB = make_tuple(get<0>(*iter2), 
                                                  iterTmp1->first, 
                                                  get<1>(*iter2), 
                                                  get<2>(*iter2), 
                                                  refChromosomeId);
                            tuple<float, long long int> outTuple= cal_score(seedingA,
                                                                            seedingB, 
                                                                            minimizerK);

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
                                    bestStrand = get<2>(*iter2);
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
                                                                bestStrand, 
                                                                refChromosomeId));
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

    // anchor_merge tuple<map<score, list<tuple(refEnd, qryEnd, length, strand, refChromosomeId)>>, mapQ>
    void anchor_merge(tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair)
    {
        // chainingOutPair��score���ñ䣬��map���ཻ��anchor�ϲ�
        // ����chainging��ϣ��
        for (auto iter = get<0>(chainingPair).begin(); iter != get<0>(chainingPair).end(); iter++)
        {
            for (auto iter1 = iter->second.begin(); iter1 != iter->second.end(); iter1++)
            {
                auto iter1Tmp = iter1;
                iter1Tmp++;
                
                while ((get<0>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<0>(*iter1) && // �ڶ�����ref��ʼλ���ڵ�һ��ref��ֹλ�õ�ǰ��
                    (get<1>(*iter1Tmp) - get<2>(*iter1Tmp) + 1) <= get<1>(*iter1) && // �ڶ�����qry��ʼλ���ڵ�һ��qry��ֹλ�õ�ǰ��
                    get<1>(*iter1Tmp) - get<1>(*iter1) - (get<0>(*iter1Tmp) - get<0>(*iter1)) == 0 && // ����anchor֮��û��gap
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
    vector<string> run_alignment(const tuple<map<float, list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>>>, int> & chainingPair, 
                                 const unordered_map<uint64_t, string> & refFastaMap, 
                                 const unordered_map<uint64_t, string> & chromosomeIdMap, 
                                 const string & qryChromosome,
                                 const string & qrySeq,
                                 const string & qryQuality)
    {
        // ��¼�ȶ�����
        int aliNum = 0;

        for (auto iter1 = get<0>(chainingPair).rbegin(); iter1 != get<0>(chainingPair).rend(); iter1++)
        {
            // ֻ���scoreǰ����
            if (aliNum > 1)
            {
                break;
            }
            
            float score = iter1->first;
            list<tuple<unsigned long long int, unsigned long long int, unsigned int, int, uint64_t>> chainingTupList = iter1->second;

            // Ⱦɫ���
            uint64_t refChromosomeId = get<4>(*chainingTupList.begin());

            // chainging����
            int strand = get<3>(*chainingTupList.begin());
            string qrySeqConvert = qrySeq;

            // ����strand�ķ���ѡ������
            if (strand == 1)
            {
                qrySeqConvert = sequence_reverse(qrySeq);
            }
            

            // �ο������������
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
                long long int length = get<2>(*iter);
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
                        string refAliSeq = iterFind1->second.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeqConvert.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

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
                        string refAliSeq = iterFind1->second.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);
                        string qryAliSeq = qrySeqConvert.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

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
                        qryAliSeqOut += qrySeqConvert.substr(aliQryStart - 1, aliQryEnd + 1 - aliQryStart);

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
                        refAliSeqOut += iterFind1->second.substr(aliRefStart - 1, aliRefEnd + 1 - aliRefStart);

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

                // ��������
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

                // ��������
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

            // ���±ȶ�����
            aliNum++;
        }
        

        vector<string> alignmentVec;
        malloc_trim(0);	// 0 is for heap memory

        return alignmentVec;
    }

    // ���ȶԵĽ��תΪcigra
    string alignment_cigar(const string & refAliSeq, 
                           const string & midAliSeq,
                           const string & qryAliSeq,
                           const unsigned long long int & softClipNumLeft,
                           const unsigned long long int & softClipNumRight,
                           const unsigned long long int & hardClipNumLeft,
                           const unsigned long long int & hardClipNumRight)
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
        unsigned long long int matchNum = 0;
        unsigned long long int delNum = 0;
        unsigned long long int insNum = 0;

        // ��¼����������
        unsigned long long int matchNumTmp = 0;
        unsigned long long int delNumTmp = 0;
        unsigned long long int insNumTmp = 0;

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
}

#endif