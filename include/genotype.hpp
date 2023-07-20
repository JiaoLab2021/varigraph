#ifndef GENOTYPE_HPP
#define GENOTYPE_HPP
#include <map>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <future>
#include "tuple"
#include <malloc.h>
#include <iomanip>
#include <bitset>

#include "save.hpp"
#include "get_time.hpp"
#include "ThreadPool.hpp"
#include "construct_index.hpp"
#include "haplotype_select.hpp"


// If the macro is not defined, define it
#ifndef DIVIDE_BY_8
    // Get the macro definition of the quotient divided by 8
    #define DIVIDE_BY_8(x) ((x) >> 3)
#endif

// If the macro is not defined, define it
#ifndef GET_LOW_3_BITS
    // Extract the macro definition of the lower 3 bits
    #define GET_LOW_3_BITS(x) ((x) & 0x7)
#endif


using namespace std;

// global variable
extern bool debugGenotype;

namespace GENOTYPE
{
    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief calculate Node k-mers Coverage
     * 
     * @param GraphKmerCovFreMap  map<kmerHash, kmerCovFre>
     * @param kmerHashHapVec      vector<kmerHashHap>   
     * @param nodeAveKmerCoverage the average k-mers coverage of node
     * @param kmerCoverageVec     the k-mers vector of node
     * 
     * @return void
    **/
    void cal_node_cov(
        const std::unordered_map<uint64_t, kmerCovFre> &GraphKmerCovFreMap, 
        const vector<kmerHashHap> & kmerHashHapVec, 
        float & nodeAveKmerCoverage,
        vector<uint8_t> & kmerCoverageVec
    );


    /**
     * @author zezhen du
     * @date 2023/07/17
     * @version v1.0.1
     * @brief forward/backward run
     * 
     * @param GraphKmerCovFreMap  map<kmerHash, kmerCovFre>
     * @param startNodeIter       node iterator
     * @param hapMap              haplotype information
     * @param threadStart         iterator start position
     * @param threadEnd           iterator end position
     * 
     * @return 0
    **/
    int for_bac_post_run(
        const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap, 
        map<string, map<uint32_t, nodeSrt> >::iterator startNodeIter, 
        const map<uint16_t, string> & hapMap, 
        uint32_t threadStart, 
        uint32_t threadEnd
    );


    /**
     * @author zezhen du
     * @date 2023/07/20
     * @version v1.0
     * @brief haplotype selection
     * 
     * @param GraphKmerCovFreMap    Coverage of all k-mers in the graph
     * @param chromosome            chromosome
     * @param newStartNodeIterL     thread left iterator
     * @param newStartNodeIterR     thread right iterator
     * @param hapMap                haplotype information
     * @param topHapVec             Haplotype index for final screening in this block
     * 
     * @return void
    **/
    void haplotype_selection(
        const std::unordered_map<uint64_t, kmerCovFre> &GraphKmerCovFreMap, 
        const string& chromosome, 
        const map<uint32_t, nodeSrt>::iterator& newStartNodeIterL, 
        const map<uint32_t, nodeSrt>::iterator& newStartNodeIterR, 
        const map<uint16_t, string> & hapMap, 
        vector<uint16_t>& topHapVec
    );


    /**
     * @author zezhen du
     * @date 2023/7/17
     * @version v1.0.1
     * @brief genotype
     * 
     * @param GraphKmerCovFreMap     map<kmerHash, kmerCovFre>
     * @param GraphMap               output of construct_index: map<string, map<uint32_t, nodeSrt>>
     * @param hapMap                 Contains all haplotypes
     * @param vcfHead                the VCF file comment lines
     * @param vcfInfoMap             VCF information
     * @param outputFileName         output filename
     * @param threads
     * @param debug
     * 
     * @return 0
    **/
    int genotype(
        const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap, 
        map<string, map<uint32_t, nodeSrt> > & GraphMap, 
        const map<uint16_t, string> & hapMap, 
        const string& vcfHead, 
        map<string, map<uint32_t, string> > & vcfInfoMap, 
        const string & outputFileName, 
        const uint32_t & threads, 
        const bool & debug
    );

    /**
     * @author zezhen du
     * @date 2023/07/20
     * @version v1.0
	 * @brief Hidden states
     * 
     * @param node         node information, output by construct
     * @param topHapVec    The haplotype finally screened out by the block
     * 
     * @return hiddenStatesMap   map<uint8_t, map<uint8_t, vector<uint8_t> > >
	**/
    map<uint16_t, map<uint16_t, vector<uint8_t> > > hidden_states(
        const nodeSrt & node, 
        const vector<uint16_t>& topHapVec
    );


    /**
     * @author zezhen du
     * @date 2022/12/25
     * @version v1.0
	 * @brief Transition probabilities
     * 
     * @param nodeDistence   node之间的距离
     * @param hapMap         单倍型信息，output by construct
     * 
     * @return pair<long double, long double>   recombProb, noRecombProb
	**/
    pair<long double, long double> transition_probabilities(const uint32_t & nodeDistence, const map<uint16_t, string> & hapMap);


    /**
     * @author zezhen du
     * @date 2022/12/25
     * @version v1.0
	 * @brief Observable states
     * 
     * @param aveKmerCoverage      平均的kmer覆盖度
     * @param hiddenStatesMap      node的所有隐藏状态，output by hidden_states
     * @param kmerCoverageVec      node的所有kmer覆盖度
     * @param node                 node所有信息
     * 
     * @return observableStatesMap   map<uint16_t, map<uint16_t, long double> >, map<genotype1, map<genotype2, result>>
	**/
    map<uint16_t, map<uint16_t, long double> > observable_states(
        const float & aveKmerCoverage, 
        const map<uint16_t, map<uint16_t, vector<uint8_t> > > & hiddenStatesMap, 
        const vector<uint8_t> & kmerCoverageVec, 
        nodeSrt & node
    );

    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief poisson
     * 
     * @param mean         kmer的覆盖度
     * @param value        kmer在隐藏状态矩阵中的数值
     * 
     * @return poisson
	**/
    long double poisson(
        long double mean, 
        unsigned int value
    );

    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief get_error_param
     * 
     * @param aveKmerCoverage         kmer的平均覆盖度
     * 
     * @return Correct rate
	**/
    double get_error_param(double aveKmerCoverage);

    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief geometric，隐藏状态为0时的返回状态
     * 
     * @param p         Correct rate, output by get_error_param
     * @param value     隐藏状态中的数值
     * 
     * @return geometric
	**/
    long double geometric(
        long double p, 
        unsigned int value
    );


    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief 向前算法
     * 
     * @param preHMMScoreVec        Alpha of previous state
     * @param recombProb            重组的概率
     * @param noRecombProb          非重组的概率
     * @param observableStatesMap   观察矩阵
     * @param HMMScoreVec           Alpha of this state
     * 
     * @return void
	**/
    void forward(
        const vector<HMMScore> & preHMMScoreVec, 
        const long double & recombProb, 
        const long double & noRecombProb, 
        const map<uint16_t, map<uint16_t, long double> > & observableStatesMap, 
        vector<HMMScore>& HMMScoreVec
    );


    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief 向后算法
     * 
     * @param preHMMScoreVec        Beta of previous state
     * @param recombProb            重组的概率
     * @param noRecombProb          非重组的概率
     * @param observableStatesMap   观察矩阵
     * @param HMMScoreVec           Beta of this state
     * 
     * @return void
	**/
    void backward(
        const vector<HMMScore> & preHMMScoreVec, 
        const long double & recombProb, 
        const long double & noRecombProb, 
        const map<uint16_t, map<uint16_t, long double> > & observableStatesMap, 
        vector<HMMScore>& HMMScoreVec
    );


   /**
     * @author zezhen du
     * @date 2023/07/04
     * @version v1.0.1
     * @brief 后验概率
     * 
     * @param nodeStart         Starting position of the node
     * @param nodePtr           Pointer to the node
     * 
     * @return 0
    **/
    int posterior(
        uint32_t nodeStart, 
        nodeSrt* nodePtr
    );


    /**
     * @author zezhen du
     * @date 2023/07/12
     * @version v1.0
     * @brief 保存结果
     * 
     * @param GraphMap            construct_index输出结果，map<string, map<uint32_t, nodeSrt>>
     * @param vcfHead             the VCF file comment lines
     * @param vcfInfoMap          vcf信息，用于输出
     * @param outputFileName      输出文件信息
     * 
     * @return 0
    **/
    int save(
        map<string, map<uint32_t, nodeSrt> > & GraphMap, 
        string vcfHead, 
        map<string, map<uint32_t, string> > & vcfInfoMap, 
        const string & outputFileName
    );
}

#endif