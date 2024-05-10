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
#include <cmath>
#include <set>

#include "save.hpp"
#include "get_time.hpp"
#include "ThreadPool.hpp"
#include "kmer.hpp"
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


struct HS {
    uint8_t h = 0;  // hidden state
    uint8_t c = 0;  // coverage: read
    uint8_t f = 0;  // frequency: genome or graph
};

struct HHS {
    vector<uint16_t> hapVec;  // all haplotype combinations
    vector<HS> HiddenStatesVec;  // All k-mer states
    long double observableScore = 0.0L;  // observable_states function calculates the score
};


namespace GENOTYPE
{
    /**
     * @author zezhen du
     * @date 2023/12/04
     * @version v1.0.2
     * @brief genotype
     * 
     * @param mFastaLenMap           chromosome length
     * @param GraphMap               output of construct_index: map<string, map<uint32_t, nodeSrt>>
     * @param hapMap                 Contains all haplotypes
     * @param vcfHead                the VCF file comment lines
     * @param vcfInfoMap             VCF information
     * @param hapIdxQRmap            map<hapIdx, tuple<quotient, remainder> >
     * @param sampleType             specify the genotype of the sample genome (hom/het)
     * @param samplePloidy           sample ploidy
     * @param hapKmerCoverage        haplotype k-mer coverage
     * @param sampleName             sample name
     * @param outputFileName         output filename
     * @param kmerLen
     * @param haploidNum             the haploid number for genotyping
     * @param chrLenThread           Chromosome granularity
     * @param transitionProType      transition probability type
     * @param svGenotypeBool         structural variation genotyping only
     * @param threads
     * @param debug
     * @param minSupportingReads     the minimum number of supporting reads for a variant
     * 
     * @return 0
    **/
    int genotype(
        const unordered_map<string, uint32_t> & mFastaLenMap, 
        map<string, map<uint32_t, nodeSrt> > & GraphMap, 
        const map<uint16_t, string> & hapMap, 
        const string& vcfHead, 
        map<string, map<uint32_t, vector<string> > > & vcfInfoMap, 
        const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap, 
        const string& sampleType, 
        const uint32_t& samplePloidy, 
        const float& hapKmerCoverage, 
        const string& sampleName, 
        const string & outputFileName, 
        const uint32_t kmerLen, 
        uint32_t haploidNum, 
        uint32_t chrLenThread, 
        const string& transitionProType, 
        const bool& svGenotypeBool, 
        uint32_t threads, 
        const bool & debug, 
        const float & minSupportingReads
    );


    /**
     * @author zezhen du
     * @date 2023/12/04
     * @version v1.0.1
     * @brief forward/backward run
     * 
     * @param startNodeIter        node iterator
     * @param hapMap               haplotype information
     * @param sampleType           specify the genotype of the sample genome (hom/het)
     * @param samplePloidy         sample ploidy
     * @param hapKmerCoverage      haplotype k-mer coverage
     * @param hapIdxQRmap          map<hapIdx, tuple<quotient, remainder> >
     * @param threadStart          iterator start position
     * @param threadEnd            iterator end position
     * @param kmerLen
     * @param haploidNum           the haploid number for genotyping
     * @param transitionProType    transition probability type
     * @param svGenotypeBool       structural variation genotyping only
     * @param vcfInfoMap           VCF information, map<chromosome, map<start, vector> >
     * 
     * @return 0
    **/
    int for_bac_post_run(
        map<string, map<uint32_t, nodeSrt> >::iterator startNodeIter, 
        const map<uint16_t, string> & hapMap, 
        const string& sampleType, 
        const uint32_t& samplePloidy, 
        const float& hapKmerCoverage,
        const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap, 
        uint32_t threadStart, 
        uint32_t threadEnd, 
        const uint32_t& kmerLen, 
        const uint32_t& haploidNum, 
        const string& transitionProType, 
        const bool& svGenotypeBool, 
        const map<string, map<uint32_t, vector<string> > >& vcfInfoMap
    );


    /**
     * @author zezhen du
     * @date 2023/12/13
     * @version v1.0
     * @brief haplotype selection
     * 
     * @param chromosome            chromosome
     * @param newStartNodeIterL     thread left iterator
     * @param newStartNodeIterR     thread right iterator
     * @param hapMap                haplotype information
     * @param haploidNum            the haploid number for genotyping
     * @param topHapVec             Haplotype index for final screening in this block
     * @param hapIdxScoreMap        Likelihood of haplotype occurrence: map<hapIdx, possibility>
     * 
     * @return void
    **/
    void haplotype_selection(
        const string& chromosome, 
        const map<uint32_t, nodeSrt>::iterator& newStartNodeIterL, 
        const map<uint32_t, nodeSrt>::iterator& newStartNodeIterR, 
        const map<uint16_t, string>& hapMap, 
        const uint32_t haploidNum, 
        vector<uint16_t>& topHapVec, 
        unordered_map<uint16_t, double>& hapIdxScoreMap
    );


    /**
     * @author zezhen du
     * @date 2023/12/04
     * @version v1.0
     * @brief Hidden states
     * 
     * @param sampleType          specify the genotype of the sample genome (hom/het)
     * @param samplePloidy        sample ploidy
     * @param chromosome          chromosome
     * @param nodeStart           Node start position
     * @param startNodeIter       node iterator, map<string, map<uint32_t, nodeSrt> >
     * @param node                node information, output by construct
     * @param topHapVec           The haplotype finally screened out by the block
     * @param hapIdxQRmap         map<hapIdx, tuple<quotient, remainder> >
     * @param lower               95% confidence interval for average coverage
     * @param upper               95% confidence interval for average coverage
     * @param kmerLen
     * @param filter              Whether to filter GraphKmerHashHapStrMapIterVec
     * 
     * @return HHSStrVec          vector<HHS>
    **/
    vector<HHS> hidden_states(
        const string& sampleType, 
        const uint32_t& samplePloidy, 
        const string& chromosome, 
        const uint32_t& nodeStart, 
        map<string, map<uint32_t, nodeSrt> >::iterator startNodeIter, 
        nodeSrt& node, 
        const vector<uint16_t>& topHapVec, 
        const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap, 
        const double& lower, 
        const double& upper, 
        const uint32_t& kmerLen, 
        bool filter
    );


    /**
     * @author zezhen du
     * @brief Gets the combination of all haplotypes
     * 
     * @param hapVec        Vector of haplotypes     
     * @param sampleType    specify the genotype of the sample genome (hom/het)
     * @param samplePloidy  sample ploidy (2-8) [2]
     * 
     * @return ComHapVec
    **/
    vector<vector<uint16_t> > increment_vector(const vector<uint16_t>& hapVec, const string& sampleType, const uint32_t& samplePloidy);


    /**
     * @author zezhen du
     * @brief Calculate the upper and lower limits of the Poisson distribution
     * 
     * @param lambda     
     * @param lower      lower limit
     * @param upper      upper limit
     * 
     * @return void
    **/
    void calculate_poisson_interval(const double& lambda, double& lower, double& upper);


    /**
     * @author zezhen du
     * @date 2023/07/29
     * @version v1.0
     * @brief Transition probabilities
     * 
     * @param nodeDistence   distance between nodes
     * @param hapMap         Haplotype information, output by construct
     * 
     * @return pair<long double, long double>   recombProb, noRecombProb
    **/
    pair<long double, long double> transition_probabilities(const uint32_t & nodeDistence, const map<uint16_t, string> & hapMap);


    /**
     * @author zezhen du
     * @date 2023/09/06
     * @version v1.0
     * @brief Observable states
     * 
     * @param aveKmerCoverage      Average kmer coverage
     * @param HHSStrVec            All hidden states of node, output by hidden_states
     * @param node                 All information about node
     * 
     * @return void
    **/
    void observable_states(
        const float & aveKmerCoverage, 
        vector<HHS>& HHSStrVec, 
        nodeSrt & node
    );
    

    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
     * @brief poisson
     * 
     * @param mean         average k-mer coverage of nodes
     * @param value        k-mer coverage
     * 
     * @return poisson
    **/
    long double poisson(
        long double mean, 
        uint8_t value
    );

    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief get_error_param
     * 
     * @param aveKmerCoverage         Average coverage of k-mer
     * 
     * @return Correct rate
	**/
    double get_error_param(double aveKmerCoverage);

    /**
     * @author zezhen du
     * @date 2022/12/23
     * @version v1.0
	 * @brief geometric, the return state when the hidden state is 0
     * 
     * @param p         Correct rate, output by get_error_param
     * @param value     coverage
     * 
     * @return geometric
	**/
    long double geometric(long double p, uint8_t value);

    // Prior probability distribution function (normal distribution)
    long double prior(long double p);

    // Likelihood function
    long double likelihood(long double p, uint8_t value);


    /**
     * @author zezhen du
     * @date 2023/12/04
     * @brief find_most_likely_depth
     * 
     * @param h                   hidden state
     * @param c                   k-mer coverage
     * @param f                   k-mer frequency
     * @param aveKmerCoverage     Average depth
     * @param lower               95% confidence interval for average coverage
     * @param upper               95% confidence interval for average coverage
     * 
     * @return void
    **/
    void find_most_likely_depth(
        const uint8_t& h, 
        uint8_t& c, 
        const uint8_t& f, 
        const float& aveKmerCoverage, 
        const double& lower, 
        const double& upper
    );


    /**
     * @author zezhen du
     * @date 2023/11/29
     * @version v1.0
     * @brief forward algorithm
     * 
     * @param preHMMScoreVec        Alpha of previous state
     * @param hapIdxScoreMap        Likelihood of haplotype occurrence: map<hapIdx, possibility>
     * @param recombProb            probability of recombination
     * @param noRecombProb          Probability of non-recombination
     * @param HHSStrVec             observation matrix
     * @param HMMScoreVec           Alpha of this state
     * 
     * @return void
    **/
    void forward(
        const vector<HMMScore>& preHMMScoreVec, 
        const unordered_map<uint16_t, double>& hapIdxScoreMap, 
        const long double& recombProb, 
        const long double& noRecombProb, 
        const vector<HHS>& HHSStrVec, 
        vector<HMMScore>& HMMScoreVec
    );


    /**
     * @author zezhen du
     * @date 2023/11/29
     * @version v1.0
     * @brief backward algorithm
     * 
     * @param preHMMScoreVec        Beta of previous state
     * @param hapIdxScoreMap        Likelihood of haplotype occurrence: map<hapIdx, possibility>
     * @param recombProb            probability of recombination
     * @param noRecombProb          Probability of non-recombination
     * @param HHSStrVec             observation matrix
     * @param HMMScoreVec           Beta of this state
     * 
     * @return void
    **/
    void backward(
        const vector<HMMScore>& preHMMScoreVec, 
        const unordered_map<uint16_t, double>& hapIdxScoreMap, 
        const long double& recombProb, 
        const long double& noRecombProb, 
        const vector<HHS>& HHSStrVec, 
        vector<HMMScore>& HMMScoreVec
    );


    /**
     * @author zezhen du
     * @date 2023/07/27
     * @version v1.0.1
     * @brief Posterior probability
     * 
     * @param nodeStart         Starting position of the node
     * @param nodePtr           Pointer to the node
     * @param topHapVec         Haplotypes screened by Dirichlet distribution
     * 
     * @return 0
    **/
    int posterior(
        uint32_t nodeStart, 
        nodeSrt* nodePtr, 
        vector<uint16_t>& topHapVec
    );


    /**
     * @author zezhen du
     * @date 2023/12/08
     * @version v1.0.1
     * @brief Number of node-unique k-mers, UK
     * 
     * @param GraphKmerHashHapStrMapIterVec      Iterator pointing to mGraphKmerHashHapStrMap, vector<iter>
     * 
     * @return uniqueKmerNum
    **/
    uint8_t get_UK(
        const vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator>& GraphKmerHashHapStrMapIterVec
    );


    /**
     * @author zezhen du
     * @date 2023/07/27
     * @version v1.0.1
     * @brief calculate Phred Scaled
     * 
     * @param value       Posterior probability
     * 
     * @return Phred Scaled
    **/
    double cal_phred_scaled(long double value);


    /**
     * @author zezhen du
     * @date 2023/07/12
     * @version v1.0
     * @brief Save the result
     * 
     * @param GraphMap            output of construct_index£¬map<string, map<uint32_t, nodeSrt>>
     * @param vcfHead             the VCF file comment lines
     * @param vcfInfoMap          vcf information for output
     * @param sampleName          sample name
     * @param outputFileName      Output file information
     * @param minSupportingReads  the minimum number of supporting reads for a variant
     * 
     * @return 0
    **/
    int save(
        map<string, map<uint32_t, nodeSrt> > & GraphMap, 
        string vcfHead, 
        map<string, map<uint32_t, vector<string> > > & vcfInfoMap, 
        const string& sampleName, 
        const string & outputFileName, 
        const float& minSupportingReads
    );
}

#endif