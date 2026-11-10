// g++ -c src/genotype.cpp -std=c++17 -O3 -march=native

#include "../include/genotype.hpp"
#include "../include/construct_index.hpp"

#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

// global variable
bool debugGenotype = false;
bool polyMixGenotype = false;
double polyMixSwitchRate = 1e-4;   // --switch-rate: per-copy switch rate per bp (poly-mix)
double polyMixMismatch = 0.05;     // --mismatch: per-copy allele mismatch vs. copied haplotype (poly-mix)

namespace {
// All size-k multisets (combinations with replacement) of the sorted haplotype list, each sorted ascending.
void multiset_rec(const vector<uint16_t>& haps, uint32_t k, size_t start, vector<uint16_t>& cur, vector<vector<uint16_t> >& out) {
    if (cur.size() == k) { out.push_back(cur); return; }
    for (size_t i = start; i < haps.size(); ++i) {
        cur.push_back(haps[i]);
        multiset_rec(haps, k, i, cur, out);
        cur.pop_back();
    }
}

vector<vector<uint16_t> > multiset_states(vector<uint16_t> haps, uint32_t k, bool homOnly) {
    std::sort(haps.begin(), haps.end());
    haps.erase(std::unique(haps.begin(), haps.end()), haps.end());
    vector<vector<uint16_t> > out;
    if (homOnly) {
        for (auto h : haps) out.push_back(vector<uint16_t>(k, h));
        return out;
    }
    vector<uint16_t> cur;
    multiset_rec(haps, k, 0, cur, out);
    return out;
}

// One-switch neighbourhood: for each state, the indices of its distinct (k-1)-sub-multisets.
// alpha'(s) = E(s) [noR^k a(s) + noR^(k-1) R sum_{x distinct in s} (A(s-x) - a(s))],  A(m) = sum over t containing m of a(t).
// This is exact for transitions that keep k-1 copies; terms with >=2 switches (R^2) are dropped.
struct MixNeighbour {
    vector<vector<uint16_t> > states;          // cached state list (signature)
    vector<vector<int32_t> > subIdx;           // per state: sub-multiset indices
    size_t nSub = 0;
};

const MixNeighbour& mix_neighbour(const vector<const vector<uint16_t>*>& stateVecs) {
    thread_local MixNeighbour cache;
    bool same = cache.states.size() == stateVecs.size();
    for (size_t i = 0; same && i < stateVecs.size(); ++i) same = (cache.states[i] == *stateVecs[i]);
    if (same) return cache;
    cache.states.clear(); cache.subIdx.clear();
    std::map<vector<uint16_t>, int32_t> subMap;
    for (const auto* sv : stateVecs) {
        cache.states.push_back(*sv);
        vector<int32_t> idx;
        for (size_t j = 0; j < sv->size(); ++j) {
            if (j > 0 && (*sv)[j] == (*sv)[j - 1]) continue;   // distinct elements only
            vector<uint16_t> sub(*sv); sub.erase(sub.begin() + j);
            auto it = subMap.emplace(sub, (int32_t)subMap.size()).first;
            idx.push_back(it->second);
        }
        cache.subIdx.push_back(std::move(idx));
    }
    cache.nSub = subMap.size();
    return cache;
}

// Shared forward/backward step. prevScore(i) is alpha (forward) or beta (backward) of the neighbouring node,
// in the same state order. Returns false if state lists differ (caller falls back to the original O(S^2) loop).
bool mix_step(
    const vector<HMMScore>& pre, bool useA,
    long double recombProb, long double noRecombProb,
    const vector<HHS>& HHSStrVec, vector<long double>& outVec
) {
    vector<const vector<uint16_t>*> cur; vector<long double> emis;
    for (const auto& h : HHSStrVec) {
        if (h.hapVec.empty() || h.HiddenStatesVec.empty()) continue;
        cur.push_back(&h.hapVec); emis.push_back(h.observableScore);
    }
    if (cur.size() != pre.size() || cur.empty()) return false;
    if (pre.front().hapVec != *cur.front() || pre.back().hapVec != *cur.back()) return false;
    const MixNeighbour& nb = mix_neighbour(cur);
    const int k = (int)cur.front()->size();
    vector<long double> A(nb.nSub, 0.0L);
    for (size_t i = 0; i < cur.size(); ++i) {
        long double p = useA ? pre[i].a : pre[i].b;
        for (int32_t m : nb.subIdx[i]) A[m] += p;
    }
    const long double stay = std::pow(noRecombProb, k);
    const long double one = std::pow(noRecombProb, k - 1) * recombProb;
    outVec.assign(cur.size(), 0.0L);
    long double sum = 0.0L;
    for (size_t i = 0; i < cur.size(); ++i) {
        long double p = useA ? pre[i].a : pre[i].b;
        long double sw = 0.0L;
        for (int32_t m : nb.subIdx[i]) sw += A[m] - p;
        outVec[i] = emis[i] * (stay * p + one * sw);
        sum += outVec[i];
    }
    if (sum > 0.0L) { for (auto& v : outVec) v /= sum; }
    else { for (auto& v : outVec) v = 1.0L / (long double)outVec.size(); }
    return true;
}

// Biallelic node under --poly-mix: dosage likelihoods L[d] (d alt copies), P(d | a alt copies in state, theta),
// and the state emission M[a] = sum_d P(d|a) L[d]. Filled in observable_states, used again in posterior
// (forward, backward and posterior of one block run on the same thread).
struct MixNode {
    int P = 0;
    vector<long double> L, M;
    vector<vector<long double> > PdA;
};
thread_local std::unordered_map<const nodeSrt*, MixNode> mixNodeCache;

vector<vector<long double> > dosage_given_alt_copies(int P, double theta) {
    vector<vector<long double> > out(P + 1, vector<long double>(P + 1, 0.0L));
    for (int a = 0; a <= P; ++a) {
        vector<long double> dist(1, 1.0L);
        for (int j = 0; j < P; ++j) {
            long double p = (j < a) ? (1.0L - theta) : theta;   // copy carries alt
            vector<long double> nd(dist.size() + 1, 0.0L);
            for (size_t d = 0; d < dist.size(); ++d) { nd[d] += dist[d] * (1.0L - p); nd[d + 1] += dist[d] * p; }
            dist.swap(nd);
        }
        out[a] = dist;
    }
    return out;
}
}  // namespace

std::mutex mtxG;


/**
 * @author zezhen du
 * @date 2023/12/04
 * @version v1.0.1
 * @brief genotype
 * 
 * @param mFastaLenMap           chromosome length
 * @param GraphMap               output of construct_index: map<string, map<uint32_t, nodeSrt>>
 * @param hapMap                 Contains all haplotypes
 * @param vcfHead                the VCF file comment lines
 * @param vcfInfoMap             VCF information
 * @param hapIdxQRmap            map<hapIdx, tuple<quotient, remainder> >
 * @param sampleType             specify the genotype of the reference genome (hom/het)
 * @param samplePloidy           sample ploidy
 * @param hapKmerCoverage        haplotype k-mer coverage
 * @param sampleName             sample name
 * @param kmerLen
 * @param haploidNum             the haploid number for genotyping
 * @param chrLenThread           Chromosome granularity
 * @param transitionProType      transition probability type
 * @param svGenotypeBool         structural variation genotyping only
 * @param threads
 * @param debug
 * @param minSupportingGQ        minimum site quality (GQ) value for genotype
 * 
 * @return 0
**/
int GENOTYPE::genotype(
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
    const uint32_t kmerLen, 
    uint32_t haploidNum, 
    uint32_t chrLenThread, 
    const string& transitionProType, 
    const bool& svGenotypeBool, 
    uint32_t threads, 
    const bool & debug, 
    const float & minSupportingGQ
) {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Genotyping ...\n";  // print log

    // debug
    debugGenotype = debug;
    if (debugGenotype) threads = 1;

    // Thread Pool
    ThreadPool pool(threads);

    // Save the results of multiple threads
    vector<future<int> > futureOutVec;

    // Initialize the thread pool
    pool.init();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Applying forward and backward algorithm ...\n";

    haploidNum = min(haploidNum, static_cast<uint32_t>(hapMap.size()));  // We need more haplotypes for genotyping than are present in the VCF file, so we will use the number of haplotypes from the VCF file instead.

    for (auto iter1 = GraphMap.begin(); iter1 != GraphMap.end(); iter1++) {  // map<string, map<uint32_t, nodeSrt> >: map<chr, map<nodeStart, nodeSrt> >
        // The number of nodes of this chromosome
        uint64_t chrNodeNum = iter1->second.size();

        // chromosome length
        string chromosome = iter1->first;
        uint32_t chrLen;

       auto findIter = mFastaLenMap.find(chromosome);  // chromosome length

        // Check if the chromosome is present in the reference genome
        if (findIter == mFastaLenMap.end()) {  // if not
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << chromosome << "' does not exist in the reference genome." << endl;
            exit(1);
        } else {  // yes
            chrLen = findIter->second;  // chromosome length
        }

        uint32_t chrLenThreadTmp = min(chrLenThread, chrLen);  // Prevent chromosome granularity from being larger than chromosome length

        uint32_t stepsNumber = ceil(double(chrLen)/chrLenThreadTmp);  // Step count

        // End position of thread
        uint32_t threadEndTmp = 0;

        // Submit tasks in a loop
        for (uint32_t i = 0; i < stepsNumber; i++) {
            uint32_t stepEnd = (i + 1) * chrLenThreadTmp;

            // End position of previous window
            uint32_t threadStart = threadEndTmp;
            
            if (threadStart >= chrNodeNum) break;  // Reach the end of the node of the chromosome and jump out of the loop

            auto threadEndIterTmp = std::next(iter1->second.begin(), threadStart);
            for (map<uint32_t, nodeSrt>::iterator iter2 = threadEndIterTmp; iter2 != iter1->second.end(); iter2++) {  // map<nodeStart, nodeSrt>
                if (iter2->first > stepEnd) break;
                threadEndTmp++;
            }
            uint32_t threadEnd = threadEndTmp;

            // Submit and save results in multiple threads
            futureOutVec.push_back(
                pool.submit(
                    for_bac_post_run, 
                    iter1, 
                    ref(hapMap), 
                    ref(sampleType), 
                    ref(samplePloidy), 
                    ref(hapKmerCoverage), 
                    ref(hapIdxQRmap), 
                    threadStart, 
                    threadEnd, 
                    ref(kmerLen), 
                    ref(haploidNum), 
                    ref(transitionProType), 
                    ref(svGenotypeBool),
                    ref(vcfInfoMap)
                )
            );
        }
    }

    // Obtain the return value of function
    for (auto& futureOut : futureOutVec) {
        int result = move(futureOut.get());
    }

    // empty vector
    futureOutVec.clear();

    malloc_trim(0); // 0 is for heap memory

    // Close the thread pool
    pool.shutdown();

    cerr << endl;

    // save
    string outputFileName = sampleName + ".varigraph.vcf.gz";
    GENOTYPE::save(GraphMap, vcfHead, vcfInfoMap, sampleName, outputFileName, minSupportingGQ);

    return 0;
}

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
int GENOTYPE::for_bac_post_run( 
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
) {
    string chromosome = startNodeIter->first;  // chromosome

    // Check if the chromosome is present in the VCF file
    auto findIter1 = vcfInfoMap.find(chromosome);
    if (findIter1 == vcfInfoMap.end()) {  // if not
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << chromosome << "' does not exist in the VCF file." << endl;
        exit(1);
    }


    // Get the forward/backward iterator
    auto newStartNodeIterL = std::next(startNodeIter->second.begin(), threadStart);
    auto newStartNodeIterR = std::next(startNodeIter->second.begin(), threadEnd);
    
    // Check if the iterator is out of range
    if (newStartNodeIterL == startNodeIter->second.end()) {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Warning: the start iterator is pointing to the end -> " << chromosome << endl;
        return 0;
    }


    /* ************************************************** Haplotype selection ************************************************** */
    vector<uint16_t> topHapVec;  // Haplotypes screened by Dirichlet distribution
    unordered_map<uint16_t, double> hapIdxScoreMap;  // Likelihood of haplotype occurrence: map<hapIdx, possibility>
    haplotype_selection(
        chromosome, 
        newStartNodeIterL, 
        newStartNodeIterR, 
        hapMap, 
        haploidNum, 
        topHapVec, 
        hapIdxScoreMap
    );
    // Sort topHapVec in ascending order
    std::sort(topHapVec.begin(), topHapVec.end());


    /* ************************************************** 95% confidence interval for average coverage ************************************************** */
    /*
     To determine whether a k-mer is present in the Bloom filter but does not have a genotype of 0, we check if it falls within the 95% confidence interval of aveKmerCoverage. 
     If the k-mer does fall within this interval, it is considered to be present in the genotype 0 haplotype of the current node; otherwise, it is not present.
    */
    double lower = 256.0f;
    double upper = -0.1f;
    GENOTYPE::calculate_poisson_interval(hapKmerCoverage, lower, upper);


    /* ************************************************** forward ************************************************** */
    vector<HMMScore> preHMMScoreVec;  // Result of previous node forward/backward algorithm
    uint32_t preNodeStart = 0;  // The start position of the previous node
    uint32_t preNodeEnd = 0;  // The ending position of the previous node
    // forward
    for (map<uint32_t, nodeSrt>::iterator iter1 = newStartNodeIterL; iter1 != newStartNodeIterR; iter1++) {  // map<nodeStart, nodeSrt>
        if (iter1->second.hapGtVec.size() <= 1) continue;  // If small than 1, skip the node

        uint32_t nodeStart = iter1->first;  // the starting position of the node
        uint32_t nodeEnd = nodeStart + iter1->second.seqVec[0].size() - 1;  // end position of the node

        // Check if only structural variants are genotyped
        if (svGenotypeBool) {
            auto findIter2 = findIter1->second.find(nodeStart);
            // Check if the chromosome is present in the VCF file
            if (findIter2 == findIter1->second.end()) {  // if not
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << chromosome << ":" << nodeStart << "' does not exist in the VCF file." << endl;
                exit(1);
            }

            // Determine whether the length of REF or ALT is greater than 50
            if (findIter2->second[3].size() < 50 && findIter2->second[4].size() < 50) {
                continue;
            }
        }

        // hidden state
        vector<HHS> HHSStrVec = GENOTYPE::hidden_states(
            sampleType, 
            samplePloidy, 
            chromosome, 
            nodeStart, 
            startNodeIter, 
            iter1->second, 
            topHapVec, 
            hapIdxQRmap, 
            lower, 
            upper, 
            kmerLen, 
            true
        );

        // The anerave coverage and vector of node k-mers, 2023/12/04
        float NodeAveKmerCoverage = hapKmerCoverage;

        if (debugGenotype) {
            cerr << "start:" << nodeStart << endl;
            for (const auto& [hapVec, HiddenStatesVec, observableScore] : HHSStrVec) {
                if (hapVec.empty() || HiddenStatesVec.empty()) {
                    continue;
                }
                vector<uint16_t> hapVecTmp = hapVec;
                cerr << "hap:" << join(hapVecTmp, "/") << " h/c/f:\n";
                // h/c/f
                for (const auto& HiddenStates : HiddenStatesVec) {
                    cerr << setw(3) << +HiddenStates.h << "/" << +HiddenStates.c << "/" << +HiddenStates.f;
                }
                cerr << endl;
            }
        }

        // transition probability
        long double recombProb = 0.0L;  // Probability of recombination 
        long double noRecombProb = 0.0L;  // Probability of non-recombination
        // Determine transition probability type
        if (transitionProType == "rec") {
            uint32_t nodeDistence = nodeStart - preNodeEnd;  // The distance from the previous node
            tie(recombProb, noRecombProb) = transition_probabilities(
                nodeDistence, 
                hapMap
            );
            if (polyMixGenotype) {  // per-copy switch: stay exp(-rate*d), switch to each candidate haplotype (1-stay)/H
                noRecombProb = std::exp(-(long double)polyMixSwitchRate * (long double)nodeDistence);
                recombProb = (1.0L - noRecombProb) / (long double)std::max<size_t>(1, topHapVec.size());
            }
        }

        // observation matrix
        GENOTYPE::observable_states(
            NodeAveKmerCoverage, 
            HHSStrVec, 
            iter1->second
        );

        if (debugGenotype) {
            cerr << "start:" << nodeStart << endl;
            for (const auto& [hapVec, HiddenStatesVec, observableScore] : HHSStrVec) {
                if (hapVec.empty() || HiddenStatesVec.empty()) {
                    continue;
                }
                vector<uint16_t> hapVecTmp = hapVec;
                cerr << "hap:" << join(hapVecTmp, "/") << " observableStates:" << observableScore << endl;
            }
        }


        // Variable binding, used to modify the score of the forward algorithm in the node corresponding to nodeStart
        auto& HMMScoreVecTmp = startNodeIter->second.at(nodeStart).HMMScoreVec;
        forward(
            preHMMScoreVec, 
            hapIdxScoreMap, 
            recombProb, 
            noRecombProb, 
            HHSStrVec, 
            HMMScoreVecTmp
        );

        if (debugGenotype) {
            cerr << "start:" << nodeStart << endl;
            for (const auto& HMMScoreTmp : HMMScoreVecTmp) {
                cerr << "hap:";
                for (size_t i = 0; i < HMMScoreTmp.hapVec.size(); i++) {
                    if (i == 0) {
                        cerr << +HMMScoreTmp.hapVec[i];
                    } else {
                        cerr << "/" << +HMMScoreTmp.hapVec[i];
                    }
                }
                cerr << " Alpha:" << HMMScoreTmp.a << endl;
            }
        }

        // Reset the information of the previous node
        preNodeStart = nodeStart;
        preNodeEnd = nodeEnd;
        preHMMScoreVec = HMMScoreVecTmp;
    }


    /* ************************************************** backward ************************************************** */
    // variable zeroing
    vector<HMMScore>().swap(preHMMScoreVec);  // Result of previous node forward/backward algorithm
    preNodeStart = 0;  // The start position of the previous node
    preNodeEnd = 0;  // The ending position of the previous node
    for (std::map<uint32_t, nodeSrt>::reverse_iterator iter1 = std::map<uint32_t, nodeSrt>::reverse_iterator(newStartNodeIterR); iter1 != std::map<uint32_t, nodeSrt>::reverse_iterator(newStartNodeIterL); ++iter1) {
        if (iter1->second.hapGtVec.size() <= 1) continue;  // If small than 1, skip the node

        uint32_t nodeStart = iter1->first;  // the starting position of the node
        uint32_t nodeEnd = nodeStart + iter1->second.seqVec[0].size() - 1;  // end position of the node

        // Check if only structural variants are genotyped
        if (svGenotypeBool) {
            auto findIter2 = findIter1->second.find(nodeStart);
            // Check if the chromosome is present in the VCF file
            if (findIter2 == findIter1->second.end()) {  // if not
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << chromosome << ":" << nodeStart << "' does not exist in the VCF file." << endl;
                exit(1);
            }

            // Determine whether the length of REF or ALT is greater than 50
            if (findIter2->second[3].size() < 50 && findIter2->second[4].size() < 50) {
                continue;
            }
        }

        // hidden state
        vector<HHS> HHSStrVec = GENOTYPE::hidden_states(
            sampleType, 
            samplePloidy, 
            chromosome, 
            nodeStart, 
            startNodeIter, 
            iter1->second, 
            topHapVec, 
            hapIdxQRmap, 
            lower, 
            upper, 
            kmerLen, 
            false
        );

        // The anerave coverage and vector of node k-mers, 2023/12/04
        float NodeAveKmerCoverage = hapKmerCoverage;

        // transition probability
        long double recombProb = 0.0L;  // Probability of recombination 
        long double noRecombProb = 0.0L;  // Probability of non-recombination
        // Determine transition probability type
        if (transitionProType == "rec") {
            uint32_t nodeDistence = preNodeStart - nodeEnd;  // The distance from the previous node
            tie(recombProb, noRecombProb) = transition_probabilities(
                nodeDistence, 
                hapMap
            );
            if (polyMixGenotype) {  // per-copy switch: stay exp(-rate*d), switch to each candidate haplotype (1-stay)/H
                noRecombProb = std::exp(-(long double)polyMixSwitchRate * (long double)nodeDistence);
                recombProb = (1.0L - noRecombProb) / (long double)std::max<size_t>(1, topHapVec.size());
            }
        }
        
        // observation matrix
        GENOTYPE::observable_states(
            NodeAveKmerCoverage, 
            HHSStrVec, 
            iter1->second
        );


        // Variable binding, used to modify the score of the backward algorithm in the node corresponding to nodeStart
        auto& HMMScoreVecTmp = startNodeIter->second.at(nodeStart).HMMScoreVec;
        backward(
            preHMMScoreVec, 
            hapIdxScoreMap, 
            recombProb, 
            noRecombProb, 
            HHSStrVec, 
            HMMScoreVecTmp
        );

        if (debugGenotype) {
            for (const auto& HMMScoreTmp : HMMScoreVecTmp) {
                cerr << "start:" << nodeStart << " genotype:";
                for (size_t i = 0; i < HMMScoreTmp.hapVec.size(); i++) {
                    if (i == 0) {
                        cerr << +HMMScoreTmp.hapVec[i];
                    } else {
                        cerr << "/" << +HMMScoreTmp.hapVec[i];
                    }
                }
                cerr << " Beta:" << HMMScoreTmp.b << endl;
            }
        }

        // Reset the information of the previous node
        preNodeStart = nodeStart;
        preNodeEnd = nodeEnd;
        preHMMScoreVec = HMMScoreVecTmp;
    }

    /* ************************************************** posterior ************************************************** */
    for (map<uint32_t, nodeSrt>::iterator iter1 = newStartNodeIterL; iter1 != newStartNodeIterR; iter1++) {  // map<nodeStart, nodeSrt>
        if (iter1->second.hapGtVec.size() <= 1) continue;  // If small than 1, skip the node

        // Check if only structural variants are genotyped
        if (svGenotypeBool) {
            auto findIter2 = findIter1->second.find(iter1->first);
            // Check if the chromosome is present in the VCF file
            if (findIter2 == findIter1->second.end()) {  // if not
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << chromosome << ":" << iter1->first << "' does not exist in the VCF file." << endl;
                exit(1);
            }

            // Determine whether the length of REF or ALT is greater than 50
            if (findIter2->second[3].size() < 50 && findIter2->second[4].size() < 50) {
                continue;
            }
        }

        // Calculate and store the result in a local variable
        posterior(iter1->first, &(iter1->second), topHapVec);
    }

    return 0;
}


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
void GENOTYPE::haplotype_selection(
    const string& chromosome, 
    const map<uint32_t, nodeSrt>::iterator& newStartNodeIterL, 
    const map<uint32_t, nodeSrt>::iterator& newStartNodeIterR, 
    const map<uint16_t, string>& hapMap, 
    const uint32_t haploidNum, 
    vector<uint16_t>& topHapVec, 
    unordered_map<uint16_t, double>& hapIdxScoreMap
) {
    // If the number of haplotypes is less than haploidNum, skip screening
    if (hapMap.size() <= haploidNum) {
        for (const auto& pair : hapMap) {
            topHapVec.push_back(pair.first);
        }
    }
    
    // k-mer count of all haplotypes
    vector<uint32_t> hapKmerCountVec(hapMap.size(), 0);

    for (map<uint32_t, nodeSrt>::iterator nodeIter = newStartNodeIterL; nodeIter != newStartNodeIterR; nodeIter++) {  // map<nodeStart, nodeSrt>

        if (nodeIter->second.hapGtVec.size() == 1) {continue;}  // If only 0, skip the node

        // Traverse all nodes kmer
        for (const auto& GraphKmerHashHapStrMapIter : nodeIter->second.GraphKmerHashHapStrMapIterVec) {

            // k-mer coverage and frequence
            const uint8_t& kmerCov = GraphKmerHashHapStrMapIter->second.c;
            const uint8_t& kmerFre = GraphKmerHashHapStrMapIter->second.f;

            // if it is 0, skip the k-mer
            // if (kmerCov == 0 || kmerFre > 1) {continue;}
            if (kmerCov <= 1 || kmerFre > 1) {continue;}  // 2023/09/27
            
            //  record the k-mer frequency of each haplotype
            // haplotype information
            vector<bitset<8> > bitsVector;
            for (const auto& hapBit : GraphKmerHashHapStrMapIter->second.BitVec) {
                bitsVector.push_back(hapBit);
            }
            
            // haplotype
            for (const auto& [hapIdx, _] : hapMap) {  // map<hapIdx, hapName>

                uint16_t quotient1 = DIVIDE_BY_8(hapIdx);  // Quotient
                uint16_t remainder1 = GET_LOW_3_BITS(hapIdx);  // get remainder

                // If the haplotype contains the k-mer then add the coverage of the k-mer
                if (bitsVector[quotient1][remainder1]) {
                    hapKmerCountVec[hapIdx] += kmerCov;
                }
            }
        }
    }

    // The Dirichlet distribution calculates the probability of occurrence of each haplotype
    HaplotypeSelect HaplotypeSelectClass(hapKmerCountVec);
    HaplotypeSelectClass.calculate_sparsity();
    HaplotypeSelectClass.generate_sparse_frequency_vector();
    hapIdxScoreMap = HaplotypeSelectClass.get_top_indices(haploidNum);  // map<hapIdx, possibility>
    if (topHapVec.empty()) {
        topHapVec = HaplotypeSelectClass.mTopHapVec;
    }

    std::lock_guard<std::mutex> mtx_locker(mtxG);
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Haplotype selection results for " << chromosome << "-" << newStartNodeIterL->first << ":";  // print log
    for (uint16_t i = 0; i < topHapVec.size(); i++) {
        uint16_t hapIdxTmp = topHapVec[i];
        if (i == 0) {
            cerr << " " << hapIdxTmp;
        } else {
            cerr << ", " << hapIdxTmp ;
        }
    }
    cerr << endl;
}


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
vector<HHS> GENOTYPE::hidden_states(
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
) {
    // all haplotype combinations
    vector<HHS> HHSStrVec;  // vector<HHS>

    // Traverse all node k-mer, construct hidden state
    vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator> GraphKmerHashHapStrMapIterVec;  // Iterator pointing to mGraphKmerHashHapStrMap, vector<iter>

    // Genotype information for each haplotype at this locus: vector<GT>
    const vector<uint16_t>& hapGtVec = node.hapGtVec;

    // Determine whether the haplotype needs to rebuild the k-mers index 
    map<uint16_t, uint32_t> hapNeedIdxNumMap;  // map<hapIdx, number>, The number of HS similar to 2/0/2


    /* ============================================ Hidden states ============================================ */

    // Build a hidden matrix
    uint16_t maxHapIdx = hapIdxQRmap.size() - 1;  // Maximum haplotype index
    vector<vector<uint16_t> > ComHapVec = increment_vector(topHapVec, sampleType, samplePloidy, maxHapIdx);
    for (const auto& hapVec : ComHapVec) {
        HHS HHSStr;
        // Record haplotypes combined into hapVec
        HHSStr.hapVec = hapVec;
        HHSStrVec.push_back(HHSStr);
    }

    // All k-mer information of this node
    for (const auto& GraphKmerHashHapStrMapIter : node.GraphKmerHashHapStrMapIterVec) {  // vector<iter>
        // k-mer information
        const uint8_t& c = GraphKmerHashHapStrMapIter->second.c;
        const uint8_t& f = GraphKmerHashHapStrMapIter->second.f;
        const vector<int8_t>& graphHapBit = GraphKmerHashHapStrMapIter->second.BitVec;

        int lastBit = construct_index::get_bit(graphHapBit.back(), 7);

        vector<bitset<8> > graphBitsVector;

        for (const auto& hapBit : graphHapBit) {
            graphBitsVector.push_back(hapBit);
        }

        // Traverse topHapVec to check whether the k-mer frequency is greater than 1
        if (filter) {
            uint64_t kmerFrq = 0;
            for (const auto& hapIdx : topHapVec) {  // vector<hapIdx>
                kmerFrq += graphBitsVector[get<0>(hapIdxQRmap.at(hapIdx))][get<1>(hapIdxQRmap.at(hapIdx))];
            }

            // If the frequency is 0, skip the k-mer
            if (kmerFrq == 0) {continue;}

            // Record the k-mer
            GraphKmerHashHapStrMapIterVec.push_back(GraphKmerHashHapStrMapIter);
        } else {
            // Record the k-mer
            GraphKmerHashHapStrMapIterVec.push_back(GraphKmerHashHapStrMapIter);
        }

        // Merge k-mer information into HHSStrVec
        for (HHS& HHSStr : HHSStrVec) {  // vector<HHS>

            // Temporary HS structure
            HS HSStr;
            uint8_t& h = HSStr.h;
            HSStr.c = c;

            for (const auto& hapIdx : HHSStr.hapVec) {  // vector<hapIdx>

                const auto& QR = hapIdxQRmap.at(hapIdx);

                // If the k-mer does fall within this interval, it is considered to be present in the genotype 0 haplotype of the current node; otherwise, it is not present.
                uint8_t hTmp = (lastBit == 1 && hapGtVec[hapIdx] == 0 && c >= lower && c <= upper) ? 1 : graphBitsVector[get<0>(QR)][get<1>(QR)];
                h += hTmp;

                // Check if the haplotype truly contains the k-mer at the current node, as it may be present in another node for the haplotype.
                if (hTmp > 0 && c < lower && f >= 2) {  // 1/0/2; 2/0/2
                    // Record the haplotype
                    auto emplacedValue = hapNeedIdxNumMap.emplace(hapIdx, 0).first;
                    emplacedValue->second++;
                }
            }

            // frequency
            uint8_t frequency = f;
            if (lastBit == 1 && frequency == 1) {
                ++frequency;
            }
            HSStr.f = frequency;

            // Record in the HiddenStatesVec of the haplotype combination
            HHSStr.HiddenStatesVec.push_back(move(HSStr));
        }
    }

    /* ============================================ Determine whether the haplotype needs to rebuild the k-mers index ============================================ */
    if (!hapNeedIdxNumMap.empty()) {

        // k-mers index
        unordered_map<uint16_t, unordered_set<uint64_t> > hapKmerHashSetMap;  // map<hapIdx, set<kmerHash> >

        const auto& seqVec = node.seqVec;  // Store sequence information of query: vector<seq>

        for (const auto& [hapIdx, _] : hapNeedIdxNumMap) {  // map<hapIdx, number>
            uint16_t gt = hapGtVec[hapIdx];

            // Check if the array is out of bounds
            if (gt >= seqVec.size()) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: Node '"<< chromosome << "-" << nodeStart << "' does not contain sequence information for haplotype " << gt << "." << endl;
                exit(1);
            }
            // Get ALT sequence
            string seqTmp = seqVec[gt];

            // Sequence with upstream and downstream 1 k-mer sequence
            pair<string, string> upDownSeq = construct_index::find_node_up_down_seq(
                hapIdx, 
                gt, 
                seqTmp, 
                kmerLen - 1, 
                startNodeIter->second.find(nodeStart), 
                startNodeIter->second
            );
            
            seqTmp = upDownSeq.first + seqTmp + upDownSeq.second;  // Add upstream and downstream sequences
            
            // k-mer indexing
            unordered_set<uint64_t> hashSet = kmerBit::kmer_sketch_genotype(seqTmp, kmerLen);
            hapKmerHashSetMap[hapIdx] = hashSet;
        }
        
        // Update the bitmap of the node
        uint32_t HiddenStatesStrVecIdx = 0;

        for (auto& GraphKmerHashHapStrMapIter : GraphKmerHashHapStrMapIterVec) {  // All k-mer information of this node
            const uint64_t& kmerHash = GraphKmerHashHapStrMapIter->first;
            const uint8_t& c = GraphKmerHashHapStrMapIter->second.c;
            const uint8_t& f = GraphKmerHashHapStrMapIter->second.f;
            const vector<int8_t>& graphHapBit = const_cast<vector<int8_t>&>(GraphKmerHashHapStrMapIter->second.BitVec);

            // When c is greater than 0 or f is equal to 1 , skip
            if (c > lower || f <= 1) {HiddenStatesStrVecIdx++; continue;}

            int lastBit = construct_index::get_bit(graphHapBit.back(), 7);

            vector<bitset<8> > graphBitsVector;

            for (const auto& hapBit : graphHapBit) {
                graphBitsVector.push_back(hapBit);
            }

            // Merge k-mer information into HHSStrVec
            for (HHS& HHSStr : HHSStrVec) {  // vector<HHS>

                // Temporary HS structure
                HS& HSStr = HHSStr.HiddenStatesVec[HiddenStatesStrVecIdx];
                uint8_t& h = HSStr.h;

                for (const auto& hapIdx : HHSStr.hapVec) {  // vector<hapIdx>

                    const auto& QR = hapIdxQRmap.at(hapIdx);

                    // If the k-mer does fall within this interval, it is considered to be present in the genotype 0 haplotype of the current node; otherwise, it is not present.
                    uint8_t hTmp = (lastBit == 1 && hapGtVec[hapIdx] == 0 && c >= lower && c <= upper) ? 1 : graphBitsVector[get<0>(QR)][get<1>(QR)];

                    // Iterator of haplotype
                    auto hapIdxFindIter = hapKmerHashSetMap.find(hapIdx);

                    // If neither hap1is in the hapKmerHashSetMap, skip the loop
                    if (hapIdxFindIter == hapKmerHashSetMap.end() || hTmp == 0) {continue;}

                    // If haplotype does not contain the k-mer, subtract 1 from h
                    if (hTmp == 1 && hapIdxFindIter !=  hapKmerHashSetMap.end() && hapIdxFindIter->second.find(kmerHash) == hapIdxFindIter->second.end()) {
                        if (h >= 1) {h--;}
                    }
                }
            }

            // Index increment
            HiddenStatesStrVecIdx++;
        }
    }

    // Node's k-mer reassignment
    if (filter) {
        std::lock_guard<std::mutex> mtx_locker(mtxG);
        node.GraphKmerHashHapStrMapIterVec = GraphKmerHashHapStrMapIterVec;
    }

    return HHSStrVec;
}


/**
 * @author zezhen du
 * @brief Gets the combination of all haplotypes
 * 
 * @param hapVec        Vector of haplotypes     
 * @param sampleType    specify the genotype of the sample genome (hom/het)
 * @param samplePloidy  sample ploidy (2-8) [2]
 * @param maxHapIdx     Maximum haplotype index
 * 
 * @return ComHapVec
**/
vector<vector<uint16_t> > GENOTYPE::increment_vector(
    const vector<uint16_t>& hapVec, 
    const string& sampleType, 
    const uint32_t& samplePloidy, 
    const uint16_t& maxHapIdx
) {
    // combination of haplotype name
    vector<vector<uint16_t> > ComHapVec;

    /* ----------------------------------------------------------- Polyploidy ----------------------------------------------------------- */
    // --poly-mix: any multiset of the selected haplotypes (the original branch below only allows the
    // ploidy haplotypes of one VCF sample together, so a new sample can only copy a panel individual).
    if (samplePloidy > 2 && polyMixGenotype) {
        return multiset_states(hapVec, samplePloidy, sampleType == "hom");
    }
    // Polyploidy, sample haplotype combination
    if (samplePloidy > 2) {
        // sample haplotype
        for (const auto& hap : hapVec) {
            vector<uint16_t> hapVecTmp(samplePloidy);

            // If it is the reference genome, then all are 0
            if (hap == 0) {
                hapVecTmp.assign(samplePloidy, 0);
            } else {
                int32_t quotient = ceil(hap / (float) samplePloidy);
                uint16_t firstHap = (quotient - 1) * samplePloidy + 1;
            
                iota(hapVecTmp.begin(), hapVecTmp.end(), firstHap);

                // Ensure haplotype index greater than maxHapIdx are set to 0
                for (auto& val : hapVecTmp) {
                    if (val > maxHapIdx) val = 0;
                }
            }
            ComHapVec.push_back(std::move(hapVecTmp));
        }

        // sort and deduplication
        std::set<std::vector<uint16_t> > ComHapVecSet(ComHapVec.begin(), ComHapVec.end());
        ComHapVec.assign(ComHapVecSet.begin(), ComHapVecSet.end());

        return ComHapVec;
    }


    /* ----------------------------------------------------------- diploid ----------------------------------------------------------- */
    // number of haplotype
	uint32_t hapNum = hapVec.size() - 1;

    // combination of haplotype index (selected vector)
	vector<vector<uint32_t> > hapIdxVecIdxVec;

	for (uint32_t hapIdx = 0; hapIdx < hapVec.size(); hapIdx++) {
        // homozygous
		std::vector<uint32_t> vec(samplePloidy, hapIdx);
		hapIdxVecIdxVec.push_back(vec);

        // if the sample is homozygous, the next loop proceeds
		if (sampleType == "hom") continue;

		auto minElement = std::min_element(vec.begin() + 1, vec.end());
		while (*minElement < hapNum) {
			uint32_t index = vec.size() - 1;
			while (vec[index] == hapNum) {
				vec[index] = *minElement + 1;
				index--;
			}
			vec[index]++;

			hapIdxVecIdxVec.push_back(vec);

			minElement = std::min_element(vec.begin() + 1, vec.end());
		}
	}

    // combination of haplotype name
    ComHapVec.reserve(hapIdxVecIdxVec.size());

	for (const auto& hapIdxVec : hapIdxVecIdxVec) {
		vector<uint16_t> hapVecTmp;
        hapVecTmp.reserve(hapIdxVec.size());
		for (const auto& hapIdx : hapIdxVec) {
			hapVecTmp.push_back(hapVec[hapIdx]);
		}
		ComHapVec.push_back(move(hapVecTmp));
	}

	return ComHapVec;
}


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
void GENOTYPE::calculate_poisson_interval(const double& lambda, double& lower, double& upper) {
    // Calculate the standard deviation
    double stdDev = std::sqrt(lambda);

    // Calculate the upper limit
    upper = lambda + 1.96 * stdDev;

    // Calculate the lower limit
    lower = lambda - 1.96 * stdDev;
}

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
pair<long double, long double> GENOTYPE::transition_probabilities(const uint32_t & nodeDistence, const map<uint16_t, string> & hapMap)
{
    // transition probability
    double effectivePopulationSize = 1e-05;  // effective population size
    uint16_t populationSize = hapMap.size();  // population size
    double recombRate = 1.26;  // recombination rate
    long double distance = nodeDistence * 0.000004L * ((long double) recombRate) * effectivePopulationSize;  // d
    long double recombProb = (1.0L - exp(-distance / (long double) populationSize) )* (1.0L / (long double) populationSize);
    long double noRecombProb = exp(-distance / (long double) populationSize) + recombProb;
    return make_pair(recombProb, noRecombProb);  // Returns the probabilities of recombination and nonrecombination
}


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
void GENOTYPE::observable_states(
    const float & aveKmerCoverage, 
    vector<HHS>& HHSStrVec, 
    nodeSrt & node
) {
    // 95% confidence interval for aveKmerCoverage
    double lower = 256.0f;
    double upper = -0.1f;
    GENOTYPE::calculate_poisson_interval(aveKmerCoverage, lower, upper);

    /* ============================================ Calculate the observation matrix ============================================ */
    // --poly-mix, biallelic node: emission marginalises the true dosage d given the state's alt-copy count a,
    // allowing each copy to differ from the copied haplotype with probability polyMixMismatch.
    if (polyMixGenotype && !debugGenotype) {
        const auto& gtv = node.hapGtVec;
        bool biallelic = true;
        for (auto g : gtv) if (g > 1) { biallelic = false; break; }
        const HHS* any = nullptr;
        for (const auto& s : HHSStrVec) if (!s.hapVec.empty() && !s.HiddenStatesVec.empty()) { any = &s; break; }
        if (biallelic && any != nullptr && any->hapVec.size() > 2) {
            const int P = (int)any->hapVec.size();
            const size_t nK = any->HiddenStatesVec.size();
            vector<double> mAlt(nK, 0.0), mRef(nK, 0.0); int nAlt = 0, nRef = 0;
            for (const auto& s : HHSStrVec) {
                if (s.hapVec.empty() || s.HiddenStatesVec.size() != nK) continue;
                bool hom = std::all_of(s.hapVec.begin(), s.hapVec.end(), [&](uint16_t x) { return x == s.hapVec[0]; });
                if (!hom || s.hapVec[0] >= gtv.size()) continue;
                bool alt = gtv[s.hapVec[0]] == 1;
                auto& m = alt ? mAlt : mRef; (alt ? nAlt : nRef)++;
                for (size_t k = 0; k < nK; ++k) m[k] += s.HiddenStatesVec[k].h / (double)P;
            }
            if (nAlt > 0 && nRef > 0) {
                for (size_t k = 0; k < nK; ++k) { mAlt[k] /= nAlt; mRef[k] /= nRef; }
                MixNode mn; mn.P = P;
                // Allele-specific k-mers (present in all alt haplotypes and no ref haplotype, or vice versa, and
                // unique in the graph): their mean coverage acts as allele read depth. Binomial in the alt fraction
                // d/P is scale-free, so it does not depend on the global k-mer coverage estimate, and zero alt
                // k-mers is evidence for low dosage (the geometric h=0 term barely differs from h=1 at low depth).
                double altSum = 0, refSum = 0; int altN = 0, refN = 0;
                for (size_t k = 0; k < nK; ++k) {
                    const auto& hs = any->HiddenStatesVec[k];
                    if (hs.f > 1) continue;
                    if (mAlt[k] > 0.999 && mRef[k] < 0.001) { altSum += hs.c; ++altN; }
                    else if (mRef[k] > 0.999 && mAlt[k] < 0.001) { refSum += hs.c; ++refN; }
                }
                vector<long double> logL(P + 1, 0.0L);
                if (altN > 0 && refN > 0) {
                    const long double ea = altSum / altN, er = refSum / refN, err = 0.01L;
                    for (int d = 0; d <= P; ++d) {
                        long double fr = (long double)d / P; fr = fr * (1.0L - err) + (1.0L - fr) * err;
                        logL[d] = ea * std::log(fr) + er * std::log(1.0L - fr);
                    }
                }
                long double mx = *std::max_element(logL.begin(), logL.end());
                mn.L.resize(P + 1);
                for (int d = 0; d <= P; ++d) mn.L[d] = std::exp(logL[d] - mx);
                mn.PdA = dosage_given_alt_copies(P, polyMixMismatch);
                mn.M.assign(P + 1, 0.0L);
                for (int a = 0; a <= P; ++a) for (int d = 0; d <= P; ++d) mn.M[a] += mn.PdA[a][d] * mn.L[d];
                for (auto& s : HHSStrVec) {
                    int a = 0;
                    for (auto x : s.hapVec) a += (x < gtv.size() && gtv[x] == 1);
                    s.observableScore = mn.M[a];
                }
                mixNodeCache[&node] = std::move(mn);
                return;
            }
        }
    }
    // --poly-mix, other nodes: every state shares the same k-mer list (same order, same c/f); only h differs.
    // Precompute the per-k-mer score for h = 0..maxH once, then each state is a product of lookups.
    if (polyMixGenotype && !debugGenotype) {
        size_t nK = 0; uint8_t maxH = 0;
        for (const auto& HHSStr : HHSStrVec) {
            if (HHSStr.HiddenStatesVec.empty()) continue;
            nK = HHSStr.HiddenStatesVec.size();
            for (const auto& hs : HHSStr.HiddenStatesVec) maxH = std::max(maxH, hs.h);
        }
        const HHS* ref = nullptr;
        for (const auto& HHSStr : HHSStrVec) if (HHSStr.HiddenStatesVec.size() == nK && nK > 0) { ref = &HHSStr; break; }
        if (ref != nullptr) {
            const size_t W = (size_t)maxH + 1;
            vector<long double> tab(nK * W);
            for (size_t k = 0; k < nK; ++k) {
                for (uint8_t h = 0; h <= maxH; ++h) {
                    uint8_t c = ref->HiddenStatesVec[k].c, f = ref->HiddenStatesVec[k].f;
                    find_most_likely_depth(h, c, f, aveKmerCoverage, lower, upper);
                    tab[k * W + h] = (h == 0) ? geometric(get_error_param(aveKmerCoverage), c) : poisson(aveKmerCoverage * h, c);
                }
            }
            for (auto& HHSStr : HHSStrVec) {
                auto& result = HHSStr.observableScore;
                result = 1.0L;
                const auto& hsv = HHSStr.HiddenStatesVec;
                for (size_t k = 0; k < hsv.size() && k < nK; ++k) result *= tab[k * W + hsv[k].h];
            }
            return;
        }
    }
    for (auto& HHSStr : HHSStrVec) {
        auto& result = HHSStr.observableScore;  // save result
        result = 1.0L;  // The initial score is 1

        for (const auto& HiddenStatesStr : HHSStr.HiddenStatesVec) {  // vector<HS>
            uint8_t h = HiddenStatesStr.h;  // k-mer hidden state
            uint8_t c = HiddenStatesStr.c;  // k-mer coverage
            uint8_t f = HiddenStatesStr.f;  // k-mer frequency

            // Most likely k-mer frequency
            find_most_likely_depth(h, c, f, aveKmerCoverage, lower, upper);

            if (h == 0) {
                if (debugGenotype) {
                    cerr << join(HHSStr.hapVec, "/") << ": " << +h << " " << aveKmerCoverage << " " << +c << " " << geometric(get_error_param(aveKmerCoverage), c) << " " << result << endl;
                }

                result *= geometric(get_error_param(aveKmerCoverage), c);
            } else {
                if (debugGenotype) {
                    cerr << join(HHSStr.hapVec, "/") << ": " << +h << " " << aveKmerCoverage << " " << +c << " " << poisson(aveKmerCoverage * h, c) << " " << result << endl;
                }
                
                result *= poisson(aveKmerCoverage * h, c);
            }
        }
    }
}

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
long double GENOTYPE::poisson(
    long double mean, 
    uint8_t value
) {
    long double sum = 0.0L;
    int v = (int) value;
    for (size_t i = 1; i <= value; ++i) sum += log(i);
    long double log_val = -mean + v * log(mean) - sum;
    return exp(log_val);
}

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
double GENOTYPE::get_error_param(double aveKmerCoverage)
{
    double cn0;
    if (aveKmerCoverage < 10.0) {
        cn0 = 0.99;
    } else if (aveKmerCoverage < 20) {
        cn0 = 0.95;
    } else if (aveKmerCoverage < 40) {
        cn0 = 0.9;
    } else {
        cn0 = 0.8;
    }
    return cn0;
}

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
// long double GENOTYPE::geometric(long double p, uint8_t value) {
//     // return pow(1.0L - p, value)*p;
//     return pow(1.0L - p, value/2)*p;
// }


/**
 * The Bayes formula is as follows:
 * P(p|data) = P(data|p) * P(p) / P(data)

 * Among them:
 * P(p|data) Is a posterior probability, representing the probability distribution given the data;
 * P(data|p) Is the likelihood function, which represents the probability of observing the data given the accuracy rate p;
 * P(p) Is the prior probability, which represents the prior belief of the correct rate p;
 * P(data) Is a standardized constant used to ensure that the sum of the posterior probabilities is 1.
**/

// Posterior probability calculation
long double GENOTYPE::geometric(long double p, uint8_t value) {
    long double prior_prob = prior(p);
    long double likelihood_prob = likelihood(p, value);
    
    // Calculate the posterior probability
    long double posterior_prob = (likelihood_prob * prior_prob);
    
    return posterior_prob;
}

// Prior probability distribution function (normal distribution)
long double GENOTYPE::prior(long double p) {
    long double mean = 0.5; // The mean of the normal distribution
    long double variance = 0.05; // Variance of the normal distribution
    long double prior_prob = (1 / (std::sqrt(2 * M_PI * variance))) * std::exp(-std::pow(p - mean, 2) / (2 * variance));
    return prior_prob;
}

// Likelihood function
long double GENOTYPE::likelihood(long double p, uint8_t value) {
    // Use the binomial distribution as the likelihood function
    long double q = 1.0 - p; // The probability of reading as sequencing error
    // long double likelihood = (std::pow(q, value/4)) * (std::pow(p, (1 - value)));  // 2023/09/14 -> better for plant's genomes
    long double likelihood = (std::pow(q, value)) * (std::pow(p, (1 - value)));  // 2023/09/14 -> better for plant's and animal's genomes
    return likelihood;
}

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
void GENOTYPE::find_most_likely_depth(
    const uint8_t& h, 
    uint8_t& c, 
    const uint8_t& f, 
    const float& aveKmerCoverage, 
    const double& lower, 
    const double& upper
) {
    // 2023/12/04 - ho
    if (f == 1) {  // h/c/f 2/21/1 21
        return;
    } else {
        if (h > 0 && c > (aveKmerCoverage * h)) {  // h/c/f 2/46/2 23 -> 2/23/2.   h/c/f 1/46/2 23 -> 2/11.5/2.
            c = aveKmerCoverage * h;
        } else if (h == 0 && c > aveKmerCoverage) {  // h/c/f 0/60/2 23 -> 0/0/2  2023/09/14 -> better for plant genomes
            c = (f > ((float)c / upper)) ? 0 : c / (float)f;
        } else if (h == 0 && c <= aveKmerCoverage) {  // h/c/f 0/11/2 21 -> 0/5.5/2  2023/09/14 -> better for plant genomes
            c /= (float)f;
        } else {
            return;
        }
    }
}

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
void GENOTYPE::forward(
    const vector<HMMScore>& preHMMScoreVec, 
    const unordered_map<uint16_t, double>& hapIdxScoreMap, 
    const long double& recombProb, 
    const long double& noRecombProb, 
    const vector<HHS>& HHSStrVec, 
    vector<HMMScore>& HMMScoreVec
)
{
    vector<long double> alphaVecTmp;  // Temporary storage
    long double alphaSumTmp = 0.0L;  // Temporary storage

    if (polyMixGenotype && !preHMMScoreVec.empty() && !(recombProb == 0.0L && noRecombProb == 0.0L) &&
        mix_step(preHMMScoreVec, true, recombProb, noRecombProb, HHSStrVec, alphaVecTmp)) {
        HMMScoreVec.resize(alphaVecTmp.size());
        uint32_t i = 0;
        std::lock_guard<std::mutex> mtx_locker(mtxG);
        for (const auto& h : HHSStrVec) {
            if (h.hapVec.empty() || h.HiddenStatesVec.empty()) continue;
            HMMScoreVec[i].a = alphaVecTmp[i]; HMMScoreVec[i].hapVec = h.hapVec; ++i;
        }
        return;
    }

    for (const auto& [hapVec, HiddenStatesVec, observableScore] : HHSStrVec) {
        if (hapVec.empty() || HiddenStatesVec.empty()) {
            continue;
        }

        int32_t hapNum = hapVec.size();  // Number of haplotypes

        long double result = 0.0L;  // save result

        if (preHMMScoreVec.empty()) {  // If the first node
            result += observableScore;
        }
        else {  // not the first node
            for (const auto& HMMScoreTmp : preHMMScoreVec) {  // vector<HMMScore>
                // haplotype frequency
                if (recombProb == 0.0L && noRecombProb == 0.0L) {
                    long double resultTmp = HMMScoreTmp.a * observableScore;  // temporary score
                    for (const auto hapIdx : hapVec) {
                        auto findIter = hapIdxScoreMap.find(hapIdx);
                        if (findIter == hapIdxScoreMap.end()) {
                            cerr << "[" << __func__ << "::" << getTime() << "] " 
                                << "'" << hapIdx << "' does not exist in 'hapIdxScoreMap'." << endl;
                            exit(1);
                        }
                        resultTmp *= findIter->second;
                    }

                    // score
                    result += resultTmp;
                } else {  // recombination rate
                    // Calculate the number of intersects
                    vector<uint16_t> intersection;
                    std::set_intersection(hapVec.begin(), hapVec.end(),
                                        HMMScoreTmp.hapVec.begin(), HMMScoreTmp.hapVec.end(),
                                        std::back_inserter(intersection));
                    
                    int32_t noRecombNum = intersection.size();  // Number of noRecomb
                    int32_t recombNum = hapNum - noRecombNum;  // Number of recomb

                    // score
                    result += HMMScoreTmp.a * pow(noRecombProb, noRecombNum)  * pow(recombProb, recombNum) * observableScore;
                }
            }
        }
        
        alphaVecTmp.push_back(result);
        alphaSumTmp += result;
    }
    
    // Standardize to prevent the value from being too small and return to zero
    if (alphaSumTmp > 0.0L) {
        transform(alphaVecTmp.begin(), alphaVecTmp.end(), alphaVecTmp.begin(), bind(divides<long double>(), placeholders::_1, alphaSumTmp));
    } else {
        long double uniform = 1.0L / (long double) alphaVecTmp.size();
        transform(alphaVecTmp.begin(), alphaVecTmp.end(), alphaVecTmp.begin(), [uniform](long double c) -> long double {return uniform;});
    }

    // index
    HMMScoreVec.resize(alphaVecTmp.size());
    uint32_t indexTmp = 0;
    for (const auto& [hapVec, HiddenStatesVec, observableScore] : HHSStrVec) {
        if (hapVec.empty() || HiddenStatesVec.empty()) {
            continue;
        }

        // Multithreaded Data Lock
        std::lock_guard<std::mutex> mtx_locker(mtxG);
        HMMScoreVec[indexTmp].a = alphaVecTmp[indexTmp];
        HMMScoreVec[indexTmp].hapVec = hapVec;
        indexTmp++;
    }
}


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
void GENOTYPE::backward(
    const vector<HMMScore>& preHMMScoreVec, 
    const unordered_map<uint16_t, double>& hapIdxScoreMap, 
    const long double& recombProb, 
    const long double& noRecombProb, 
    const vector<HHS>& HHSStrVec, 
    vector<HMMScore>& HMMScoreVec
)
{
    vector<long double> betaVecTmp;  // Temporary storage
    long double betaSumTmp = 0.0L;  // Temporary storage

    if (polyMixGenotype && !preHMMScoreVec.empty() && !(recombProb == 0.0L && noRecombProb == 0.0L) &&
        mix_step(preHMMScoreVec, false, recombProb, noRecombProb, HHSStrVec, betaVecTmp)) {
        std::lock_guard<std::mutex> mtx_locker(mtxG);
        for (size_t i = 0; i < betaVecTmp.size() && i < HMMScoreVec.size(); ++i) HMMScoreVec[i].b = betaVecTmp[i];
        return;
    }

    for (const auto& [hapVec, HiddenStatesVec, observableScore] : HHSStrVec) {
        if (hapVec.empty() || HiddenStatesVec.empty()) {
            continue;
        }

        int32_t hapNum = hapVec.size();  // Number of haplotypes

        long double result = 0.0L;  // save result

        if (preHMMScoreVec.empty()) {  // If the first node
            result += observableScore;
        }
        else {  // not the first node
            for (const auto& HMMScoreTmp : preHMMScoreVec) {  // vector<HMMScore>
                // haplotype frequency
                if (recombProb == 0.0L && noRecombProb == 0.0L) {
                    long double resultTmp = HMMScoreTmp.b * observableScore;  // temporary score
                    for (const auto hapIdx : hapVec) {
                        auto findIter = hapIdxScoreMap.find(hapIdx);
                        if (findIter == hapIdxScoreMap.end()) {
                            cerr << "[" << __func__ << "::" << getTime() << "] " 
                                << "'" << hapIdx << "' does not exist in 'hapIdxScoreMap'." << endl;
                            exit(1);
                        }
                        resultTmp *= findIter->second;
                    }
                    
                    // score
                    result += resultTmp;
                } else {  // recombination rate
                    // Calculate the number of intersects
                    vector<uint16_t> intersection;
                    std::set_intersection(hapVec.begin(), hapVec.end(),
                                        HMMScoreTmp.hapVec.begin(), HMMScoreTmp.hapVec.end(),
                                        std::back_inserter(intersection));
                    
                    int32_t noRecombNum = intersection.size();  // Number of noRecomb
                    int32_t recombNum = hapNum - noRecombNum;  // Number of recomb

                    // score
                    result += HMMScoreTmp.b * pow(noRecombProb, noRecombNum) * pow(recombProb, recombNum) * observableScore;
                }
            }
        }

        betaVecTmp.push_back(result);
        betaSumTmp += result;
     }

    // Standardize to prevent the value from being too small and return to zero
    if (betaSumTmp > 0.0L) {
        transform(betaVecTmp.begin(), betaVecTmp.end(), betaVecTmp.begin(), bind(divides<long double>(), placeholders::_1, betaSumTmp));
    } else {
        long double uniform = 1.0L / (long double) betaVecTmp.size();
        transform(betaVecTmp.begin(), betaVecTmp.end(), betaVecTmp.begin(), [uniform](long double c) -> long double {return uniform;});
    }

    // index
    uint32_t indexTmp = 0;
    for (const auto& [hapVec, HiddenStatesVec, observableScore] : HHSStrVec) {
        if (hapVec.empty() || HiddenStatesVec.empty()) {
            continue;
        }

        // Multithreaded Data Lock
        std::lock_guard<std::mutex> mtx_locker(mtxG);
        HMMScoreVec[indexTmp].b = betaVecTmp[indexTmp];
        indexTmp++;
    }
}

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
int GENOTYPE::posterior(
    uint32_t nodeStart, 
    nodeSrt* nodePtr, 
    vector<uint16_t>& topHapVec
) {
    if (nodePtr == nullptr) {return 0;}

    // Number of node-unique k-mers, UK
    uint8_t uniqueKmerNum = get_UK((*nodePtr).GraphKmerHashHapStrMapIterVec);

    // haplotype genotype information
    const auto& hapGtVec = (*nodePtr).hapGtVec;

    // map<hapIdx, tuple<k-mer Number, k-mer total coverage> >
    unordered_map<uint16_t, tuple<uint64_t, uint64_t> > hapIdxKmerInfoMap;

    for (const auto& GraphKmerHashHapStrMapIter : (*nodePtr).GraphKmerHashHapStrMapIterVec) {

        const vector<int8_t>& hapBitTmp = GraphKmerHashHapStrMapIter->second.BitVec;

        vector<bitset<8> > bitsVector;

        for (const auto& hapBit : hapBitTmp) {
            bitsVector.push_back(hapBit);
        }

        for (const auto& hapIdx : topHapVec) {  // vector<hapIdx>
            uint16_t quotient = DIVIDE_BY_8(hapIdx);  // Quotient
            uint16_t remainder = GET_LOW_3_BITS(hapIdx);  // get remainder

            auto it = hapIdxKmerInfoMap.find(hapIdx);

            if (bitsVector[quotient][remainder]) {
                if (it == hapIdxKmerInfoMap.end()) {
                    hapIdxKmerInfoMap[hapIdx] = {1, GraphKmerHashHapStrMapIter->second.c};
                } else {
                    ++get<0>(it->second);
                    get<1>(it->second) += GraphKmerHashHapStrMapIter->second.c;
                }
            } else if (it == hapIdxKmerInfoMap.end()) {
                hapIdxKmerInfoMap[hapIdx] = {0, 0};
            }
        }
    }


    /* ------------------------------------------------- Posterior probability ------------------------------------------------- */
    // denominator
    long double denominator = 0.0L;

    for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
        denominator += HMMScoreTmp.a * HMMScoreTmp.b;

        if (debugGenotype) {
            cerr << "start:" << nodeStart << " haplotype:";
            for (const auto& hap : HMMScoreTmp.hapVec) {
                cerr << +hap << " ";
            }
            cerr << "genotype:";
            for (const auto& hap : HMMScoreTmp.hapVec) {
                cerr << +hapGtVec[hap] << " ";
            }
            cerr <<"alpha:" << HMMScoreTmp.a << " beta:" << HMMScoreTmp.b << " alpha*beta:" << HMMScoreTmp.a * HMMScoreTmp.b << endl;
        }
    }

    if (debugGenotype) {
        cerr << "denominator: " << denominator << endl;
    }

    // --poly-mix biallelic node: posterior over the true dosage d.
    // a*b carries the node emission twice (both alpha and beta include it), so gamma(s) = a*b/M[a_s]
    // and P(d | s, data) = P(d|a_s) L[d] / M[a_s]  ->  P(d) proportional to sum_s a*b P(d|a_s) L[d] / M[a_s]^2.
    if (polyMixGenotype && !debugGenotype) {
        auto mit = mixNodeCache.find(nodePtr);
        if (mit != mixNodeCache.end() && !(*nodePtr).HMMScoreVec.empty()) {
            const MixNode& mn = mit->second; const int P = mn.P;
            vector<long double> pd(P + 1, 0.0L);
            for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
                int a = 0;
                for (auto x : HMMScoreTmp.hapVec) a += (x < hapGtVec.size() && hapGtVec[x] == 1);
                if (a > P || mn.M[a] <= 0.0L) continue;
                long double w = HMMScoreTmp.a * HMMScoreTmp.b / (mn.M[a] * mn.M[a]);
                for (int d = 0; d <= P; ++d) pd[d] += w * mn.PdA[a][d] * mn.L[d];
            }
            long double tot = 0.0L; for (auto v : pd) tot += v;
            int best = 0;
            for (int d = 1; d <= P; ++d) if (pd[d] > pd[best]) best = d;
            // haplotypes carrying alt / ref here (prefer the selected ones, which have k-mer info)
            int altHap = -1, refHap = -1;
            for (auto h : topHapVec) {
                if (h >= hapGtVec.size()) continue;
                if (hapGtVec[h] == 1 && altHap < 0) altHap = h;
                if (hapGtVec[h] == 0 && refHap < 0) refHap = h;
            }
            for (size_t h = 0; h < hapGtVec.size() && (altHap < 0 || refHap < 0); ++h) {
                if (hapGtVec[h] == 1 && altHap < 0) altHap = (int)h;
                if (hapGtVec[h] == 0 && refHap < 0) refHap = (int)h;
            }
            if (tot > 0.0L && altHap >= 0 && refHap >= 0) {
                auto& posteriorInfo = (*nodePtr).posteriorInfo;
                std::lock_guard<std::mutex> mtx_locker(mtxG);
                posteriorInfo.probability = pd[best] / tot;
                posteriorInfo.hapVec.assign(P - best, (uint16_t)refHap);
                posteriorInfo.hapVec.insert(posteriorInfo.hapVec.end(), best, (uint16_t)altHap);
                posteriorInfo.kmerNumVec.clear(); posteriorInfo.kmerAveCovVec.clear();
                for (const auto& hapIdx : posteriorInfo.hapVec) {
                    uint64_t kmerNum = get<0>(hapIdxKmerInfoMap[hapIdx]);
                    float kmerAveCov = (kmerNum != 0) ? static_cast<float>(get<1>(hapIdxKmerInfoMap[hapIdx])) / (float)kmerNum : 0.0;
                    posteriorInfo.kmerNumVec.push_back(kmerNum);
                    posteriorInfo.kmerAveCovVec.push_back(kmerAveCov);
                }
                posteriorInfo.uniqueKmerNum = uniqueKmerNum;
                mixNodeCache.erase(mit);
                vector<HMMScore>().swap((*nodePtr).HMMScoreVec);
                return 0;
            }
        }
        if (mit != mixNodeCache.end()) mixNodeCache.erase(mit);
    }

    // --poly-mix fast path: integer genotype keys (sorted allele indices < 256, ploidy <= 8) instead of strings
    bool smallGt = polyMixGenotype && !debugGenotype;
    for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
        if (!smallGt) break;
        if (HMMScoreTmp.hapVec.size() > 8) { smallGt = false; break; }
        for (const auto& hap : HMMScoreTmp.hapVec) if (hapGtVec[hap] > 255) { smallGt = false; break; }
    }
    if (smallGt && denominator > 0.0L) {
        auto gtKey = [&](const vector<uint16_t>& hv) {
            uint16_t g[8]; size_t n = hv.size();
            for (size_t i = 0; i < n; ++i) g[i] = hapGtVec[hv[i]];
            std::sort(g, g + n);
            uint64_t key = n;
            for (size_t i = 0; i < n; ++i) key = (key << 8) | g[i];
            return key;   // all states at a node share n, so the packed sorted alleles identify the genotype
        };
        std::unordered_map<uint64_t, long double> gtSum;
        vector<uint64_t> keys; keys.reserve((*nodePtr).HMMScoreVec.size());
        for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
            if (HMMScoreTmp.hapVec.empty()) { keys.push_back(UINT64_MAX); continue; }
            uint64_t key = gtKey(HMMScoreTmp.hapVec);
            keys.push_back(key);
            gtSum[key] += (HMMScoreTmp.a * HMMScoreTmp.b) / denominator;
        }
        uint64_t bestKey = 0; long double bestSum = -1.0L;
        for (const auto& [key, s] : gtSum) if (s > bestSum) { bestSum = s; bestKey = key; }
        auto& posteriorInfo = (*nodePtr).posteriorInfo;
        long double maxPosteriorScore = 0.0L; size_t bestIdx = SIZE_MAX;
        for (size_t i = 0; i < (*nodePtr).HMMScoreVec.size(); ++i) {
            if (keys[i] != bestKey) continue;
            long double ps = ((*nodePtr).HMMScoreVec[i].a * (*nodePtr).HMMScoreVec[i].b) / denominator;
            if (bestIdx == SIZE_MAX || maxPosteriorScore < ps) { maxPosteriorScore = ps; bestIdx = i; }
        }
        std::lock_guard<std::mutex> mtx_locker(mtxG);
        posteriorInfo.probability = bestSum;
        if (bestIdx != SIZE_MAX) {
            posteriorInfo.hapVec = (*nodePtr).HMMScoreVec[bestIdx].hapVec;
            posteriorInfo.kmerNumVec.clear(); posteriorInfo.kmerAveCovVec.clear();
            for (const auto& hapIdx : posteriorInfo.hapVec) {
                uint64_t kmerNum = get<0>(hapIdxKmerInfoMap[hapIdx]);
                float kmerAveCov = (kmerNum != 0) ? static_cast<float>(get<1>(hapIdxKmerInfoMap[hapIdx])) / (float)kmerNum : 0.0;
                posteriorInfo.kmerNumVec.push_back(kmerNum);
                posteriorInfo.kmerAveCovVec.push_back(kmerAveCov);
            }
            posteriorInfo.uniqueKmerNum = uniqueKmerNum;
        }
        vector<HMMScore>().swap((*nodePtr).HMMScoreVec);
        return 0;
    }

    // Record the sum probability of each genotype
    map<string, long double> genotypeProbabilitySumMap;  // map<genotype, probability>

    // First pass: calculate probability sums for each genotype
    for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
        long double posteriorScore = (HMMScoreTmp.a * HMMScoreTmp.b)/(long double)denominator;

        // Genotype information
        vector<string> genotypeVec;
        for (const auto& hap : HMMScoreTmp.hapVec) {
            genotypeVec.push_back(to_string(hapGtVec[hap]));
        }
        if (genotypeVec.empty()) continue;  // Skip empty genotypes
        std::sort(genotypeVec.begin(), genotypeVec.end());  // sort
        string genotypeStr = join(genotypeVec, "/");

        // Record the probability of each genotype
        genotypeProbabilitySumMap[genotypeStr] += posteriorScore;
    }

    // Find the genotype with the highest score
    string highestSumScoreGenotype;
    long double highestScore = -1.0;

    for (const auto& entry : genotypeProbabilitySumMap) {
        if (entry.second > highestScore) {
            highestScore = entry.second;
            highestSumScoreGenotype = entry.first;
        }
    }

    // Second pass: find the haplotype information for the highest scoring genotype
    auto& posteriorInfo = (*nodePtr).posteriorInfo;  // Posterior Probabilities
    long double maxPosteriorScore = 0.0L;  // Record the maximum posterior probability (highestSumScoreGenotype)

    for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
        long double posteriorScore = (HMMScoreTmp.a * HMMScoreTmp.b)/(long double)denominator;

        // Genotype information
        vector<string> genotypeVec;
        for (const auto& hap : HMMScoreTmp.hapVec) {
            genotypeVec.push_back(to_string(hapGtVec[hap]));
        }
        if (genotypeVec.empty()) continue;  // Skip empty genotypes
        std::sort(genotypeVec.begin(), genotypeVec.end());  // sort
        string genotypeStr = join(genotypeVec, "/");

        if (genotypeStr != highestSumScoreGenotype) continue;  // Only the highestSumScoreGenotype is retained

        // Multithreaded Data Lock
        std::lock_guard<std::mutex> mtx_locker(mtxG);

        posteriorInfo.probability = highestScore;  // Posterior probability

        if (maxPosteriorScore < posteriorScore) {
            maxPosteriorScore = posteriorScore;  //  maximum posterior probability

            posteriorInfo.hapVec = HMMScoreTmp.hapVec;  // haplotype information

            posteriorInfo.kmerNumVec.clear();  // The number of k-mers corresponding to the haplotype
            posteriorInfo.kmerAveCovVec.clear();  // // The k-mer average depth corresponding to the haplotype

            for (const auto& hapIdx : posteriorInfo.hapVec) {
                uint64_t kmerNum = get<0>(hapIdxKmerInfoMap[hapIdx]);
                float kmerAveCov = (kmerNum != 0) ? static_cast<float>(get<1>(hapIdxKmerInfoMap[hapIdx])) / (float)kmerNum : 0.0;
                posteriorInfo.kmerNumVec.push_back(kmerNum);
                posteriorInfo.kmerAveCovVec.push_back(kmerAveCov);
            }

            // Number of node-unique k-mers, UK
            posteriorInfo.uniqueKmerNum = uniqueKmerNum;
        }
    }

    // Multithreaded Data Lock
    std::lock_guard<std::mutex> mtx_locker(mtxG);

    // Memory Deallocation
    vector<HMMScore>().swap((*nodePtr).HMMScoreVec);

    return 0;
}


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
uint8_t GENOTYPE::get_UK(
    const vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator>& GraphKmerHashHapStrMapIterVec
) {
    uint8_t uniqueKmerNum = 0;
    
    for (const auto& GraphKmerHashHapStrMapIter : GraphKmerHashHapStrMapIterVec) {  // use structured bindings to access only the key
        if (GraphKmerHashHapStrMapIter->second.f > 1) continue;
        
        if (uniqueKmerNum < UINT8_MAX) uniqueKmerNum++;
    }
    return uniqueKmerNum;
}


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
double GENOTYPE::cal_phred_scaled(long double value) {
    return (value >= 1.0) ? 99 : (-10 * log10(1.0 - value));
}


/**
 * @author zezhen du
 * @date 2023/07/12
 * @version v1.0
 * @brief save result
 * 
 * @param GraphMap            construct_index output result, map<string, map<uint32_t, nodeSrt>>
 * @param vcfHead             the VCF file comment lines
 * @param vcfInfoMap          vcf information, for output
 * @param sampleName          sample name
 * @param outputFileName      output file information
 * @param minSupportingGQ     minimum site quality (GQ) value for genotype
 * 
 * @return 0
**/
int GENOTYPE::save(
    map<string, map<uint32_t, nodeSrt> > & GraphMap, 
    string vcfHead, 
    map<string, map<uint32_t, vector<string> > > & vcfInfoMap, 
    const string& sampleName, 
    const string & outputFileName, 
    const float& minSupportingGQ
) {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Wrote genotyped variants to '" << outputFileName << "' \n\n";

    std::stringstream oss;
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10;  // Cache size is 10mb
    oss.str().reserve(CACHE_SIZE);
    // Format the output string with two decimal places
    oss << std::fixed << std::setprecision(1);

    // save the VCF file comment lines
    oss << vcfHead + "\t" + sampleName + "\n";

    SAVE SAVEClass(outputFileName);

    for (const auto& [chromosome, startVcfInfoMap] : vcfInfoMap) {
        // Check if the chromosome is in the posterior probability table
        auto findIter1 = GraphMap.find(chromosome);
        if (findIter1 == GraphMap.end()) {
            continue;
        }

        for (const auto& [nodeStart, vcfInfo] : startVcfInfoMap) {
            // Check if the node is in the posterior probability table
            auto findIter2 = findIter1->second.find(nodeStart);
            if (findIter2 == findIter1->second.end()) {
                continue;
            }

            // Variable binding
            const auto& hapGtVec = findIter2->second.hapGtVec;

            // Check if it is in the posterior probability table, if yes, add the corresponding genotype, otherwise ./.
            const auto& posteriorInfo = findIter2->second.posteriorInfo;

            if (posteriorInfo.hapVec.size() > 0) {
                // genotyping result
                vector<string> gtTxtVec;
                for (size_t i = 0; i < posteriorInfo.hapVec.size(); i++) {
                    gtTxtVec.push_back(to_string(hapGtVec[posteriorInfo.hapVec[i]]));
                }

                // If the genotype is empty, skip this site
                if (all_of(gtTxtVec.begin(), gtTxtVec.end(), [](string c) { return c == "0" || c == "."; })) {
                    continue;
                }

                // Original information
                for (size_t i = 0; i < 9; i++) {
                    if (i == 0) {
                        oss << vcfInfo[i];
                    } else if (i < 8) {
                        if (i == 6) {  // FILTER
                            oss << "\tPASS";
                        } else {
                            oss << "\t" << vcfInfo[i];
                        }
                    } else if (i == 8) {
                        oss << "\t" << "GT:GQ:GPP:NAK:CAK:UK";
                    }
                }

                // Genotype quality (phred-scaled 1 - max(GPP)) (GQ)
                float GQ = cal_phred_scaled(posteriorInfo.probability);

                // missing site
                if (GQ < minSupportingGQ) {
                    std::fill(gtTxtVec.begin(), gtTxtVec.end(), ".");
                }

                // Genotype
                oss << "\t" << join(gtTxtVec, "/");
                
                // Genotype quality , Genotype posterior probabilities, Number of allele k-mers
                oss << ":" << GQ << ":" 
                    << posteriorInfo.probability << ":" 
                    << join(posteriorInfo.kmerNumVec, ",") << ":" ;
                
                // Coverage of allele k-mers
                for (size_t i = 0; i < posteriorInfo.kmerAveCovVec.size(); i++) {
                    if (i == 0) {
                        oss << posteriorInfo.kmerAveCovVec[i];
                    } else {
                        oss << "," << posteriorInfo.kmerAveCovVec[i];
                    }
                }
                
                // Number of node-unique k-mers, UK
                oss << ":" << +posteriorInfo.uniqueKmerNum << "\n";
            }
            
            // Check the buffer size, if it exceeds the maximum, write the buffer content to the file and empty the buffer
            if (oss.tellp() > CACHE_SIZE) {
                string outTxt = oss.str();
                SAVEClass.save(outTxt);
                // empty oss
                oss.str("");
                oss.clear();
            }
        }
    }

    if (oss.tellp() > 0) {
        string outTxt = oss.str();
        SAVEClass.save(outTxt);
        // empty oss
        oss.str("");
        oss.clear();
    }

    return 0;
}