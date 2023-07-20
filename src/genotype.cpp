// g++ -c src/genotype.cpp -std=c++17 -O3 -march=native

#include "../include/genotype.hpp"
#include "../include/construct_index.hpp"

using namespace std;

// global variable
bool debugGenotype = false;

std::mutex mtxG;


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
int GENOTYPE::genotype(
    const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap, 
    map<string, map<uint32_t, nodeSrt> > & GraphMap, 
    const map<uint16_t, string> & hapMap, 
    const string& vcfHead, 
    map<string, map<uint32_t, string> > & vcfInfoMap, 
    const string & outputFileName, 
    const uint32_t & threads, 
    const bool & debug
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Genotyping ...\n";  // print log

    // debug
    debugGenotype = debug;

    // 进程池
    ThreadPool pool(threads);

    // 保存多线程的结果
    vector<future<int> > futureOutVec;

    // 初始化线程池
    pool.init();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Applying forward and backward algorithm ...\n";
    for (auto iter1 = GraphMap.begin(); iter1 != GraphMap.end(); iter1++) {  // map<string, map<uint32_t, nodeSrt> >: map<chr, map<nodeStart, nodeSrt> >
        // 每个线程处理的元素数量
        uint32_t elementsPerThread = iter1->second.size() / threads;
        for (int i = 0; i < threads; ++i) {
            uint32_t threadStart = i * elementsPerThread;
            uint32_t threadEnd = (i == threads - 1) ? iter1->second.size() - 1 : (threadStart + elementsPerThread - 1);

            // Submit and save results in multiple threads
            futureOutVec.push_back(
                pool.submit(
                    for_bac_post_run, 
                    ref(GraphKmerCovFreMap), 
                    iter1,  
                    ref(hapMap), 
                    threadStart, 
                    threadEnd
                )
            );
        }
    }

    // To obtain the return value of a function
    if (futureOutVec.size() >= 0) {
        // Save the result of the forward-backward algorithm into the graph
        for (size_t i = 0; i < futureOutVec.size(); i++) {
            int result = move(futureOutVec[i].get());
        }

        // empty vector
        futureOutVec.clear();
    }

    malloc_trim(0); // 0 is for heap memory

    // 关闭线程池
    pool.shutdown();

    cerr << endl << endl;

    // save
    GENOTYPE::save(GraphMap, vcfHead, vcfInfoMap, outputFileName);

    return 0;
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0
 * @brief calculate Node k-mers Coverage
 * 
 * @param GraphKmerCovFreMap map<kmerHash, coverage>
 * @param kmerHashHapVec      vector<kmerHashHap>   
 * @param nodeAveKmerCoverage the average k-mers coverage of node
 * @param kmerCoverageVec     the k-mers vector of node
 * 
 * @return void
**/
void GENOTYPE::cal_node_cov(
    const std::unordered_map<uint64_t, kmerCovFre> &GraphKmerCovFreMap, 
    const vector<kmerHashHap> & kmerHashHapVec, 
    float & nodeAveKmerCoverage,
    vector<uint8_t> & kmerCoverageVec
)
{
    kmerCoverageVec.reserve(kmerHashHapVec.size());

    uint64_t nodeKmerNum = 0;
    uint64_t nodeSumKmerCoverage = 0;
    
    for (const auto& [kmerHash, _] : kmerHashHapVec) {  // use structured bindings to access only the key
        ++nodeKmerNum;

        if (auto findIter = GraphKmerCovFreMap.find(kmerHash); findIter != GraphKmerCovFreMap.end()) {
            nodeSumKmerCoverage += findIter->second.c;
            kmerCoverageVec.emplace_back(findIter->second.c);
        } else {
            kmerCoverageVec.emplace_back(0);
        }
    }

    nodeAveKmerCoverage = (nodeKmerNum > 0) ? static_cast<float>(nodeSumKmerCoverage) / nodeKmerNum : 0.0f;
}

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
int GENOTYPE::for_bac_post_run(
    const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap, 
    map<string, map<uint32_t, nodeSrt> >::iterator startNodeIter, 
    const map<uint16_t, string> & hapMap, 
    uint32_t threadStart, 
    uint32_t threadEnd
)
{
    string chromosome = startNodeIter->first;  // chromosome

    // Get the forward/backward iterator
    auto newStartNodeIterL = std::next(startNodeIter->second.begin(), threadStart);
    auto newStartNodeIterR = std::next(startNodeIter->second.begin(), threadEnd + 1);
    
    // Check if the iterator is out of range
    if (newStartNodeIterL == startNodeIter->second.end())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Warning: the start iterator is pointing to the end. Please check the code -> " << chromosome << endl;
        return 0;
    }

    /* ************************************************** Haplotype selection ************************************************** */
    vector<uint16_t> topHapVec;
    haplotype_selection(
        GraphKmerCovFreMap, 
        chromosome, 
        newStartNodeIterL, 
        newStartNodeIterR, 
        hapMap, 
        topHapVec
    );

    /* ************************************************** forward ************************************************** */
    vector<HMMScore> preHMMScoreVec;  // Result of previous node forward/backward algorithm
    uint32_t preNodeEnd = 0;  // The ending position of the previous node
    // forward
    for (map<uint32_t, nodeSrt>::iterator iter1 = newStartNodeIterL; iter1 != newStartNodeIterR; iter1++) {  // map<nodeStart, nodeSrt>
        uint32_t nodeStart = iter1->first;  // the starting position of the node
        uint32_t nodeEnd = nodeStart + iter1->second.seqMap[0].size() - 1;  // end position of the node

        if (iter1->second.hapGtVec.size() == 1) {continue;}  // If only 0, skip the node

        // The Average coverage and vector of node k-mers
        float nodeAveKmerCoverage;
        vector<uint8_t> kmerCoverageVec;
        cal_node_cov(
            GraphKmerCovFreMap, 
            iter1->second.kmerHashHapVec, 
            nodeAveKmerCoverage, 
            kmerCoverageVec
        );

        // hidden state
        map<uint16_t, map<uint16_t, vector<uint8_t> > > hiddenStatesMap = GENOTYPE::hidden_states(
            iter1->second, 
            topHapVec
        );

        if (debugGenotype) {
            // for (const auto& KmerHashHap : iter1->second.kmerHashHapVec) {
            //     cerr << KmerHashHap.kmerHash << endl;
            // }
            
            for (auto itTmp1 : hiddenStatesMap) {
                for (auto itTmp2 : itTmp1.second) {
                    cerr << "start:" << nodeStart << " genotype:" << +itTmp1.first << " genotype:" << +itTmp2.first << " hiddenStates/coverage:\n";
                    for (auto itTmp3 : itTmp2.second) {
                        cerr << setw(4) << +itTmp3;
                    }
                    cerr << endl;

                    // coverage
                    for (auto itTmp1 : kmerCoverageVec) {
                        cerr << setw(4) << +itTmp1;
                    }
                    cerr << endl;
                }
            }
        }

        // transition probability
        uint32_t nodeDistence = nodeStart - preNodeEnd;  // The distance from the previous node
        long double recombProb, noRecombProb;  // Probability of recombination and non-recombination
        tie(recombProb, noRecombProb) = transition_probabilities(
            nodeDistence, 
            hapMap
        );

        // observation matrix
        map<uint16_t, map<uint16_t, long double> > observableStatesMap = GENOTYPE::observable_states(
            nodeAveKmerCoverage, 
            hiddenStatesMap, 
            kmerCoverageVec, 
            iter1->second
        );

        if (debugGenotype) {
            for (auto itTmp1 : observableStatesMap) {
                for (auto itTmp2 : itTmp1.second) {
                    cerr << "start:" << nodeStart << " genotype:" << +itTmp1.first << " genotype:" << +itTmp2.first << " observableStates:" << itTmp2.second << endl;
                }
            }
        }


        // Variable binding, used to modify the score of the forward algorithm in the node corresponding to nodeStart
        auto& HMMScoreVecTmp = startNodeIter->second.at(nodeStart).HMMScoreVec;
        forward(
            preHMMScoreVec, 
            recombProb, 
            noRecombProb, 
            observableStatesMap, 
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
                cerr << " Alpha:" << HMMScoreTmp.a << endl;
            }
        }

        // Reset the information of the previous node
        preNodeEnd = nodeEnd;
        preHMMScoreVec = HMMScoreVecTmp;
    }


    /* ************************************************** backward ************************************************** */
    // variable zeroing
    vector<HMMScore>().swap(preHMMScoreVec);  // Result of previous node forward/backward algorithm
    preNodeEnd = 0;  // The ending position of the previous node
    for (std::map<uint32_t, nodeSrt>::reverse_iterator iter1 = std::map<uint32_t, nodeSrt>::reverse_iterator(newStartNodeIterR); iter1 != std::map<uint32_t, nodeSrt>::reverse_iterator(newStartNodeIterL); ++iter1) {
        uint32_t nodeStart = iter1->first;  // the starting position of the node
        uint32_t nodeEnd = nodeStart + iter1->second.seqMap[0].size() - 1;  // end position of the node

        if (iter1->second.hapGtVec.size() == 1) {continue;}  // If only 0, skip the node

        // The Average coverage of node kmer
        float nodeAveKmerCoverage;
        vector<uint8_t> kmerCoverageVec;
        cal_node_cov(
            GraphKmerCovFreMap, 
            iter1->second.kmerHashHapVec, 
            nodeAveKmerCoverage, 
            kmerCoverageVec
        );


        // hidden state
        map<uint16_t, map<uint16_t, vector<uint8_t> > > hiddenStatesMap = GENOTYPE::hidden_states(
            iter1->second, 
            topHapVec
        );


        // transition probability
        uint32_t nodeDistence = preNodeEnd - nodeStart;  // The distance from the previous node
        long double recombProb, noRecombProb;  // Probability of recombination and non-recombination
        tie(recombProb, noRecombProb) = transition_probabilities(
            nodeDistence, 
            hapMap
        );


        // observation matrix
        map<uint16_t, map<uint16_t, long double> > observableStatesMap = GENOTYPE::observable_states(
            nodeAveKmerCoverage, 
            hiddenStatesMap, 
            kmerCoverageVec, 
            iter1->second
        );


        // Variable binding, used to modify the score of the backward algorithm in the node corresponding to nodeStart
        auto& HMMScoreVecTmp = startNodeIter->second.at(nodeStart).HMMScoreVec;
        backward(
            preHMMScoreVec, 
            recombProb, 
            noRecombProb, 
            observableStatesMap, 
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
                cerr << " Beta:" << HMMScoreTmp.a << endl;
            }
        }

        // Reset the information of the previous node
        preNodeEnd = nodeEnd;
        preHMMScoreVec = HMMScoreVecTmp;
    }

    /* ************************************************** posteriorTup ************************************************** */
    for (map<uint32_t, nodeSrt>::iterator iter1 = newStartNodeIterL; iter1 != newStartNodeIterR; iter1++) {  // map<nodeStart, nodeSrt>
        if (iter1->second.hapGtVec.size() == 1) {continue;}  // If only 0, skip the node

        // Calculate and store the result in a local variable
        posterior(iter1->first, &(iter1->second));
    }

    return 0;
}


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
void GENOTYPE::haplotype_selection(
    const std::unordered_map<uint64_t, kmerCovFre> &GraphKmerCovFreMap, 
    const string& chromosome, 
    const map<uint32_t, nodeSrt>::iterator& newStartNodeIterL, 
    const map<uint32_t, nodeSrt>::iterator& newStartNodeIterR, 
    const map<uint16_t, string> & hapMap, 
    vector<uint16_t>& topHapVec
)
{
    // If the number of haplotypes is less than 20, skip directly
    if (hapMap.size() < 20) {
        for (const auto& pair : hapMap) {
            topHapVec.push_back(pair.first);
        }
        return;
    }
    
    // k-mer count of all haplotypes
    vector<uint32_t> hapKmerCountVec(hapMap.size(), 0);

    for (map<uint32_t, nodeSrt>::iterator nodeIter = newStartNodeIterL; nodeIter != newStartNodeIterR; nodeIter++) {  // map<nodeStart, nodeSrt>

        if (nodeIter->second.hapGtVec.size() == 1) {continue;}  // If only 0, skip the node

        // Traverse all nodes kmer
        for (const auto& [kmerHash, hapBitTmp] : nodeIter->second.kmerHashHapVec) {
            // Check the coverage of k-mer
            uint8_t kmerCov = 0;
            if (auto findIter = GraphKmerCovFreMap.find(kmerHash); findIter != GraphKmerCovFreMap.end()) {
                kmerCov = findIter->second.c;
            }
            // if it is 0, skip the k-mer
            if (kmerCov == 0) {continue;}
            
            // haplotype information
            vector<bitset<8> > bitsVector;
            for (const auto& hapBit : hapBitTmp) {
                bitsVector.push_back(hapBit);
            }
            
            // haplotype1
            for (const auto& [hapIdx, _] : hapMap) {  // map<hapIdx, hapName>
                uint16_t quotient1 = DIVIDE_BY_8(hapIdx);  // acquirer
                uint16_t remainder1 = GET_LOW_3_BITS(hapIdx); // get remainder
                // If the haplotype contains the k-mer then add the frequency of the k-mer
                if (bitsVector[quotient1][remainder1]) {
                    hapKmerCountVec[hapIdx]++;
                }
            }
        }
    }

    // Calculation of haplotype probabilities from the Dirichlet distribution
    HaplotypeSelect HaplotypeSelectClass(hapKmerCountVec);
    HaplotypeSelectClass.calculate_sparsity();
    HaplotypeSelectClass.generate_sparse_frequency_vector();
    HaplotypeSelectClass.get_top_indices(15);
    topHapVec = HaplotypeSelectClass.mTopHapVec;

    std::lock_guard<std::mutex> mtx_locker(mtxG);
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Haplotype selection results for " << chromosome << "-" << newStartNodeIterL->first << ":";  // print log
    for (uint16_t i = 0; i < topHapVec.size(); i++) {
        if (i == 0) {
            cerr << " " << topHapVec[i];
        } else {
            cerr << ", " << topHapVec[i];
        }
    }
    cerr << endl;
}


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
map<uint16_t, map<uint16_t, vector<uint8_t> > > GENOTYPE::hidden_states(
    const nodeSrt & node, 
    const vector<uint16_t>& topHapVec
)
{
    // 所有的单倍型组合
    map<uint16_t, map<uint16_t, vector<uint8_t> > > hiddenStatesMap;  // map<sample1, map<sample2, vector<k-mers coverage> > >

    const vector<kmerHashHap> & kmerHashHapVecTmp = node.kmerHashHapVec;  // 该节点所有的kmer信息

    // 遍历所有的节点kmer，构造隐藏状态
    for (const auto& [kmerHash, hapBitTmp] : kmerHashHapVecTmp) {

        vector<bitset<8> > bitsVector;

        for (const auto& hapBit : hapBitTmp) {
            bitsVector.push_back(hapBit);
        }
        
        // haplotype1
        for (const auto& hapIdx1 : topHapVec) {  // vector<hapIdx>

            uint16_t quotient1 = DIVIDE_BY_8(hapIdx1);  // acquirer
            uint16_t remainder1 = GET_LOW_3_BITS(hapIdx1); // get remainder
            // haplotype2
            for (const auto& hapIdx2 : topHapVec) {  // vector<hapIdx>

                if (hapIdx2 < hapIdx1) {continue;}  // If haplotype 2 is smaller than haplotype 1, skip this combination

                uint16_t quotient2 = DIVIDE_BY_8(hapIdx2);  // acquirer
                uint16_t remainder2 = GET_LOW_3_BITS(hapIdx2); // get remainder

                hiddenStatesMap[hapIdx1][hapIdx2].emplace_back(bitsVector[quotient1][remainder1] + bitsVector[quotient2][remainder2]);
            }
        }
    }

    return hiddenStatesMap;
}


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
pair<long double, long double> GENOTYPE::transition_probabilities(const uint32_t & nodeDistence, const map<uint16_t, string> & hapMap)
{
    // 转移概率
    double effectivePopulationSize = 1e-05;  // 有效种群大小
    uint16_t populationSize = hapMap.size();  // 群体大小
    double recombRate = 1.26;  // 重组率
    long double distance = nodeDistence * 0.000004L * ((long double) recombRate) * effectivePopulationSize;  // d
    long double recombProb = (1.0L - exp(-distance / (long double) populationSize) )* (1.0L / (long double) populationSize);
    long double noRecombProb = exp(-distance / (long double) populationSize) + recombProb;
    return make_pair(recombProb, noRecombProb);  // 返回重组和不重组的概率
}


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
map<uint16_t, map<uint16_t, long double> > GENOTYPE::observable_states(
    const float & aveKmerCoverage, 
    const map<uint16_t, map<uint16_t, vector<uint8_t> > > & hiddenStatesMap, 
    const vector<uint8_t> & kmerCoverageVec, 
    nodeSrt & node
)
{
    map<uint16_t, map<uint16_t, long double> > observableStatesMap;
    
    for (auto& it1 : hiddenStatesMap) {  // map<uint16_t, map<uint16_t, vector<uint8_t> > >
        for (auto& it2 : it1.second) {  // map<uint16_t, vector<uint8_t> >
            long double result = 1.0L;  // save result

            // Check if the hidden state is all o, if yes, assign 0
            bool allZero = std::all_of(it2.second.begin(), it2.second.end(), [](uint8_t value) {
                return value == 0;
            });

            if (allZero) {
                observableStatesMap[it1.first][it2.first] = 0.0;
                continue;
            }
            
            for (size_t i = 0; i < it2.second.size(); i++) {  // vector<uint8_t>
                uint8_t kmerCoverage = kmerCoverageVec[i];  // kmer coverage
                
                uint8_t kmerCount = it2.second[i];  // value in hidden state

                if (kmerCount == 2) {
                    if (debugGenotype) {
                        cerr << +it1.first << " " << +it2.first << ": " << +kmerCount << " " << aveKmerCoverage << " " << +kmerCoverage << " " << poisson(aveKmerCoverage, kmerCoverage) << " " << result << endl;
                    }
                    result *= poisson(aveKmerCoverage, kmerCoverage);
                } else if (kmerCount == 1) {
                    if (debugGenotype) {
                        cerr << +it1.first << " " << +it2.first << ": " << +kmerCount << " " << aveKmerCoverage << " " << +kmerCoverage << " " << poisson(aveKmerCoverage / 2.0, kmerCoverage) << " " << result << endl;
                    }
                    result *= poisson(aveKmerCoverage / 2.0, kmerCoverage);
                } else {
                    if (debugGenotype) {
                        cerr << +it1.first << " " << +it2.first << ": " << +kmerCount << " " << aveKmerCoverage << " " << +kmerCoverage << " " << geometric(get_error_param(aveKmerCoverage), kmerCoverage) << " " << result << endl;
                    }
                    result *= geometric(get_error_param(aveKmerCoverage), kmerCoverage);
                }
            }
            observableStatesMap[it1.first][it2.first] = result;
        }
    }
    return observableStatesMap;
}

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
long double GENOTYPE::poisson(
    long double mean, 
    unsigned int value
)
{
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
 * @param aveKmerCoverage         kmer的平均覆盖度
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
 * @brief geometric，隐藏状态为0时的返回状态
 * 
 * @param p         Correct rate, output by get_error_param
 * @param value     隐藏状态中的数值
 * 
 * @return geometric
**/
long double GENOTYPE::geometric(
    long double p, 
    unsigned int value
)
{
    return pow(1.0L - p, value)*p;
}


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
 * @param @param HMMScoreVec    Alpha of this state
 * 
 * @return void
**/
void GENOTYPE::forward(
    const vector<HMMScore> & preHMMScoreVec, 
    const long double & recombProb, 
    const long double & noRecombProb, 
    const map<uint16_t, map<uint16_t, long double> > & observableStatesMap, 
    vector<HMMScore>& HMMScoreVec
)
{
    vector<long double> alphaVecTmp;  // 临时保存
    long double alphaSumTmp = 0.0L;  // 临时保存

    for (const auto& [genotype1, observableStatesMapForGenotype1] : observableStatesMap) {
        for (const auto& [genotype2, observableStates] : observableStatesMapForGenotype1) {
            long double result = 0.0L;  // 保存结果

            if (preHMMScoreVec.empty()) {  // 如果是第一个节点
                result += observableStates;
            }
            else {  // 非第一个节点
                for (const auto& HMMScoreTmp : preHMMScoreVec) {  // vector<HMMScore>
                    // genotype
                    if (HMMScoreTmp.hapVec[0] != genotype1 && HMMScoreTmp.hapVec[0] != genotype2 && HMMScoreTmp.hapVec[1] != genotype1 && HMMScoreTmp.hapVec[1] != genotype2) {  // 都不一样
                        result += HMMScoreTmp.a * recombProb * recombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype1 && HMMScoreTmp.hapVec[1] != genotype2) {  // 有一个一样
                        result += HMMScoreTmp.a * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] != genotype1 && HMMScoreTmp.hapVec[1] == genotype2) {  // 有一个一样
                        result += HMMScoreTmp.a * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype2 && HMMScoreTmp.hapVec[1] != genotype1) {  // 有一个一样
                        result += HMMScoreTmp.a * recombProb * noRecombProb * observableStates;
                    }  else if (HMMScoreTmp.hapVec[0] != genotype2 && HMMScoreTmp.hapVec[1] == genotype1) {  // 有一个一样
                        result += HMMScoreTmp.a * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype1 && HMMScoreTmp.hapVec[1] == genotype2) {  // 都一样
                        result += HMMScoreTmp.a * noRecombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype2 && HMMScoreTmp.hapVec[1] == genotype1) {  // 都一样 
                        result += HMMScoreTmp.a * noRecombProb * noRecombProb * observableStates;
                    } else {
                        cerr << "[" << __func__ << "::" << getTime() << "] " 
                            << "Error: Unknown situation: " 
                            << +genotype1 << " " 
                            << +genotype2 << " " 
                            << +HMMScoreTmp.hapVec[0] << " " 
                            << +HMMScoreTmp.hapVec[1] << endl;
                        exit(1);
                    }
                }
            }
            alphaVecTmp.push_back(result);
            alphaSumTmp += result;
        }
    }
    
    // 标准化，防止数值过小归零
    if (alphaSumTmp > 0.0L) {
        transform(alphaVecTmp.begin(), alphaVecTmp.end(), alphaVecTmp.begin(), bind(divides<long double>(), placeholders::_1, alphaSumTmp));
    } else {
        long double uniform = 1.0L / (long double) alphaVecTmp.size();
        transform(alphaVecTmp.begin(), alphaVecTmp.end(), alphaVecTmp.begin(), [uniform](long double c) -> long double {return uniform;});
    }

    // index
    HMMScoreVec.resize(alphaVecTmp.size());
    uint32_t indexTmp = 0;
    for (const auto& [genotype1, observableStatesMapForGenotype1] : observableStatesMap) {  // map<uint16_t, map<uint16_t, long double> >
        for (const auto& [genotype2, _] : observableStatesMapForGenotype1) {  // map<uint16_t, long double>
            // Multithreaded Data Lock
            std::lock_guard<std::mutex> mtx_locker(mtxG);
            HMMScoreVec[indexTmp].a = alphaVecTmp[indexTmp];
            HMMScoreVec[indexTmp].hapVec = {genotype1, genotype2};
            indexTmp++;
        }
    }
}


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
void GENOTYPE::backward(
    const vector<HMMScore> & preHMMScoreVec, 
    const long double & recombProb, 
    const long double & noRecombProb, 
    const map<uint16_t, map<uint16_t, long double> > & observableStatesMap, 
    vector<HMMScore>& HMMScoreVec
)
{
    vector<long double> betaVecTmp;  // 临时保存
    long double betaSumTmp = 0.0L;  // 临时保存

    for (const auto& [genotype1, observableStatesMapForGenotype1] : observableStatesMap) {
        for (const auto& [genotype2, observableStates] : observableStatesMapForGenotype1) {
            long double result = 0.0L;  // 保存结果

            if (preHMMScoreVec.empty()) {  // 如果是第一个节点
                result += observableStates;
            }
            else {  // 非第一个节点
                for (const auto& HMMScoreTmp : preHMMScoreVec) {  // vector<HMMScore>
                    if (HMMScoreTmp.hapVec[0] != genotype1 && HMMScoreTmp.hapVec[0] != genotype2 && HMMScoreTmp.hapVec[1] != genotype1 && HMMScoreTmp.hapVec[1] != genotype2) {  // 都不一样
                        result += HMMScoreTmp.b * recombProb * recombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype1 && HMMScoreTmp.hapVec[1] != genotype2) {  // 有一个一样
                        result += HMMScoreTmp.b * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] != genotype1 && HMMScoreTmp.hapVec[1] == genotype2) {  // 有一个一样
                        result += HMMScoreTmp.b * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype2 && HMMScoreTmp.hapVec[1] != genotype1) {  // 有一个一样
                        result += HMMScoreTmp.b * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] != genotype2 && HMMScoreTmp.hapVec[1] == genotype1) {  // 有一个一样
                        result += HMMScoreTmp.b * recombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype2 && HMMScoreTmp.hapVec[1] == genotype1) {  // 都一样
                        result += HMMScoreTmp.b * noRecombProb * noRecombProb * observableStates;
                    } else if (HMMScoreTmp.hapVec[0] == genotype1 && HMMScoreTmp.hapVec[1] == genotype2) {  // 都一样
                        result += HMMScoreTmp.b * noRecombProb * noRecombProb * observableStates;
                    } else {
                        cerr << "[" << __func__ << "::" << getTime() << "] " 
                            << "Error: Unknown situation: " 
                            << +genotype1 << " " 
                            << +genotype2 << " " 
                            << +HMMScoreTmp.hapVec[0] << " " 
                            << +HMMScoreTmp.hapVec[1] << endl;
                        exit(1);
                    }
                }
            }
            betaVecTmp.push_back(result);
            betaSumTmp += result;
        }
    }

    // 标准化，防止数值过小归零
    if (betaSumTmp > 0.0L) {
        transform(betaVecTmp.begin(), betaVecTmp.end(), betaVecTmp.begin(), bind(divides<long double>(), placeholders::_1, betaSumTmp));
    } else {
        long double uniform = 1.0L / (long double) betaVecTmp.size();
        transform(betaVecTmp.begin(), betaVecTmp.end(), betaVecTmp.begin(), [uniform](long double c) -> long double {return uniform;});
    }

    // index
    size_t indexTmp = 0;
    for (const auto& [genotype1, observableStatesMapForGenotype1] : observableStatesMap)  // map<uint16_t, map<uint16_t, long double> >
    {
        for (const auto& [genotype2, _] : observableStatesMapForGenotype1)  // map<uint16_t, long double>
        {
            // Multithreaded Data Lock
            std::lock_guard<std::mutex> mtx_locker(mtxG);
            HMMScoreVec[indexTmp].b = betaVecTmp[indexTmp];
            indexTmp++;
        }
    }
}

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
int GENOTYPE::posterior(
    uint32_t nodeStart, 
    nodeSrt* nodePtr
)
{
    if (nodePtr == nullptr) {return 0;}

    // 分母
    long double denominator = 0.0L;

    for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
        denominator += HMMScoreTmp.a * HMMScoreTmp.b;

        if (debugGenotype) {
            cerr << "start:" << nodeStart << " ";
            for (const auto& hap : HMMScoreTmp.hapVec) {
                cerr << +hap << " ";
            }
            cerr << "alpha:" << HMMScoreTmp.a << " beta:" << HMMScoreTmp.b << " alpha*beta:" << HMMScoreTmp.a * HMMScoreTmp.b << endl;
        }
    }

    if (debugGenotype) {
        cerr << "denominator: " << denominator << endl;
    }

    for (const auto& HMMScoreTmp : (*nodePtr).HMMScoreVec) {
        long double posteriorScore = (HMMScoreTmp.a * HMMScoreTmp.b)/(long double)denominator;

        // Multithreaded Data Lock
        std::lock_guard<std::mutex> mtx_locker(mtxG);

        // Recording Posterior Probabilities
        auto& [currentScore, currentGenotypeVec] = (*nodePtr).posteriorTup;
        if (currentScore < posteriorScore) {
            currentScore = posteriorScore;
            currentGenotypeVec = HMMScoreTmp.hapVec;
        }
    }

    // Multithreaded Data Lock
    std::lock_guard<std::mutex> mtx_locker(mtxG);

    // Memory Deallocation
    (*nodePtr).HMMScoreVec = vector<HMMScore>();

    return 0;
}


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
int GENOTYPE::save(
    map<string, map<uint32_t, nodeSrt> > & GraphMap, 
    string vcfHead, 
    map<string, map<uint32_t, string> > & vcfInfoMap, 
    const string & outputFileName
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Wrote genotyped variants to '" << outputFileName << "' \n\n\n";

    SAVE SAVEClass(outputFileName);

    // save the VCF file comment lines
    SAVEClass.save(vcfHead);

    for (const auto& [chromosome, startVcfInfoMap] : vcfInfoMap) {
        for (const auto& [nodeStart, vcfInfo] : startVcfInfoMap) {
            string outTxt;
            // Check if it is in the posterior probability table, if yes, add the corresponding genotype, otherwise ./.
            vector<uint16_t> hapIdxVec = get<1>(GraphMap[chromosome][nodeStart].posteriorTup);
            if (hapIdxVec.size() > 0) {
                outTxt = vcfInfo + "\t" + 
                        to_string(+GraphMap[chromosome][nodeStart].hapGtVec[hapIdxVec[0]]) + "/" +
                        to_string(+GraphMap[chromosome][nodeStart].hapGtVec[hapIdxVec[1]]) + "\n";
            } else {
                outTxt = vcfInfo + "\t" + "./." + "\n";
            }
            SAVEClass.save(outTxt);
        }
    }

    return 0;
}