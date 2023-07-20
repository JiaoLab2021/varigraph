// g++ -c construct_index.cpp -std=c++17 -lz -O3 -march=native

#include "../include/kmer.hpp"
#include "../include/construct_index.hpp"

using namespace std;

std::mutex mtxCI;

// kseq.h
KSEQ_INIT(gzFile, gzread)


ConstructIndex::ConstructIndex(
    const string& refFileName, 
    const string& vcfFileName, 
    const string& inputMbfFileName, 
    const string& outputMbfFileName, 
    const uint32_t& kmerLen, 
    const string& prefix, 
    const uint32_t& ploidy, 
    const uint32_t& threads, 
    const bool& debug
) : refFileName_(refFileName), vcfFileName_(vcfFileName), inputMbfFileName_(inputMbfFileName), outputMbfFileName_(outputMbfFileName), 
    kmerLen_(kmerLen), prefix_(prefix), ploidy_(ploidy), threads_(threads), debug_(debug)
{
    mHapMap[0] = "reference";
}

ConstructIndex::~ConstructIndex()
{
    if (mbf != nullptr)
    {
        delete mbf;  // Release the memory occupied by the Counting Bloom filter
        mbf = nullptr;
    }
}

/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief Free memory
 * 
 * @return void
**/
void ConstructIndex::clear_memory()
{
    unordered_map<string, string>().swap(mFastaMap);
    if (mbf != nullptr) {
        delete mbf;
        mbf = nullptr;
    }
    malloc_trim(0);	// 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief Building the k-mers index of reference genome
 * 
 * @return void
**/
void ConstructIndex::build_fasta_index()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building refgenome index: " << refFileName_ << endl << endl << endl;

    // open fasta file
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    // ���ļ�
    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << refFileName_ << "': No such file or directory." << endl;
        exit(1);
    }
    else {
        kseq_t *ks;
        ks = kseq_init(gzfp);

        // Define local variables to avoid frequent application and release of memory inside the loop
        string chromosome;
        string sequence;

        while( kseq_read(ks) >= 0 ) {
            // ks->name.s  name
            // ks->seq.s   sequence
            chromosome = ks->name.s;
            sequence = ks->seq.s;

            // record the length of sequence
            genomeSize_ += ks->seq.l;

            // build fasta index
            mFastaMap.emplace(chromosome, sequence);
        }

        // free memory and close file
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    malloc_trim(0);	// 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/07/18
 * @version v1.0.1
 * @brief Making Counting Bloom Filter
 * 
 * @return void
**/
void ConstructIndex::make_mbf()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating k-mer frequency in reference genome ...\n";

    /* *************************************************** making or load *************************************************** */
    uint64_t bfSize = genomeSize_ - kmerLen_ + 1;
    double errorRate = 0.01;
    mbf = new BloomFilter(bfSize, errorRate);

    // making
    if (inputMbfFileName_.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Making Counting Bloom Filter with a false positive rate of " << errorRate << " ...\n";

        for (const auto& [chromosome, sequence] : mFastaMap)  // map<chromosome, sequence>
        {
            // Constructing k-mer index
            kmerBit::kmer_sketch_fasta(sequence, kmerLen_, mbf);

            cerr << "[" << __func__ << "::" << getTime() << "] " << "Successfully processed chromosome '" << chromosome << "' ...\n";
        }

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom Filter successfully constructed ..." << endl;

        // save to file
        if (!outputMbfFileName_.empty()) {
            mbf->save(outputMbfFileName_);
        }
        
        cerr << endl;
    } else {  // load from file
        mbf->load(inputMbfFileName_);

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom Filter successfully loaded ..." << endl << endl;
    }
    
    cerr << "           - " << "Size of Counting Bloom Filter: " << mbf->get_size() << endl;
    cerr << "           - " << "Number of hash functions: " << mbf->get_num() << endl;
    cerr << "           - " << "Usage rate of the Counting Bloom Filter: " << mbf->get_cap() << endl << endl << endl;

    malloc_trim(0);	// 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0
 * @brief Graph construction
 * 
 * @return void
**/
void ConstructIndex::construct ()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Constructing ...\n";  // print log

    // number
    uint32_t snpNum = 0;
    uint32_t indelNum = 0;
    uint32_t insNum = 0;
    uint32_t delNum = 0;
    uint32_t invNum = 0;
    uint32_t dupNum = 0;
    uint32_t otherNum = 0;

    // Record the position of the previous node
    uint32_t tmpRefEnd = 0;
    string tmpChromosome;

    // open file
    GzChunkReader GzChunkReaderClass(vcfFileName_);

    // read line
    string line;
    while (GzChunkReaderClass.read_line(line)) {
        // skip empty lines
        if (line.empty()) {
            continue;
        }

        line = strip(line, '\n');  // remove trailing newline

        // First judge whether it is a comment line, if yes, skip
        if (line.find("#") != string::npos && line.find("#CHROM") == string::npos) {
            mVcfHead += line + "\n";  // Store VCF file comment line
            continue;
        }

        // Split the line into tokens based on whitespace
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

        // Check if the file is correct, if not, jump out of the code
        if (lineVec.size() < 10)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "'" << vcfFileName_  << "': Error -> number of columns in the VCF file is less than 10.: " << lineVec.size() << endl;
            exit(1);
        }

        // store haplotype information
        if (line.find("#CHROM") != string::npos) {
            mVcfHead += line + "\t" + prefix_ + "\n";  // Store VCF file comment line

            uint16_t hapIdx = 1;  // haplotype index
            for (size_t i = 9; i < lineVec.size(); i++) {
                for (size_t j = 0; j < ploidy_; j++) {
                    mHapMap[hapIdx] = lineVec[i];  // assignment
                
                    if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                        hapIdx++;  // Haplotype index +1
                    } else {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: Haplotype number >= " << UINT16_MAX << endl;
                        exit(1);
                    }
                }
            }

            // haplotype number
            mHapNum = mHapMap.size();
        } else {  // ��ͼ
            string chromosome = lineVec[0];  // 0
            uint32_t refStart = stoul(lineVec[1]);  // 1
            string refSeq = lineVec[3];  // 3
            uint32_t refLen = refSeq.size();  // the length of reference sequence
            uint32_t refEnd = refStart + refSeq.size() - 1;  // end position
            string qrySeq = lineVec[4];  // 4
            vector<string> qrySeqVec = split(qrySeq, ",");  // the vector of ALT
            mVcfInfoMap[chromosome][refStart] = line;  // Store vcf file information for output
            
            // ��FORMAT�ֶ���gt��λ��
            vector<string> formatVec = split(strip(lineVec[8], '\n'), ":");
            vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

            int gtIndex = distance(formatVec.begin(), gtItera);
            if (gtIndex == formatVec.size()) {  // FORMAT��û��GT���˳�����
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: no GT information in FORMAT: " << line << endl;
                exit(1);
            }

            // Check if the chromosome is in the reference genome
            auto mFastaMapFindIter = mFastaMap.find(chromosome);
            if (mFastaMapFindIter == mFastaMap.end()) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "Error: no '" << chromosome << "' found in reference genome."<< endl;
                exit(1);
            }
            
            // ���vcf��refSeq�ͻ������Ƿ�һ��
            if (mFastaMapFindIter->second.substr(refStart-1, refSeq.length()) != refSeq) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "Error: base difference between reference and vcf: "
                        << mFastaMapFindIter->second.substr(refStart-1, refSeq.length()) 
                        << " -> " << refSeq 
                        << endl;
                exit(1);
            }

            // ͼ�еĽڵ���
            uint32_t nodeNum = mGraphMap.size();
            
            // ��� genotype=0 �Ľڵ�
            // ���Ⱦɫ���һ��vcfλ�ò���1����ѻ�����ǰ�ߵ����й�����һ��node
            uint16_t genotype = 0;
            if (chromosome != tmpChromosome && refStart > 0) {
                // �����Ⱦɫ���ǰ����
                string preRefSeq;
                uint32_t preRefStart;
                uint32_t preRefEnd;
                uint32_t preRefLen;
                
                // ��Ӿ�Ⱦɫ��ĺ�벿��
                if (nodeNum >= 1 && tmpRefEnd < mFastaMap[tmpChromosome].length()) {  // �ж����һ�������ڲ���Ⱦɫ�����
                    preRefStart = tmpRefEnd + 1;  // ��һ����ʼλ��
                    preRefEnd = mFastaMap[tmpChromosome].length();  // ��һ����ֹλ�ã�Ⱦɫ���ĩβ
                    preRefLen = preRefEnd - preRefStart + 1;  // ��һ�����еĳ���
                    preRefSeq = mFastaMap[tmpChromosome].substr(preRefStart-1, preRefLen);  // ��һ��Ⱦɫ����󲿷ֵ�������Ϣ

                    // Variable Binding
                    nodeSrt& mGraphMapForChrForStart = mGraphMap[tmpChromosome][preRefStart];

                    mGraphMapForChrForStart.seqMap[genotype] = preRefSeq;  // ��ӵ�ͼ���У�0
                    mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                }

                // �����Ⱦɫ���ǰ����
                preRefStart = 1;  // ��Ⱦɫ�����ʼ��1
                preRefEnd = refStart - 1;  // ��һ�������ǰһ��λ��
                preRefLen = preRefEnd - preRefStart + 1;  // ����
                preRefSeq = mFastaMapFindIter->second.substr(0, preRefLen);  // ����
                
                // Variable Binding
                nodeSrt& mGraphMapForChrForStart = mGraphMap[chromosome][preRefStart];

                mGraphMapForChrForStart.seqMap[genotype] = preRefSeq;  // ��ӵ�ͼ���У�0
                mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
            } else {  // �����vcf�м�����й���node
                string preRefSeq;
                uint32_t preRefStart;
                uint32_t preRefEnd;
                uint32_t preRefLen;

                preRefStart = tmpRefEnd + 1;
                preRefEnd = refStart - 1;
                preRefLen = preRefEnd - preRefStart + 1;

                // �ж��ǲ������ڣ����ڵĻ�preRefLenС�ڵ���0������
                if (preRefLen > 0 && preRefStart < preRefEnd) {
                    preRefSeq = mFastaMapFindIter->second.substr(preRefStart-1, preRefLen);

                    // Variable Binding
                    nodeSrt& mGraphMapForChrForStart = mGraphMap[chromosome][preRefStart];

                    mGraphMapForChrForStart.seqMap[genotype] = preRefSeq;  // ��ӵ�ͼ���У�0
                    mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                }
            }
            
            // ���vcf�ڵ�
            // Variable Binding
            nodeSrt& mGraphMapForChrForStart = mGraphMap[chromosome][refStart];

            // �����ref�Ľڵ�
            mGraphMapForChrForStart.seqMap[genotype] = refSeq;  // ��ӵ�ͼ���У�0
            mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information

            for (auto it : qrySeqVec) {  // ��qrySeqVec�б����ѭ��
                if (genotype < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                    genotype++;  // Haplotype index +1
                } else {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: Haplotype number >= " << UINT16_MAX << endl;
                    exit(1);
                }

                mGraphMapForChrForStart.seqMap[genotype] = it;  // ��ӵ�ͼ���У�genotype

                // Record the amount of variation
                uint32_t qryLen = it.size();
                int32_t svLen = qryLen - refLen;
                double lengthRatio = qryLen / float(refLen);

                if (svLen == 0 && refLen == 1 && qryLen == 1) {
                    snpNum++;
                } else if (svLen <= 49 && svLen >= -49 && refLen <= 49 && qryLen <= 49) {
                    indelNum++;
                } else if (svLen >= -2 && svLen <= 2 && refLen > 49 && qryLen > 49) {
                    invNum++;
                } else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && refLen > 49 && qryLen > 49) {
                    dupNum++;
                } else if (svLen < 0) {
                    delNum++;
                } else if (svLen > 0) {
                    insNum++;
                } else {
                    otherNum++;
                }
            }

            uint16_t hapIdx = 1;  // ��������Ϣ����
            // ���ÿ�������͵ķ�����Ϣ
            for (size_t i = 9; i < lineVec.size(); i++) {
                vector<string> gtVec = construct_index::gt_split(split(lineVec[i], ":")[gtIndex]);  // ��ϵ��Ӧ�ķ����б�

                // ���gt_split���ؿ��б������λ��Ϊ '.'��Ϊδ�������ϵ
                if (gtVec.empty()) {
                    for (size_t j = 0; j < ploidy_; j++) {
                        if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                            mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                            hapIdx++;  // Haplotype index +1
                        } else {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                << "Error: Haplotype number >= " << UINT16_MAX << endl;
                            exit(1);
                        }
                    }
                } else {  // ������λ��
                    // ��鱶���Ƿ�Ͳ���һ��
                    if (gtVec.size() > ploidy_) {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                                << "Error: " <<  "(" << lineVec[i] << ").gtVec.size()" << " != " << ploidy_ << "." << endl;
                        exit(1);
                    }

                    // �������ٱ���
                    for (size_t j = 0; j < gtVec.size(); j++) {
                        if (gtVec[j] == ".") {
                            if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                                mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                                hapIdx++;  // Haplotype index +1
                            } else {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: Haplotype number >= " << UINT16_MAX << endl;
                                exit(1);
                            }
                        } else {
                            // uint16_t gtTmp = stoi(gtVec[j]);
                            mGraphMapForChrForStart.hapGtVec.push_back(stoi(gtVec[j]));  // haplotype information

                            if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                                hapIdx++;  // Haplotype index +1
                            } else {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: Haplotype number >= " << UINT16_MAX << endl;
                                exit(1);
                            }
                        }
                    }

                    // ����б���С��ploidy���ٵĲ�����0����
                    if (gtVec.size() < ploidy_) {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                                << "Error: haplotype number: " << gtVec.size() << " < " << ploidy_ << ", instead by '.'.\n";

                        for (size_t j = 0; j < ploidy_ - gtVec.size(); j++) {
                            // mGraphMap[chromosome][refStart].undefinedHapMap[hapIdx] = 1;  // ��¼Ϊδ�����λ��
                            mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                        
                            if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                                hapIdx++;  // Haplotype index +1
                            } else {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                        << "Error: Haplotype number >= " << UINT16_MAX << endl;
                                exit(1);
                            }
                        }
                    }
                }
            }

            tmpRefEnd = refEnd;  // ��������
            tmpChromosome = chromosome;  // ��������
        }
    }

    // ���һ��Ⱦɫ��ĺ�벿��
    if (tmpRefEnd < mFastaMap[tmpChromosome].length()) {  // �ж����һ�������ڲ���Ⱦɫ�����
        uint32_t preRefStart = tmpRefEnd + 1;
        uint32_t preRefEnd = mFastaMap[tmpChromosome].length();
        uint32_t preRefLen = preRefEnd - preRefStart + 1;

        string preRefSeq = mFastaMap[tmpChromosome].substr(preRefStart-1, preRefLen);

        // Variable Binding
        nodeSrt& mGraphMapForChrForStart = mGraphMap[tmpChromosome][preRefStart];

        mGraphMapForChrForStart.seqMap[0] = preRefSeq;  // ��ӵ�ͼ���У�0
        mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Parsed " << snpNum + indelNum + insNum + delNum + invNum + dupNum + otherNum << " alternative alleles ...\n\n";  // print log
    cerr << "           - " << "SNP: " << snpNum << endl;
    cerr << "           - " << "InDels: " << indelNum << endl;
    cerr << "           - " << "Insertion: " << insNum << endl;
    cerr << "           - " << "Deletion: " << delNum << endl;
    cerr << "           - " << "Inversion: " << invNum << endl;
    cerr << "           - " << "Duplication: " << dupNum << endl;
    cerr << "           - " << "Other: " << otherNum << endl << endl << endl;

    malloc_trim(0);	// 0 is for heap memory
}

/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0
 * @brief building the k-mer index of graph
 * 
 * @return void
**/
void ConstructIndex::index()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building the graph index ...\n";  // print log

    // ������̵߳Ľ��
    vector<future<tuple<map<uint32_t, nodeSrt>::iterator, vector<kmerHashHap> > > > futureVec;  // nodeIter, kmerHashHap.kmerHash  kmerHashHap.BitVec

    // ���̳�
    ThreadPool pool(threads_);

    // ��ʼ���̳߳�
    pool.init();

    for(auto& [chromosome, startNodeMap] : mGraphMap) {  // map<chr, map<nodeStart, nodeSrt> >
        for(auto iter = startNodeMap.begin(); iter != startNodeMap.end(); iter++) {  // map<nodeStart, nodeSrt>
            if (iter->second.hapGtVec.size() == 1) {continue;}  // Skip the node if there is only one GT (0)

            // ���߳�
            futureVec.push_back(
                pool.submit(
                    construct_index::index_run,
                    chromosome, 
                    iter, 
                    ref(startNodeMap), 
                    ref(kmerLen_), 
                    mbf, 
                    ref(ploidy_), 
                    ref(debug_)
                )
            );
        }
    }

    // �����߳̽�����浽ͼ��// �����߳̽�����浽ͼ��
    for (auto& futureResult : futureVec) {  // vector<future<{nodeIter, vector<kmerHashHap>}> >
        auto [nodeIter, tmpkmerHashHapVec] = move(futureResult.get());  // {nodeIter, kmerHashHapVec}
        nodeIter->second.kmerHashHapVec = move(tmpkmerHashHapVec);

        // ��kmerHash�ӵ�GraphKmerCovMap��GraphKmerFreMap��
        for (const auto& kmerHashHapTmp : nodeIter->second.kmerHashHapVec)  // vector<kmerHashHap>
        {
            // record frequency of the k-mers informations
            auto& value = mGraphKmerCovFreMap[kmerHashHapTmp.kmerHash];
            if (value.f < UINT8_MAX)  // ��ֹ����Խ��
            {
                value.f++;
            }
        }
    }

    // free memory
    vector<future<tuple<map<uint32_t, nodeSrt>::iterator, vector<kmerHashHap> > > >().swap(futureVec);

    malloc_trim(0);	// 0 is for heap memory

    // �ر��̳߳�
    pool.shutdown();
}


/**
 * @author zezhen du
 * @date 2023/07/14
 * @version v1.0
 * @brief Remove duplicate K-mers from graph genome
 * 
 * @return void
**/
void ConstructIndex::kmer_deduplication()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Deduplicating the graph index ...\n";  // print log

    /* ******************************** Delete k-mers with frequencies greater than 1 ******************************** */
    for (auto it = mGraphKmerCovFreMap.begin(); it != mGraphKmerCovFreMap.end(); ) {
        if (it->second.f > 1) {
            it = mGraphKmerCovFreMap.erase(it);
        } else {
            ++it;
        }
    }

    /* ******************************** Remove k-mers with frequencies greater than 1 from the nodes ******************************** */
    // Save the results of multiple threads
    vector<future<tuple<vector<kmerHashHap>&, vector<kmerHashHap> > > > futureVec;  // vector<tuple<vector<kmerHashHap>&, vector<kmerHashHap> > >

    // ���̳�
    ThreadPool pool(threads_);

    // ��ʼ���̳߳�
    pool.init();

    for (auto& [chromosome, startNodeMap] : mGraphMap) {  // map<chr, map<nodeStart, nodeSrt> >
        for (auto& [nodeStart, node] : startNodeMap) {  // map<nodeStart, nodeSrt>
            if (node.hapGtVec.size() == 1) {continue;}  // ���ֻ��0�������ýڵ�

            // ���߳�
            futureVec.push_back(
                pool.submit(
                    construct_index::kmer_deduplication_run,
                    ref(node.kmerHashHapVec), 
                    ref(mGraphKmerCovFreMap)
                )
            );
        }
    }

    // �����߳̽�����浽ͼ��// �����߳̽�����浽ͼ��
    for (auto& futureResult : futureVec) {  // vector<tuple<vector<kmerHashHap>&, vector<kmerHashHap> > >
        auto [kmerHashHapVec, kmerHashHapVecTmp] = move(futureResult.get());
        kmerHashHapVec = move(kmerHashHapVecTmp);
    }

    // free memory
    vector<future<tuple<vector<kmerHashHap>&, vector<kmerHashHap> > > >().swap(futureVec);
    
    malloc_trim(0);	// 0 is for heap memory

    // �ر��̳߳�
    pool.shutdown();
}


/**
 * @brief graph index for kmer (threads)
 * 
 * @date 2023/07/14
 * 
 * @param chromosome            mGraphMap output by construct��map<chr, map<start, nodeSrt> >
 * @param nodeIter              node iterator
 * @param startNodeMap          Chromosome all nodes
 * @param kmerLen               the length of kmer
 * @param bf                    Kmer frequency in the reference genome: Counting Bloom Filter
 * @param ploidy                Ploidy of vcf file
 * @param debug                 debug code
 * 
 * @return {nodeIter, tmpKmerHapBitMap}     kmer: map<kmerHash, vector<int8_t> >
**/
tuple<map<uint32_t, nodeSrt>::iterator, vector<kmerHashHap> > construct_index::index_run(
    string chromosome, 
    map<uint32_t, nodeSrt>::iterator nodeIter, 
    const map<uint32_t, nodeSrt>& startNodeMap,
    const uint32_t& kmerLen, 
    BloomFilter* bf, 
    const uint32_t& ploidy, 
    const bool& debug
) {
    vector<kmerHashHap> tmpkmerHashHapVec;  // Whether the storage haplotype contains the corresponding kmer: map<kmerHash, vector<int8_t> >:  0000 0000, Each bits represents a haplotype, 0->False 1->True

    unordered_map<uint64_t, vector<int8_t> > KmerHapBitMap;  // kmer: map<kmerHash, vector<int8_t> >:  0000 0000, Each bits represents a haplotype, 0->False 1->True

    const auto& seqMap = nodeIter->second.seqMap;
    const auto& hapGtVec = nodeIter->second.hapGtVec;

    uint16_t haplotype = 0;  // Index of the haplotype

    for (const auto& gt : hapGtVec) {  // Iterate over the genotypes
        // calculate whether the genotype of the corresponding sample is empty or zero. If it is, skip all its haplotypes
        // if (haplotype > 0 && gt == (uint16_t)0) {
        if (haplotype > 0) {
            uint16_t groupIdx = (haplotype - 1) / ploidy;
            uint16_t hapIdxL = groupIdx * ploidy + 1;
            uint16_t hapIdxR = (groupIdx + 1) * ploidy;

            uint16_t gtSum = std::accumulate(
                hapGtVec.begin() + hapIdxL,
                hapGtVec.begin() + hapIdxR + 1,
                0
            );
            
            if (gtSum == 0) {
                continue;
            }
        }
        
        string seqTmp;

        try {
            seqTmp = seqMap.at(gt);
        } catch (const std::out_of_range& e) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: sequence information not found -> " << chromosome << " " << nodeIter->first << " GT:" << gt << endl;
            exit(1);
        }
        
        pair<string, string> upDownSeq = construct_index::find_node_up_down_seq(haplotype, kmerLen, nodeIter, startNodeMap);  // Sequence with upstream and downstream 1 k-mer sequence
        // if (debug) {
        //     cerr << "nodeStart:" << nodeIter->first() << " haplotype:" << haplotype << " gt:" << +gt << " up:" << upDownSeq.first << " this:" << seqTmp << " down:" << upDownSeq.second << endl;
        // }
        
        seqTmp = upDownSeq.first + seqTmp + upDownSeq.second;  // Add upstream and downstream sequences
        
        // Kmer indexing
        unordered_map<uint64_t, uint8_t> kmerHashMapTmp;
        kmerBit::kmer_sketch(seqTmp, kmerLen, kmerHashMapTmp, bf);

        uint64_t quotient = DIVIDE_BY_8(haplotype);  // Variable to store the quotient
        uint64_t remainder = GET_LOW_3_BITS(haplotype); // Get the remainder

        // Add the kmerHash to the node's information
        for (const auto& [kmerHash, _] : kmerHashMapTmp) {
            // Record the haplotype information of k-mer
            auto findIter = KmerHapBitMap.emplace(kmerHash, vector<int8_t>(DIVIDE_BY_8(hapGtVec.size()) + 1, 0)).first;

            int8_t& hapBitTmp = findIter->second[quotient];
            construct_index::set_bit_to_one(hapBitTmp, remainder);  // Set the corresponding haplotype to 1

            // if (debug) {
            //     cerr << construct_index::print_bits(hapBitTmp) << " bit:" << +hapBitTmp << endl;
            // }
        }

        // Haplotype index increment
        ++haplotype;
    }

    // Convert map to vector
    construct_index::convert_map2vec(KmerHapBitMap, tmpkmerHashHapVec);

    return {nodeIter, tmpkmerHashHapVec};
}


/**
 * @author zezhen du
 * @date 2023/07/13
 * @version v1.0.1
 * @brief find the sequence information corresponding to the haplotype upstream and downstream of the node
 * 
 * @param haplotype       ����������
 * @param seqLen          kmer�ĳ���
 * @param nodeIter        startNodeMap�нڵ����ڵĵ�����
 * @param startNodeMap    Ⱦɫ�����еĽڵ���ϢstartNodeMap��index��������
 * 
 * @return pair<upSeq, downSeq>
**/
pair<string, string> construct_index::find_node_up_down_seq(
    const uint16_t & haplotype, 
    const uint32_t & seqLen,
    const map<uint32_t, nodeSrt>::iterator & nodeIter, 
    const map<uint32_t, nodeSrt> & startNodeMap
)
{
    uint32_t nodeEnd = nodeIter->first + nodeIter->second.seqMap[0].size() - 1;  // Node start position

    // ��¼�ڵ�����
    string upSeq;
    string downSeq;

    // ��ʱ������
    auto iterTmp = nodeIter;

    // ��������һ���ڵ��ĩβ����
    while (upSeq.size() < seqLen && iterTmp != startNodeMap.begin()) {
        // ������ǰ��
        iterTmp--;

        const uint32_t& nodeStartTmp = iterTmp->first;  // �ڵ����ʼλ��
        const nodeSrt& nodeTmp = iterTmp->second;  // �ڵ���Ϣ
        uint32_t nodeEndTmp = nodeStartTmp + nodeTmp.seqMap.at(0).size() - 1;  // �ڵ����ֹλ��

        uint16_t gt = (haplotype < nodeTmp.hapGtVec.size()) ? nodeTmp.hapGtVec[haplotype] : 0;

        // �����Ͷ�Ӧ������
        string seqTmp;
        try {
            seqTmp = nodeTmp.seqMap.at(gt);
        } catch (const std::out_of_range& e) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: sequence information not found -> " << nodeStartTmp << " GT:" << gt << endl;
            exit(1);
        }

        /*Ŀǰ��Ƕ��ͼ�汾��ֻ��ref�����꣬���Խڵ��ཻʱ�����Ҫ��ref�����У����ж�����ȡ�������qry�����У�����һ���ڵ�*/
        // �жϽڵ��Ƿ��ཻ
        if (nodeEndTmp == nodeIter->first) {  // ��һ���ڵ����ֹ�͸ýڵ���ʼһ��
            if (gt == 0) {  // ��Ҫ�����ref����ʱ����������
            
                seqTmp = seqTmp.substr(0, nodeIter->first - nodeStartTmp);  // �ض�ref�����У���Ҫ���һ�����
            } else {  // ��Ҫ��qry����ʱ������û��qry�����ֱ꣬����һ��ѭ��
                continue;
            }
        } else if (nodeEndTmp > nodeIter->first) {  // ��һ���ڵ���ֹ�ڸýڵ���ʼ�ĺ�
            if (gt == 0) {  // ��Ҫ�����ref����ʱ����������
                seqTmp = seqTmp.substr(0, nodeIter->first - nodeStartTmp);  // �ض�ref������
            } else {  // ��Ҫ��qry����ʱ������û��qry�����ֱ꣬����һ��ѭ��
                continue;
            }
        }
        
        if (seqTmp.size() >= seqLen - upSeq.size()) {
            upSeq = seqTmp.substr(seqTmp.size() - (seqLen - upSeq.size()), seqLen - upSeq.size()) + upSeq;
        } else {
            upSeq = seqTmp + upSeq;
        }
    }

    // ���õ�����
    iterTmp = nodeIter;

    // ��������һ���ڵ�Ŀ�ʼ���У����������Ƽӵ��ж������
    while ((downSeq.size() < seqLen) && (++iterTmp != startNodeMap.end())) {
        const uint32_t& nodeStartTmp = iterTmp->first;  // �ڵ����ʼλ��
        const nodeSrt& nodeTmp = iterTmp->second;  // �ڵ���Ϣ
        uint32_t nodeEndTmp = nodeStartTmp + nodeTmp.seqMap.at(0).size() - 1;  // �ڵ����ֹλ��

        uint16_t gt = (haplotype < nodeTmp.hapGtVec.size()) ? nodeTmp.hapGtVec[haplotype] : 0;  // the genotype in this node

        // �����Ͷ�Ӧ������
        string seqTmp;
        try {
            seqTmp = nodeTmp.seqMap.at(gt);
        } catch (const std::out_of_range& e) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: sequence information not found -> " << nodeStartTmp << " GT:" << gt << endl;
            exit(1);
        }

        /*Ŀǰ��Ƕ��ͼ�汾��ֻ��ref�����꣬���Խڵ��ཻʱ�����Ҫ��ref�����У����ж�����ȡ�������qry�����У�����һ���ڵ�*/
        // �жϽڵ��Ƿ��ཻ
        if (nodeEndTmp <= nodeEnd) {  // �����һ���ڵ���ֹС�ڸýڵ���ֹ��ֱ����һ��ѭ��
            continue;
        } else if (nodeEnd == nodeStartTmp) {  // �ýڵ����ֹ����һ���ڵ���ʼһ��
            if (gt == 0) {  // ��Ҫ�����ref����ʱ����������
                seqTmp = seqTmp.substr(1, seqTmp.size() - 1);  // �ض�ref�����У���Ҫ���һ�����
            } else {  // ��Ҫ��qry����ʱ������û��qry�����ֱ꣬����һ��ѭ��
                continue;
            }
        } else if (nodeEnd > nodeStartTmp) {  // ��һ���ڵ���ֹ�ڸýڵ���ʼ�ĺ�
            if (gt == 0) {  // ��Ҫ�����ref����ʱ����������
                seqTmp = seqTmp.substr(nodeEnd - nodeStartTmp + 1, nodeEndTmp - nodeIter->first);  // �ض�ref������
            } else {  // ��Ҫ��qry����ʱ������û��qry�����ֱ꣬����һ��ѭ��
                continue;
            }
        }

        
        if (seqTmp.size() >= seqLen - downSeq.size()) {
            downSeq = downSeq + seqTmp.substr(0, seqLen - downSeq.size());
        } else {
            downSeq = downSeq + seqTmp;
        }
    }

    return make_pair(upSeq, downSeq);
}


/**
 * @brief Remove duplicate K-mers from graph genome (threads)
 * 
 * @date 2023/07/13
 * 
 * @param kmerHashHapVec            Node k-mers information
 * @param GraphKmerCovFreMap        Total k-mers coverage and frequency information
 * 
 * @return tuple<vector<kmerHashHap>&, vector<kmerHashHap> >: tuple<ref(kmerHashHapVec), kmerHashHapVecTmp>
**/
tuple<vector<kmerHashHap>&, vector<kmerHashHap> > construct_index::kmer_deduplication_run(
    vector<kmerHashHap>& kmerHashHapVec, 
    unordered_map<uint64_t, kmerCovFre>& GraphKmerCovFreMap
)
{
    // record the k-mer which need to be deleted
    vector<uint64_t> kmerHashDelVec;

    vector<kmerHashHap> kmerHashHapVecTmp = kmerHashHapVec;  // �ýڵ����е�kmer��Ϣ

    for (auto iter1 = kmerHashHapVecTmp.begin(); iter1 != kmerHashHapVecTmp.end();) {
        auto iter2 = GraphKmerCovFreMap.find(iter1->kmerHash);  // ��kmer��ͼ�е�Ƶ��

        // If the k-mer is not found in the GraphKmerCovFreMap, it means that the frequency of this k-mer is greater than 1, so delete.
        if (iter2 == GraphKmerCovFreMap.end()) {
            // record the k-mer which need to be deleted
            kmerHashDelVec.emplace_back(iter1->kmerHash);

            iter1 = kmerHashHapVecTmp.erase(iter1);  // ͼ��ɾ��
        } else {
            iter1++;
        }
    }

    return {kmerHashHapVec, kmerHashHapVecTmp};
}


/**
 * @author zezhen du
 * @date 2023/07/13
 * @version v1.0.1
 * @brief Convert map to vector
 * 
 * @param inputMap            unordered_map<uint64_t, vector<int8_t> >
 * @param tmpkmerHashHapVec   Converted vector: vector<kmerHashHap>
 * 
 * @return void
**/
void construct_index::convert_map2vec(const unordered_map<uint64_t, vector<int8_t> >& inputMap, vector<kmerHashHap>& tmpkmerHashHapVec) {
    tmpkmerHashHapVec.reserve(inputMap.size());
    
    for (auto& kvp : inputMap) {
        tmpkmerHashHapVec.emplace_back(move(kvp.first), move(kvp.second));
    }
}

/**
 * @author zezhen du
 * @date 2023/06/11
 * @version v1.0
 * @brief split the GT information
 * 
 * @param gtTxt     GT information
 * 
 * @return vector<string> gtVecTmp
**/
vector<string> construct_index::gt_split(const string & gtTxt)
{
    // ��ʱgt�б�
    vector<string> gtVecTmp;

    // ����� '.' ֱ�ӷ��ؿ��б������Թ�
    if (gtTxt == ".")
    {
        return gtVecTmp;
    }
    
    // ��gt���в��
    if (gtTxt.find("/") != string::npos)
    {
        gtVecTmp = split(gtTxt, "/");
    }
    else if (gtTxt.find("|") != string::npos)
    {
        gtVecTmp = split(gtTxt, "|");
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "GT is not separated by '/' or '|': " << gtTxt << endl;
        exit(1);
    }

    return gtVecTmp;
}


/**
 * @author zezhen du
 * @date 2023/06/29
 * @version v1.0.1
 * @brief Sets a specific bit of a given number to 1.
 * 
 * @param num         The value that needs to be modified.
 * @param bitIndex    The index of the bit to be set to 1.
 * 
 * @return void
**/
template <typename T>
void construct_index::set_bit_to_one(T& num, int bitIndex) {
    T mask = static_cast<T>(1) << bitIndex;
    num |= mask;
}

/**
 * @author zezhen du
 * @date 2023/06/29
 * @version v1.0.1
 * @brief Sets a specific bit of a given number to 0.
 * 
 * @param num         The value that needs to be modified.
 * @param bitIndex    The index of the bit to be set to 0.
 * 
 * @return void
**/
template <typename T>
void construct_index::set_bit_to_zero(T& num, int bitIndex) {
    T mask = ~(static_cast<T>(1) << bitIndex);
    num &= mask;
}

/**
 * @author zezhen du
 * @date 2023/06/29
 * @version v1.0.1
 * @brief Retrieves the value of a specific bit in a given number.
 * 
 * @param num         The value from which to retrieve the bit.
 * @param bitIndex    The index of the bit to be queried.
 * 
 * @return int
**/
template <typename T>
int construct_index::get_bit(T num, int bitIndex) {
    T mask = static_cast<T>(1) << bitIndex;
    return ((num & mask) ? 1 : 0);
}

/**
 * @author zezhen du
 * @date 2023/07/03
 * @version v1.0.1
 * @brief print bitmap.
 * 
 * @param num         The value from which to print.
 * 
 * @return int
**/
string construct_index::print_bits(int8_t num) {
    stringstream ss;
    for (int i = 7; i >= 0; i--) {
        ss << ((num >> i) & 1);
    }
    return ss.str();
}


// Explicitly instantiate the template
// int8_t
template void construct_index::set_bit_to_one<int8_t>(int8_t&, int);
template void construct_index::set_bit_to_zero<int8_t>(int8_t&, int);
template int construct_index::get_bit<int8_t>(int8_t, int);

// int32_t
template void construct_index::set_bit_to_one<int32_t>(int32_t&, int);
template void construct_index::set_bit_to_zero<int32_t>(int32_t&, int);
template int construct_index::get_bit<int32_t>(int32_t, int);

// int64_t
template void construct_index::set_bit_to_one<int64_t>(int64_t&, int);
template void construct_index::set_bit_to_zero<int64_t>(int64_t&, int);
template int construct_index::get_bit<int64_t>(int64_t, int);