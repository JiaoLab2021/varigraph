// g++ main.cpp src/*.cpp -o varigraph -lz -lpthread -lstdc++fs -std=c++17 -O3 -march=native
// cmake -DCMAKE_CXX_FLAGS="-march=native" .
// cmake -DCMAKE_INSTALL_PREFIX=$path .
#include <iostream>
#include <vector>
#include "zlib.h"
#include <getopt.h>
#include <malloc.h>
#include <thread>

#include "include/varigraph.hpp"
#include "include/sys.hpp"

using namespace std;


// define data
#define PROGRAM_DATA "2024/06/19"
// define version
#define PROGRAM_VERSION "1.0.5"
// define author
#define PROGRAM_AUTHOR "Zezhen Du"
// define E-mail
#define PROGRAM_E_MAIL "dzz0539@gmail.com or dzz0539@163.com"


int main_construct(int argc, char** argv);
int main_genotype(int argc, char** argv);

void help(char** argv);
void help_construct(char** argv);
void help_genotype(char** argv);

int main(int argc, char** argv) {
    // Print help document
	if (argc == 1) {
        help(argv);
        return 1;
    }

    // Select the subfunction
	string subcommand = argv[1];

	if (subcommand == "-h" || subcommand == "--help") {
		help(argv);
        return 1;
	} else if (subcommand == "construct") {
		main_construct(argc, argv);
	} else if (subcommand == "genotype") {
        main_genotype(argc, argv);
    } else {
        cerr << "Error: ["<< argv[0] << "] command " << subcommand << " not found" << endl;
		help(argv);
        return 1;
    }

    return 0;
}

// help document
void help(char** argv) {
    cerr << "Usage: " << argv[0] << " <command> [options]" << endl
         << "Perform genotyping and phasing based on k-mer counting." << endl
         << endl
         << "Data: " << PROGRAM_DATA << endl
         << "Version: " << PROGRAM_VERSION << endl
         << "Author: " << PROGRAM_AUTHOR << endl
         << endl
         << "Subcommands:" << endl
         << "  construct        Construct a genome graph from the reference genome and variants." << endl
         << "  genotype         Perform genotyping and phasing based on k-mer counting." << endl
         << endl
         << "  -h, --help       Display this help document." << endl;
}

int main_construct(int argc, char** argv) {
    // initial time
    double realtime0 = realtime();

    /* -------------------------------------------------- input/output -------------------------------------------------- */
    // reference genome
    string refFileName;  // r

    // VCF file
	string vcfFileName;  // v

    // Genome Graph
    string outputGraphFileName = "graph.bin";  // 1

    /* -------------------------------------------------- ploidy -------------------------------------------------- */
    uint32_t vcfPloidy = 2;  // 2

    /* -------------------------------------------------- algorithm arguments -------------------------------------------------- */
	// k-mer size
	uint32_t kmerLen = 27;  // k

    // fast mode
    bool fastMode = false;  // 3

    bool useUniqueKmers = false; // 4: use only unique k-mers for indexing
    
    /* -------------------------------------------------- optional arguments -------------------------------------------------- */
    // Debug code
    bool debug = false;  // D

	// thread
	uint32_t threads = 10;  // t

	// Input parameter
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
            {"reference", required_argument, 0, 'r'},
			{"vcf", required_argument, 0, 'v'},
            {"save-graph", required_argument, 0, 1},

            {"vcf-ploidy", required_argument, 0, 2},

            {"kmer", required_argument, 0, 'k'},
            {"fast", no_argument, 0, 3},
            {"use-unique-kmers", no_argument, 0, 4},

            {"debug", no_argument, 0, 'D'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:v:1:2:k:34Dt:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'r':
            refFileName = optarg;
            break;
        case 'v':
            vcfFileName = optarg;
            break;
        case 1:
            outputGraphFileName = optarg;
            break;

        case 2:
            vcfPloidy = max(stoi(optarg), 2);
            break;

        case 'k':
            kmerLen = max(stoi(optarg), 5);
            break;
        case 3:
            fastMode = true;
            break;
        case 4:
            useUniqueKmers = true;
            break;

        case 'D':
            debug = true;
            break;
        case 't':
            threads = max(stoi(optarg), 1);
            break;

        case 'h':
        case '?':
            help_construct(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc < 3 || argv[2] == "-h" || argv[2] == "--help") {
        help_construct(argv);
        return 1;
    }

    // Validate the parameters
    if (refFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -r. The reference genome file cannot be empty.\n\n";
        help_construct(argv);
        return 1;
    }

    if (vcfFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -v. The VCF file cannot be empty.\n\n";
        help_construct(argv);
        return 1;
    }

    if (outputGraphFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --save-graph. The Genome Graph file cannot be empty.\n\n";
        help_construct(argv);
        return 1;
    }

    if (vcfPloidy <= 0 || vcfPloidy > 8) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --vcf-ploidy. The provided value must be between 2 and 8 (inclusive).\n\n";
        help_construct(argv);
        return 1;
    }

    if (kmerLen <= 0 || kmerLen > 28) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -k. The provided value must be between 1 and 28 (inclusive).\n\n";
        help_construct(argv);
        return 1;
    }
    
    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are now running varigraph (v" << PROGRAM_VERSION << ").\n\n\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Execution started ..." << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of threads: " << threads << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "k-mer size: " << kmerLen << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Reference file path: " << refFileName << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Variants file path: " << vcfFileName << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Ploidy of genotypes in the VCF file: " << vcfPloidy << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug mode: " << (debug ? "Enabled" : "Disabled") << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Fast mode: " << (fastMode ? "Enabled" : "Disabled") << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Use only unique k-mers for indexing: " << (useUniqueKmers ? "Enabled" : "Disabled") << endl << endl << endl;

    // construct
    string samplesConfigFileName = "";
    string inputGraphFileName = "";
    string sampleType = "het";
    uint32_t samplePloidy = 2;
    uint32_t haploidNum = 15;
    uint32_t chrLenThread = 1 * 1000 * 1000;
    string transitionProType = "fre";
    bool svGenotypeBool = false;
    float minSupportingReads = 0.0;
    bool useDepth = false;  // 6: use sequencing depth as the depth for homozygous k-mers

    Varigraph VarigraphClass(
        refFileName, 
        vcfFileName, 
        samplesConfigFileName, 
        inputGraphFileName, 
        outputGraphFileName, 
        fastMode, 
        kmerLen, 
        sampleType, 
        samplePloidy, 
        vcfPloidy, 
        haploidNum, 
        chrLenThread, 
        transitionProType, 
        svGenotypeBool, 
        debug, 
        threads, 
        minSupportingReads, 
        useUniqueKmers, 
        useDepth
    );

    // build the kmer index of reference and construct graph
    VarigraphClass.construct();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n\n\n";

    // output resource usage
    fprintf(stderr, "[varigraph::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return 0;
}

// help document
void help_construct(char** argv) {
    cerr << "Usage: " << argv[0] << " " << argv[1] << " -r FILE -v FILE [options]" << endl
         << "Construct a genome graph from the reference genome and variants." << endl
         << endl
         << "Input/Output Options:" << endl
         << "    -r, --reference    FILE     input FASTA reference file" << endl
         << "    -v, --vcf          FILE     VCF file for index building" << endl
         << "    --save-graph       FILE     save Genome Graph index to file [graph.bin]" << endl
         << endl
         << "ploidy:" << endl
         << "    --vcf-ploidy       INT      ploidy of genotypes in VCF file (2-8) [2]" << endl
         << endl
         << "algorithm arguments:" << endl
         << "    -k, --kmer         INT      k-mer size (maximum: 28) [27]" << endl
         << "    --fast                      enable 'fast mode' for increased speed at the cost of slightly reduced genotyping accuracy" << endl
         << "    --use-unique-kmers          use only unique k-mers for indexing" << endl
         << endl
         << "optional arguments:" << endl
         << "    -D, --debug                 enable debug code" << endl
         << "    -t, --threads      INT      number of compute threads to use [10]" << endl
         << endl
         << "    -h, --help                  display this help document" << endl;
}


int main_genotype(int argc, char** argv) {
    // initial time
    double realtime0 = realtime();

    /* -------------------------------------------------- input/output -------------------------------------------------- */
    // Genome Graph
    string inputGraphFileName = "graph.bin";  // 1

    // Sample configuration file
    string samplesConfigFileName = "";  // s

    /* -------------------------------------------------- sample type -------------------------------------------------- */
    string sampleType = "het";  // g

    /* -------------------------------------------------- ploidy -------------------------------------------------- */
    uint32_t samplePloidy = 2;  // 2

    /* -------------------------------------------------- genotyping arguments -------------------------------------------------- */
    uint32_t haploidNum = 15;  // n
    uint32_t chrLenThread = 1 * 1000 * 1000;  // 3
    string transitionProType = "fre";  // m: Transition probability type
    bool svGenotypeBool = false;  // 4: structural variation genotyping only
    float minSupportingReads = 0.0;  // 5: min-support
    bool useDepth = false;  // 6: use sequencing depth as the depth for homozygous k-mers

    /* -------------------------------------------------- optional arguments -------------------------------------------------- */
    // Debug code
    bool debug = false;

	// thread
	uint32_t threads = 10;

	// Input parameter
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
            {"load-graph", required_argument, 0, 1},
            {"sample", required_argument, 0, 's'},

            {"genotype", required_argument, 0, 'g'},

            {"sample-ploidy", required_argument, 0, 2},

            {"number", required_argument, 0, 'n'},
            {"granularity", required_argument, 0, 3},
            {"mode", required_argument, 0, 'm'},
            {"sv", no_argument, 0, 4},
            {"min-support", required_argument, 0, 5},
            {"use-depth", no_argument, 0, 6},
            
            {"debug", no_argument, 0, 'D'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "1:s:g:2:n:3:m:45:6Dt:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 1:
            inputGraphFileName = optarg;
            break;
        case 's':
            samplesConfigFileName = optarg;
            break;

        case 'g':
            sampleType = optarg;
            break;

        case 2:
            samplePloidy = max(stoi(optarg), 2);
            break;

        case 'n':
            haploidNum = stoull(optarg);
            break;
        case 3:
            chrLenThread = stof(optarg) * 1e6;
            break;
        case 'm':
            transitionProType = optarg;
            break;
        case 4:
            svGenotypeBool = true;
            break;
        case 5:
            minSupportingReads = stof(optarg);
            break;
        case 6:
            useDepth = true;
            break;

        case 'D':
            debug = true;
            break;
        case 't':
            threads = max(stoi(optarg), 1);
            break;

        case 'h':
        case '?':
            help_genotype(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc < 3 || argv[2] == "-h" || argv[2] == "--help") {
        help_genotype(argv);
        return 1;
    }

    // Validate the parameters
    if (inputGraphFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --load-graph. The genome graph file cannot be empty.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (samplesConfigFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -s. The sample configuration file cannot be empty.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (sampleType != "hom" && sampleType != "het") {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -g. The provided value must be either 'hom' or 'het'.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (samplePloidy == 0 || samplePloidy > 8) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --sample-ploidy. The provided value must be between 2 and 8 (inclusive).\n\n";
        help_genotype(argv);
        return 1;
    }

    if (haploidNum == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -n. The provided value must be greater than 0.\n\n";
        help_genotype(argv);
        return 1;
    } else if (haploidNum < 10) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: -n. The number of haploids for genotyping is relatively low, which may affect the accuracy of genotyping.\n\n";
    }

    if (chrLenThread < 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --granularity. The chromosome granularity must be greater than 1 (" << chrLenThread << ").\n\n";
        help_genotype(argv);
        return 1;
    } else if (chrLenThread < 1000) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: --granularity. The chromosome granularity is less than 1000 (" << chrLenThread << ").\n\n";
    }

    if (transitionProType != "fre" && transitionProType != "rec") {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -m. The transition probability type must be either 'fre' or 'rec'.\n\n";
        help_genotype(argv);
        return 1;
    }

    // Print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are now running varigraph (v" << PROGRAM_VERSION << ").\n\n\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Execution started ..." << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of threads: " << threads << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Genome graph file: " << inputGraphFileName << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample configuration file: " << samplesConfigFileName << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample genome status: " << sampleType << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample ploidy: " << samplePloidy << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of haploids for genotyping: " << haploidNum << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Chromosome granularity: " << chrLenThread << " bp" << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Transition probability type: " << transitionProType << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Structural variation genotyping only: " << (svGenotypeBool ? "Enabled" : "Disabled") << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Minimum site support: " << minSupportingReads << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Use sequencing depth as the depth for homozygous k-mers: " << (useDepth ? "Enabled" : "Disabled") << endl;

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug mode: " << (debug ? "Enabled" : "Disabled") << endl << endl << endl;

    // construct, index and genotype
    string refFileName = "";
    string vcfFileName = "";
    string outputGraphFileName = "";
    bool fastMode = false;
    uint32_t kmerLen = 27;
    uint32_t vcfPloidy = 2;
    bool useUniqueKmers = false;

    Varigraph VarigraphClass(
        refFileName, 
        vcfFileName, 
        samplesConfigFileName, 
        inputGraphFileName, 
        outputGraphFileName, 
        fastMode, 
        kmerLen, 
        sampleType, 
        samplePloidy, 
        vcfPloidy, 
        haploidNum, 
        chrLenThread, 
        transitionProType, 
        svGenotypeBool, 
        debug, 
        threads, 
        minSupportingReads, 
        useUniqueKmers, 
        useDepth
    );

    // parse the sample configuration file
    VarigraphClass.parse_sample_config();

    // load the genome graph from file
    VarigraphClass.load();

    // fastq and genotype
    VarigraphClass.fastq_genotype();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n\n\n";

    // output resource usage
    fprintf(stderr, "[varigraph::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    
    return 0;
}

// help document
void help_genotype(char** argv) {
    cerr << "Usage: " << argv[0] << " " << argv[1] << " --load-graph FILE -s FILE [options]" << endl
         << "Perform genotyping and phasing based on k-mer counting." << endl
         << endl
         << "Input/Output Options:" << endl
         << "    --load-graph       FILE     load Genome Graph index from file generated by '" << argv[0] << " construct' [graph.bin]" << endl
         << "    -s, --samples      FILE     Specify the samples configuration file, format: sample read1.fq.gz read2.fq.gz" << endl
         << endl
         << "Examples:"
         << endl
         << "  # To specify the samples configuration, use the following format in your samples file:" << endl
         << "      sample1 sample1.1.fq.gz sample1.2.fq.gz" << endl
         << "      sample2 sample2.1.fq.gz sample2.2.fq.gz" << endl
         << "      ..." << endl
         << endl
         << "genome type:" << endl
         << "    -g, --genotype     STRING   specify the sample genotype as either 'hom' for homozygous or 'het' for heterozygous [het]" << endl
         << endl
         << "ploidy:" << endl
         << "    --sample-ploidy    INT      sample ploidy (2-8) [2]" << endl
         << endl
         << "genotyping arguments:" << endl
         << "    -n, --number       INT      the haploid number for genotyping [15]" << endl
         << "    --granularity      FLOAT    control the chromosome length processed by each thread (unit: Mb) [1]" << endl
         << "    -m, --mode         STRING   using haplotype frequency (fre) or recombination rate (rec) as transition probability [fre]" << endl
         << "    --sv                        structural variation genotyping only" << endl
         << "    --min-support      FLOAT    minimum site support (N) for genotype [0]" << endl
         << "    --use-depth                 use sequencing depth as the depth for homozygous k-mers" << endl
         << endl
         << "optional arguments:" << endl
         << "    -D, --debug                 enable debug code" << endl
         << "    -t, --threads      INT      number of compute threads to use [10]" << endl
         << endl
         << "    -h, --help                  display this help document" << endl;
}