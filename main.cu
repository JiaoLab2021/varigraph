// nvcc --extended-lambda -G -o varigraph-gpu main.cu src/*.cpp src/*.cu -lz -lstdc++fs -std=c++17 -O3 -Xcompiler "-march=native"
#include <iostream>
#include <vector>
#include "zlib.h"
#include <getopt.h>
#include <malloc.h>
#include <thread>

#include "include/varigraph.cuh"
#include "include/sys.hpp"

using namespace std;


// define data
#define PROGRAM_DATA "2024/09/25"
// define version
#define PROGRAM_VERSION "1.0.8"
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

    // Varigraph configuration
    VarigraphKernelConfig config;

	// Input parameter
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
            /* -------------------------------------------------- input/output -------------------------------------------------- */
            {"reference", required_argument, 0, 'r'},    // refFileName
			{"vcf", required_argument, 0, 'v'},          // vcfFileName
            {"save-graph", required_argument, 0, 1},     // outputGraphFileName
            /* -------------------------------------------------- ploidy -------------------------------------------------- */
            {"vcf-ploidy", required_argument, 0, 2},     // vcfPloidy
            /* -------------------------------------------------- algorithm arguments -------------------------------------------------- */
            {"kmer", required_argument, 0, 'k'},         // kmerLen
            {"fast", no_argument, 0, 3},                 // fastMode
            {"use-unique-kmers", no_argument, 0, 4},     // useUniqueKmers
            /* -------------------------------------------------- gpu arguments -------------------------------------------------- */
            {"gpu", required_argument, 0, 5},            // gpu
            {"buffer", required_argument, 0, 6},         // buffer
            /* -------------------------------------------------- optional arguments -------------------------------------------------- */
            {"debug", no_argument, 0, 'D'},              // debug
            {"threads", required_argument, 0, 't'},      // threads
            {"help", no_argument, 0, 'h'},               // help
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:v:1:2:k:345:6:Dt:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'r':
            config.refFileName = optarg;
            break;
        case 'v':
            config.vcfFileName = optarg;
            break;
        case 1:
            config.outputGraphFileName = optarg;
            break;

        case 2:
            config.vcfPloidy = max(stoi(optarg), 2);
            break;

        case 'k':
            config.kmerLen = max(stoi(optarg), 5);
            break;
        case 3:
            config.fastMode = true;
            break;
        case 4:
            config.useUniqueKmers = true;
            break;

        case 5:
            config.gpu = stoi(optarg);
            break;
        case 6:
            config.buffer = stoi(optarg);
            break;

        case 'D':
            config.debug = true;
            break;
        case 't':
            config.threads = max(stoi(optarg), 1);
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
    if (config.refFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -r. The reference genome file cannot be empty.\n\n";
        help_construct(argv);
        return 1;
    }

    if (config.vcfFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -v. The VCF file cannot be empty.\n\n";
        help_construct(argv);
        return 1;
    }

    if (config.outputGraphFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --save-graph. The Genome Graph file cannot be empty.\n\n";
        help_construct(argv);
        return 1;
    }

    if (config.vcfPloidy <= 0 || config.vcfPloidy > 8) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --vcf-ploidy. The provided value must be between 2 and 8 (inclusive).\n\n";
        help_construct(argv);
        return 1;
    }

    if (config.kmerLen <= 0 || config.kmerLen > 28) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -k. The provided value must be between 1 and 28 (inclusive).\n\n";
        help_construct(argv);
        return 1;
    }

    if (config.gpu < 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --gpu. The provided value must be greater than or equal to 0.\n\n";
        help_construct(argv);
        return 1;
    }

    if (config.buffer <= 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --buffer. The provided value must be greater than 0.\n\n";
        help_construct(argv);
        return 1;
    } else if (config.buffer > 1000) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: --buffer. The buffer size is relatively large, which may consume a large amount of GPU memory.\n\n";
    }
    
    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are now running varigraph (v" << PROGRAM_VERSION << ").\n\n\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Execution started ..." << endl;
    config.logConstructionConfigKernel();  // print configuration log

    // select the GPU device
    cudaSetDevice(config.gpu);
    // get the GPU device properties
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, config.gpu);
    cerr << endl;
    cerr << "           - " << "Selected GPU: " << deviceProp.name << endl;
    cerr << "           - " << "Compute capability: " << deviceProp.major << "." << deviceProp.minor << endl;
    cerr << "           - " << "GPU memory: " << deviceProp.totalGlobalMem / 1024 / 1024 / 1024 << " GB" << endl << endl << endl;

    // construct
    VarigraphKernel VarigraphKernelClass(config);

    // build the kmer index of reference and construct graph
    VarigraphKernelClass.construct_kernel();

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
         << "GPU arguments:" << endl
         << "    --gpu              INT      specify which GPU to use [0]" << endl
         << "    --buffer           INT      specify the size of the GPU buffer in MB (larger buffer results in faster processing but higher memory usage) [100]" << endl
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

    // Varigraph configuration
    VarigraphKernelConfig config;

	// Input parameter
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
            {"load-graph", required_argument, 0, 1},       // inputGraphFileName
            {"sample", required_argument, 0, 's'},         // samplesConfigFileName

            {"genotype", required_argument, 0, 'g'},       // sampleType

            {"sample-ploidy", required_argument, 0, 2},    // samplePloidy

            {"number", required_argument, 0, 'n'},         // haploidNum
            {"granularity", required_argument, 0, 3},      // chrLenThread
            {"mode", required_argument, 0, 'm'},           // transitionProType
            {"sv", no_argument, 0, 4},                     // svGenotypeBool
            {"min-support", required_argument, 0, 5},      // minSupportingGQ
            {"use-depth", no_argument, 0, 6},              // useDepth

            {"gpu", required_argument, 0, 7},              // gpu
            {"buffer", required_argument, 0, 8},           // buffer
            
            {"debug", no_argument, 0, 'D'},                // debug
            {"threads", required_argument, 0, 't'},        // threads
            {"help", no_argument, 0, 'h'},                 // help
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "1:s:g:2:n:3:m:45:67:8:Dt:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 1:
            config.inputGraphFileName = optarg;
            break;
        case 's':
            config.samplesConfigFileName = optarg;
            break;

        case 'g':
            config.sampleType = optarg;
            break;

        case 2:
            config.samplePloidy = max(stoi(optarg), 2);
            break;

        case 'n':
            config.haploidNum = stoull(optarg);
            break;
        case 3:
            config.chrLenThread = stof(optarg) * 1e6;
            break;
        case 'm':
            config.transitionProType = optarg;
            break;
        case 4:
            config.svGenotypeBool = true;
            break;
        case 5:
            config.minSupportingGQ = stof(optarg);
            break;
        case 6:
            config.useDepth = true;
            break;

        case 7:
            config.gpu = stoi(optarg);
            break;
        case 8:
            config.buffer = stoi(optarg);
            break;

        case 'D':
            config.debug = true;
            break;
        case 't':
            config.threads = max(stoi(optarg), 1);
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
    if (config.inputGraphFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --load-graph. The genome graph file cannot be empty.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (config.samplesConfigFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -s. The sample configuration file cannot be empty.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (config.sampleType != "hom" && config.sampleType != "het") {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -g. The provided value must be either 'hom' or 'het'.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (config.samplePloidy == 0 || config.samplePloidy > 8) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --sample-ploidy. The provided value must be between 2 and 8 (inclusive).\n\n";
        help_genotype(argv);
        return 1;
    }

    if (config.haploidNum == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -n. The provided value must be greater than 0.\n\n";
        help_genotype(argv);
        return 1;
    } else if (config.haploidNum < 10) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: -n. The number of haploids for genotyping is relatively low, which may affect the accuracy of genotyping.\n\n";
    }

    if (config.chrLenThread < 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --granularity. The chromosome granularity must be greater than 1 (" << config.chrLenThread << ").\n\n";
        help_genotype(argv);
        return 1;
    } else if (config.chrLenThread < 1000) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: --granularity. The chromosome granularity is less than 1000 bp (" << config.chrLenThread << " bp).\n\n";
    }

    if (config.transitionProType != "fre" && config.transitionProType != "rec") {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -m. The transition probability type must be either 'fre' or 'rec'.\n\n";
        help_genotype(argv);
        return 1;
    }

    if (config.buffer <= 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: --buffer. The provided value must be greater than 0.\n\n";
        help_construct(argv);
        return 1;
    } else if (config.buffer > 1000) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: --buffer. The buffer size is relatively large, which may consume a large amount of GPU memory.\n\n";
    }

    // Print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are now running varigraph (v" << PROGRAM_VERSION << ").\n\n\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Execution started ..." << endl;
    config.logGenotypeConfigKernel();  // print configuration log

    // select the GPU device
    cudaSetDevice(config.gpu);
    // get the GPU device properties
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, config.gpu);
    cerr << endl;
    cerr << "           - " << "Selected GPU: " << deviceProp.name << endl;
    cerr << "           - " << "Compute capability: " << deviceProp.major << "." << deviceProp.minor << endl;
    cerr << "           - " << "GPU memory: " << deviceProp.totalGlobalMem / 1024 / 1024 / 1024 << " GB" << endl << endl << endl;

    // construct, index and genotype
    VarigraphKernel VarigraphKernelClass(config);

    // parse the sample configuration file
    VarigraphKernelClass.parse_sample_config();

    // load the genome graph from file
    VarigraphKernelClass.load();

    // fastq and genotype
    VarigraphKernelClass.fastq_genotype_kernel();

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
         << "    -m, --mode         STRING   using haplotype frequency (fre) or recombination rate (rec) as transition probability [rec]" << endl
         << "    --sv                        structural variation genotyping only" << endl
         << "    --min-support      FLOAT    minimum site quality (GQ) value for genotype [0]" << endl
         << "    --use-depth                 use sequencing depth as the depth for homozygous k-mers" << endl
         << endl
         << "GPU arguments:" << endl
         << "    --gpu              INT      specify which GPU to use [0]" << endl
         << "    --buffer           INT      specify the size of the GPU buffer in MB (larger buffer results in faster processing but higher memory usage) [500]" << endl
         << endl
         << "optional arguments:" << endl
         << "    -D, --debug                 enable debug code" << endl
         << "    -t, --threads      INT      number of compute threads to use [10]" << endl
         << endl
         << "    -h, --help                  display this help document" << endl;
}