// rm src/alignment*; g++ -c src/*.cpp -std=c++17 -lz -lpthread -O3 -march=native; g++ main.cpp *.o -o varigraph -lz -lpthread -std=c++17 -O3 -march=native
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
#define PROGRAM_DATA "2023/08/31"
// define version
#define PROGRAM_VERSION "1.0.1"
// define author
#define PROGRAM_AUTHOR "Zezhen Du"
// define E-mail
#define PROGRAM_E_MAIL "dzz0539@gmail.com or dzz0539@163.com"


void help(char** argv);

int main(int argc, char** argv)
{
    // initial time
    double realtime0 = realtime();

    /* -------------------------------------------------- input/output -------------------------------------------------- */
    // reference genome
    string refFileName;

    // sequencing data
    vector<string> fastqFileNameVec;

    // VCF file
	string vcfFileName;

    // output file
    string outputFileName = "";

    // sample name
    string sampleName = "out";

    /* -------------------------------------------------- genome type -------------------------------------------------- */
    string genomeType = "homozygous";

    /* -------------------------------------------------- ploidy -------------------------------------------------- */
    // reference -> 1
    uint32_t refPloidy = 2;

    // vcf -> 2
    uint32_t vcfPloidy = 2;

    /* -------------------------------------------------- genotyping arguments -------------------------------------------------- */
    uint32_t haploidNum = 15;  // H

    /* -------------------------------------------------- algorithm arguments -------------------------------------------------- */
	// k-mer size
	uint32_t kmerLen = 27;

    // fast mode -> 3
    bool fastMode = false;

    /* -------------------------------------------------- storing/loading index -------------------------------------------------- */
    // Counting Bloom Filter
    string outputMbfFileName = "";  // 4
    string inputMbfFileName = "";  // 5

    // Genome Graph
    string outputGraphFileName = "";  // 6
    string inputGraphFileName = "";  // 7

    // k-mer index
    string outputFastqKmerFileName = "";  // 8
    string inputFastqKmerFileName = "";  // 9
    
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
            {"reference", required_argument, 0, 'r'},
            {"fastq", required_argument, 0, 'f'},
			{"vcf", required_argument, 0, 'v'},
            {"out", required_argument, 0, 'o'},
            {"name", required_argument, 0, 'n'},

            {"genotype", required_argument, 0, 'g'},

            {"genome-ploidy", required_argument, 0, '1'},
            {"vcf-ploidy", required_argument, 0, '2'},

            {"haploid", required_argument, 0, 'H'},

            {"kmer", required_argument, 0, 'k'},
            {"fast", no_argument, 0, '3'},

            {"save-cbf", required_argument, 0, '4'},
            {"load-cbf", required_argument, 0, '5'},
            {"save-graph", required_argument, 0, '6'},
            {"load-graph", required_argument, 0, '7'},
            {"save-reads", required_argument, 0, '8'},
            {"load-reads", required_argument, 0, '9'},

            {"debug", no_argument, 0, 'D'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:f:v:o:n:g:1:2:H:k:34:5:6:7:8:9:Dt:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'r':
            refFileName = optarg;
            break;
        case 'f':
            fastqFileNameVec.push_back(optarg);
            break;
        case 'v':
            vcfFileName = optarg;
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 'n':
            sampleName = optarg;
            break;

        case 'g':
            genomeType = optarg;
            break;

        case '1':
            refPloidy = max(stoi(optarg), 2);
            break;
        case '2':
            vcfPloidy = max(stoi(optarg), 2);
            break;

        case 'H':
            haploidNum = stoull(optarg);
            break;

        case 'k':
            kmerLen = max(stoi(optarg), 5);
            break;
        case '3':
            fastMode = true;
            break;

        case '4':
            outputMbfFileName = optarg;
            break;
        case '5':
            inputMbfFileName = optarg;
            break;
        case '6':
            outputGraphFileName = optarg;
            break;
        case '7':
            inputGraphFileName = optarg;
            break;
        case '8':
            outputFastqKmerFileName = optarg;
            break;
        case '9':
            inputFastqKmerFileName = optarg;
            break;

        case 'D':
            debug = true;
            break;
        case 't':
            threads = max(stoi(optarg), 1);
            break;

        case 'h':
        case '?':
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help(argv);
        return 1;
    }

    // Determine whether the parameters are correct
    if (refFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -r. The reference genome file is empty.\n\n";
		help(argv);
        return 1;
    }

    if (fastqFileNameVec.empty() && inputFastqKmerFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: No sequencing files or sequencing file index provided.\n\n";
		help(argv);
        return 1;
    }

    if (vcfFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -v. The VCF file is empty.\n\n";
		help(argv);
        return 1;
    }

    if (genomeType != "homozygous" && genomeType != "heterozygous") {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -g. The submitted parameter must be either 'homozygous' or 'heterozygous'.\n\n";
		help(argv);
        return 1;
    }

    if (refPloidy == 0 || refPloidy > 8) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -1. The submitted parameter must be between 2 and 8 (inclusive).\n\n";
		help(argv);
        return 1;
    }

    if (vcfPloidy == 0 || vcfPloidy > 8) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -2. The submitted parameter must be between 2 and 8 (inclusive).\n\n";
		help(argv);
        return 1;
    }

    if (haploidNum == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -H. The submitted parameter must be greater than 0.\n\n";
		help(argv);
        return 1;
    }

    if (haploidNum < 10) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter warning: -H. The number of haploid for genotyping is relatively small, which may affect the accuracy of genotyping..\n\n";
    }

    if (kmerLen <= 0 || kmerLen > 28) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -k. The submitted parameter must be between 1 and 28 (inclusive).\n\n";
		help(argv);
        return 1;
    }
    

    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are using varigraph (v" << PROGRAM_VERSION << ")\n\n\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Threads: " << threads << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "k-mer size: " << kmerLen << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Reference file: " << refFileName << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Variants file: " << vcfFileName << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Status of the reference genome: " << genomeType << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Genome ploidy: " << refPloidy << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Ploidy of genotypes in VCF file: " << vcfPloidy << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Haploid number for genotyping: " << haploidNum << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug: " << debug << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Fast mode: " << fastMode << endl << endl << endl;


    // construct, index and genotype
    Varigraph VarigraphClass(
        refFileName, 
        fastqFileNameVec, 
        vcfFileName, 
        inputMbfFileName, 
        inputGraphFileName, 
        inputFastqKmerFileName, 
        outputMbfFileName, 
        outputGraphFileName, 
        outputFastqKmerFileName, 
        outputFileName, 
        fastMode, 
        kmerLen, 
        sampleName, 
        genomeType, 
        refPloidy, 
        vcfPloidy, 
        haploidNum, 
        debug, 
        threads
    );

    // build the kmer index of reference and construct graph
    VarigraphClass.ref_idx_construct();

    // build the kmer index of files
    VarigraphClass.kmer_read();

    // genotype
    VarigraphClass.genotype();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n\n\n";

    // output resource usage
    fprintf(stderr, "[varigraph::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    
    return 0;
}

// help document
void help(char** argv)
{
  cerr << "usage: " << argv[0] << " -r FILE -v FILE -f FILE [options]" << endl
       << "Genotyping and phasing based on k-mer counting" << endl
       << endl
       << "data: " << PROGRAM_DATA << endl
       << "version: " << PROGRAM_VERSION << endl
       << "author: " << PROGRAM_AUTHOR << endl
       << endl
       << "input/output:" << endl
       << "    -r, --reference    FILE     input FASTA reference file" << endl
	   << "    -f, --fastq        FILE     fastq files for index building, two are allowed, one for each mate" << endl
	   << "    -v, --vcf          FILE     VCF file for index building" << endl
       << "    -o, --out          FILE     output genotyping to FILE [stdout]" << endl
       << "    -n, --name         STRING   sample name for VCF annotation line [out]" << endl
       << endl
       << "genome type:" << endl
       << "    -g, --genotype     STRING   specify the genotype of the reference genome (homozygous/heterozygous) [homozygous]" << endl
       << endl
       << "ploidy:" << endl
       << "    --genome-ploidy    INT      genome ploidy (2-8) [2]" << endl
       << "    --vcf-ploidy       INT      ploidy of genotypes in VCF file (2-8) [2]" << endl
       << endl
       << "genotyping arguments:" << endl
       << "    -H, --haploid      INT      the haploid number for genotyping [15]" << endl
       << endl
       << "algorithm arguments:" << endl
       << "    -k, --kmer         INT      k-mer size (maximum: 28) [27]" << endl
       << "    --fast                      fast mode with slightly decreased genotyping accuracy" << endl
       << endl
       << "storing/loading index:" << endl
       << "    --save-cbf         FILE     save Counting Bloom Filter index information to a file (optional)" << endl
       << "    --load-cbf         FILE     load Counting Bloom Filter index information from a file (optional)" << endl
       << "    --save-graph       FILE     save Genome Graph index information to a file (optional)" << endl
       << "    --load-graph       FILE     load Genome Graph index information from a file (optional)" << endl
       << "    --save-reads       FILE     save fastq index information to a file (optional)" << endl
       << "    --load-reads       FILE     load fastq index information from a file (optional)" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -D, --debug                 enable debug code" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}