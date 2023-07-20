// g++ main.cpp *.o -o varigraph -lz -lpthread -std=c++17 -O3 -march=native
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
#define PROGRAM_DATA "2023/06/06"
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

    // 输入测序文件
    vector<string> fastqFileNameVec;

    // 输入测序kmer信息的文件
    string inputFastqKmerFileName = "";

    // 输出测序kmer信息的文件
    string outputFastqKmerFileName = "";

    // Counting Bloom Filter
    string inputMbfFileName = "";
    string outputMbfFileName = "";

    // 参考基因组
    string refFileName;

    // 变异文件
	string vcfFileName;

    // 输出文件
    string outputFileName = "";
	
	// kmer长度
	uint32_t kmerLen = 27;

	// 线程数
	uint32_t threads = 10;

    // 前缀
    string prefix = "out";

    // 调试代码
    bool debug = false;

    // vcf中的倍型
    uint32_t ploidy = 2;

	// 输入参数
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
            {"reference", required_argument, 0, 'r'},
            {"fastq", required_argument, 0, 'f'},
			{"vcf", required_argument, 0, 'v'},
            {"out", required_argument, 0, 'o'},

            {"save-cbf", required_argument, 0, '3'},
            {"load-cbf", required_argument, 0, '4'},
            {"save", required_argument, 0, '5'},
            {"load", required_argument, 0, '6'},

            {"kmer", required_argument, 0, 'k'},
            {"prefix", required_argument, 0, 'p'},
            {"", required_argument, 0, '1'},
            {"debug", no_argument, 0, '2'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:f:v:o:3:4:5:6:k:p:1:2t:h", long_options, &option_index);

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
        case '3':
            outputMbfFileName = optarg;
            break;
        case '4':
            inputMbfFileName = optarg;
            break;
        case '5':
            inputFastqKmerFileName = optarg;
            break;
        case '6':
            outputFastqKmerFileName = optarg;
            break;
        case 'k':
            kmerLen = max(stoi(optarg), 5);
            break;
        case 't':
            threads = max(stoi(optarg), 1);
            break;
        case 'p':
            prefix = optarg;
            break;
        case '1':
            ploidy = max(stoi(optarg), 1);
            break;
        case '2':
            debug = true;
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

    // 判断参数是否正确
	if ((fastqFileNameVec.empty() && inputFastqKmerFileName.empty()) || refFileName.empty() || vcfFileName.empty())
	{
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -r -f -v\n\n";
		help(argv);
        return 1;
	}


    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are using varigraph (v" << PROGRAM_VERSION << ")\n\n\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Threads: " << threads << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "kmer size: " << kmerLen << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Reference file: " << refFileName << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Variants file: " << vcfFileName << endl << endl << endl;


    // construct, index and genotype
    Varigraph VarigraphClass(
        refFileName, 
        fastqFileNameVec, 
        vcfFileName, 
        inputFastqKmerFileName, 
        inputMbfFileName, 
        outputFastqKmerFileName, 
        outputMbfFileName, 
        outputFileName, 
        kmerLen, 
        prefix, 
        ploidy, 
        threads, 
        debug
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

// 帮助文档
void help(char** argv)
{
  cerr << "usage: " << argv[0] << " -r FILE -v FILE -f FILE [options]" << endl
       << "Genotyping and phasing based on kmer-counting" << endl
       << endl
       << "data: " << PROGRAM_DATA << endl
       << "version: " << PROGRAM_VERSION << endl
       << "author: " << PROGRAM_AUTHOR << endl
       << endl
       << "input/output:" << endl
       << "    -r, --reference    FILE     input FASTA reference" << endl
	   << "    -f, --fastq        FILE     fastq files for index building, two are allowed, one for each mate" << endl
	   << "    -v, --vcf          FILE     VCF file for index building" << endl
       << "    -o, --out          FILE     output genotyping to FILE [stdout]" << endl
       << endl
       << "storing/loading Index:" << endl
       << "    --save-cbf         FILE     save index information of Counting Bloom Filter to file (optional)" << endl
       << "    --load-cbf         FILE     load index information of Counting Bloom Filter from file (optional)" << endl
       << "    --save-reads       FILE     save index information of fastq to file (optional)" << endl
       << "    --load-reads       FILE     load index information of fastq from file (optional)" << endl
       << "    --save-graph       FILE     save index information of graph to file (optional)" << endl
       << "    --load-graph       FILE     load index information of graph from file (optional)" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -k, --kmer         INT      k-mer size (no larger than 28) [27]" << endl
       << "    -p, --prefix       STRING   for output [out]" << endl
       << "    -1                 INT      ploidy of the vcf file [2]" << endl
       << "    --debug                     debug code" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}