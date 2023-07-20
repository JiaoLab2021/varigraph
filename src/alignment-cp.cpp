// g++ alignment.cpp -o alignment -lpthread -lz
#include <iostream>
#include <vector>
#include "zlib.h"
#include "map"
#include <unordered_map>
#include <getopt.h>
#include "src/get_time.hpp"
#include "src/split.hpp"
#include "src/alignment-minimizer.hpp"
#include <malloc.h>

using namespace std;
void help_alignment(char** argv);

int main(int argc, char** argv)
{
    // 输入文件
    string refgenome = "";
    vector<string> readFile;

    // 输出文件
    string outputFileName;
	
	// kmer和minimizer的长度
	int32_t kmerLen = 31;
	int minimizerK = 15;
	int minimizerW = 10;

	// 线程数
	uint32_t threads = 10;

	// 输入参数
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
			{"refgenome", required_argument, 0, 'r'},
            {"fastq", required_argument, 0, 'f'},

            {"kmerLen", required_argument, 0, 'l'},
            {"minimizerK", required_argument, 0, 'k'},
            {"minimizerW", required_argument, 0, 'w'},
            {"threads", required_argument, 0, 't'},
            {"output", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:f:l:k:w:t:o:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'r':
            refgenome = optarg;
            break;
        case 'f':
            readFile.push_back(optarg);
            break;
        case 'l':
            kmerLen = stoi(optarg);
            break;
        case 'k':
            minimizerK = stoi(optarg);
            break;
        case 'w':
            minimizerW = stoi(optarg);
            break;
        case 't':
            threads = stoi(optarg);
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 'h':
        case '?':
            help_alignment(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_alignment(argv);
        return 1;
    }

    // 检查参数是否正确
    if (refgenome.size() == 0 || readFile.size() == 0)
    {
        help_alignment(argv);
        return 1;
    }

    // 记录minimizer索引
    unordered_map<string, unordered_map<long long int, vector<tuple<long long int, int>>>> minimizerRefMap;
    // map<染色体, 序列>
    map<string, string> refFastaMap = alignment::read_fasta(refgenome, 
                                                            kmerLen, 
                                                            minimizerK, 
                                                            minimizerW, 
                                                            threads, 
                                                            minimizerRefMap);

    // 测序文件比对
    // SR
    if (readFile.size() == 1)
    {
        alignment::s_read_fastq(readFile[0], 
                            kmerLen, 
                            minimizerK, 
                            minimizerW, 
                            threads,
                            minimizerRefMap, 
                            refFastaMap);
    }
    // PR
    else if (readFile.size() == 2)
    {
        alignment::p_read_fastq(readFile[0], 
                                readFile[1], 
                                kmerLen, 
                                minimizerK, 
                                minimizerW, 
                                threads,
                                minimizerRefMap, 
                                refFastaMap);
    }
    else
    {
        cerr << "[" << getTime() << "] " 
                << "The number of input files should be one or two: " << readFile.size()
                << endl;
            exit(1);
    }

    // 释放内存
    minimizerRefMap.clear();
    unordered_map<string, unordered_map<long long int, vector<tuple<long long int, int>>>>().swap(minimizerRefMap);

    refFastaMap.clear();
    map<string, string>().swap(refFastaMap);
    sleep(10);

    malloc_trim(0);	// 0 is for heap memory

    return 0;
}

// 帮助文档
void help_alignment(char** argv)
{
  cerr << "usage: " << argv[0] << " alignment [options] -r FILE -f FILE [-F FILE]" << endl
       << "alignment from refgenome and qrygenome file:" << endl
       << endl
       << "required arguments:" << endl
	   << "    -r, --refgenome    FILE     input FASTA refgenome" << endl
       << "    -f, --fastq       FILE      input fastq or fasta, possibly compressed, two are allowed, one for each mate" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -l, --kmerLen      INT      k-mer size [31]" << endl
       << "    -k, --minimizerK   INT      minimizer-mer size [15]" << endl
       << "    -w, --minimizerW   INT      minimizer window size [10]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << "    -o, --output       FILE     output filename [stdout]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}