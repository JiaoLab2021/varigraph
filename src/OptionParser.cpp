#include "../include/OptionParser.hpp"

#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/HelpPrinter.hpp"
#include "../include/ProgramMetadata.hpp"
#include "../include/logger.hpp"

namespace {

constexpr int max_gpu_buffer_mb = std::numeric_limits<int>::max() / 1024 / 1024;

enum LongOption {
    SaveGraph = 1000,
    LoadGraph,
    VcfPloidy = 2000,
    SamplePloidy,
    Fast,
    UseUniqueKmers,
    Granularity,
    StructuralVariants,
    MinSupport,
    UseDepth,
    Gpu = 3000,
    Buffer
};

// Print the commands and program information.
void print_help(const char* executable, bool advanced = false) {
    diagnostic_stream() << "Usage: " << executable << " <command> [options]\n\n"
              << program::description << "\n\n"
              << "Version: " << program::version << '\n'
              << "Date: " << program::build_date << '\n';
    if (advanced) {
        diagnostic_stream() << "Author: " << program::author << '\n'
                            << "Contact: " << program::email << '\n';
    }
    diagnostic_stream() << '\n';

    HelpPrinter help;
    help.section("Commands (run in this order)");
    help.line("construct", "", "build a genome graph index from a reference FASTA and population VCF");
    help.line("genotype", "", "genotype samples from their reads using a saved graph index");
    help.blank();
    help.section("General options");
    help.line("-h, --help", "", "show common help");
    help.line("-H, --advanced", "", "show complete help");
    help.line("--version", "", "show the program version");
    help.print(diagnostic_stream());
}

// Print command options using the configuration defaults.
void print_command_help(const char* executable, Command command, bool gpu_enabled, bool advanced = false) {
    const AppConfig defaults;
    const auto& config = defaults.common;
    HelpPrinter help;

    if (command == Command::Construct) {
        diagnostic_stream() << "Usage: " << executable << " construct -r FILE -v FILE [options]\n\n"
                  << "Construct a genome graph from the reference genome and variants.\n\n";
        help.section("Input/output");
        help.line("-r, --reference", "FILE", "input reference FASTA");
        help.line("-v, --vcf", "FILE", "population VCF for index building");
        help.line("--save-graph", "FILE", "save the graph index [" + config.outputGraphFileName + "]");
        help.blank();
        help.section("Ploidy");
        help.line("--vcf-ploidy", "INT", "ploidy of genotypes in the VCF, 2-8 [" + std::to_string(config.vcfPloidy) + "]");
        help.blank();
        help.section("Indexing");
        help.line("-k, --kmer", "INT", "k-mer size, values below 5 use 5; maximum 28 [" + std::to_string(config.kmerLen) + "]");
        if (advanced) {
            help.line("--fast", "", "build faster with a possible reduction in genotyping accuracy");
            help.line("--use-unique-kmers", "", "use only unique k-mers for indexing");
        }
    } else {
        diagnostic_stream() << "Usage: " << executable << " genotype --load-graph FILE -s FILE [options]\n\n"
                  << "Perform genotyping and phasing based on k-mer counting.\n\n";
        help.section("Input/output");
        help.line("--load-graph", "FILE", "load a graph index made by construct [" + config.inputGraphFileName + "]");
        help.line("-s, --sample", "FILE", "sample configuration; each line: SAMPLE READS ...");
        help.blank();
        help.section("Sample");
        help.line("-g, --genotype", "STR", "sample genotype: hom (homozygous) or het (heterozygous) [" + config.sampleType + "]");
        help.line("--sample-ploidy", "INT", "sample ploidy, 2-8 [" + std::to_string(config.samplePloidy) + "]");
        help.blank();
        help.section("Genotyping");
        if (advanced) {
            std::ostringstream granularity;
            granularity << config.chrLenThread / 1e6;
            std::ostringstream min_support;
            min_support << config.minSupportingGQ;
            help.line("-n, --number", "INT", "number of haplotypes used for genotyping [" + std::to_string(config.haploidNum) + "]");
            help.line("--granularity", "FLOAT", "chromosome length processed per thread, in Mb [" + granularity.str() + "]");
            help.line("-m, --mode", "STR", "transition probability: fre (haplotype frequency) or rec (recombination rate) [" + config.transitionProType + "]");
            help.line("--min-support", "FLOAT", "minimum genotype quality (GQ) [" + min_support.str() + "]");
        }
        help.line("--sv", "", "genotype structural variants only");
        help.line("--use-depth", "", "use sequencing depth as the depth of homozygous k-mers");
    }

    if (gpu_enabled) {
        help.blank();
        help.section("GPU");
        help.line("--gpu", "INT", "GPU device ID [" + std::to_string(defaults.gpu.gpu) + "]");
        help.line("--buffer", "INT", "GPU buffer size in MB, maximum " + std::to_string(max_gpu_buffer_mb) + " [" + std::to_string(defaults.gpu.buffer) + "]");
    }

    help.blank();
    help.section("General options");
    help.line("-t, --threads", "INT", "number of compute threads [" + std::to_string(config.threads) + "]");
    if (advanced) help.line("-D, --debug", "", "enable debug output");
    help.line("-h, --help", "", "show common options");
    help.line("-H, --advanced", "", "show all options");
    help.print(diagnostic_stream());
}

// Validate parameter constraints before starting the workflow.
std::string validate_options(Command command, const AppConfig& options, bool gpu_enabled) {
    const auto& config = options.common;

    if (command == Command::Construct) {
        if (config.refFileName.empty()) return "-r. The reference genome file cannot be empty.";
        if (config.vcfFileName.empty()) return "-v. The VCF file cannot be empty.";
        if (config.outputGraphFileName.empty()) return "--save-graph. The genome graph file cannot be empty.";
        if (config.vcfPloidy == 0 || config.vcfPloidy > 8) {
            return "--vcf-ploidy. The provided value must be between 2 and 8 (inclusive).";
        }
        if (config.kmerLen == 0 || config.kmerLen > 28) {
            return "-k. The k-mer size must not exceed 28; values below 5 use 5.";
        }
        if (gpu_enabled && options.gpu.gpu < 0) {
            return "--gpu. The provided value must be greater than or equal to 0.";
        }
    } else {
        if (config.inputGraphFileName.empty()) return "--load-graph. The genome graph file cannot be empty.";
        if (config.samplesConfigFileName.empty()) return "-s. The sample configuration file cannot be empty.";
        if (config.sampleType != "hom" && config.sampleType != "het") {
            return "-g. The provided value must be either 'hom' or 'het'.";
        }
        if (config.samplePloidy == 0 || config.samplePloidy > 8) {
            return "--sample-ploidy. The provided value must be between 2 and 8 (inclusive).";
        }
        if (config.haploidNum == 0) return "-n. The provided value must be greater than 0.";
        if (config.chrLenThread < 1) {
            return "--granularity. The chromosome granularity must be at least 1 bp (" + std::to_string(config.chrLenThread) + ").";
        }
        if (config.transitionProType != "fre" && config.transitionProType != "rec") {
            return "-m. The transition probability type must be either 'fre' or 'rec'.";
        }
    }

    if (gpu_enabled) {
        if (options.gpu.buffer <= 0) return "--buffer. The provided value must be greater than 0.";
        if (options.gpu.buffer > max_gpu_buffer_mb) {
            return "--buffer. The provided value must not exceed " + std::to_string(max_gpu_buffer_mb) + " MB.";
        }
    }
    return {};
}

}

// Parse the selected command without starting computation.
ParseResult parse_command(int argc, char** argv, bool gpu_enabled) {
    ParseResult result;
    if (argc < 2) {
        print_help(argv[0]);
        result.status = 1;
        return result;
    }

    const std::string command = argv[1];
    if (command == "-h" || command == "--help" || command == "-H" || command == "--advanced") {
        print_help(argv[0], command == "-H" || command == "--advanced");
        return result;
    }
    if (command == "--version") {
        diagnostic_stream() << program::version << '\n';
        return result;
    }
    if (command == "construct") {
        result.command = Command::Construct;
    } else if (command == "genotype") {
        result.command = Command::Genotype;
    } else {
        error_stream() << "Error: [" << argv[0] << "] command " << command << " not found\n";
        print_help(argv[0]);
        result.status = 1;
        return result;
    }

    // ------ Options ------
    std::vector<option> long_options;
    if (result.command == Command::Construct) {
        long_options = {
            {"reference", required_argument, nullptr, 'r'},
            {"vcf", required_argument, nullptr, 'v'},
            {"save-graph", required_argument, nullptr, SaveGraph},
            {"vcf-ploidy", required_argument, nullptr, VcfPloidy},
            {"kmer", required_argument, nullptr, 'k'},
            {"fast", no_argument, nullptr, Fast},
            {"use-unique-kmers", no_argument, nullptr, UseUniqueKmers}
        };
    } else {
        long_options = {
            {"load-graph", required_argument, nullptr, LoadGraph},
            {"sample", required_argument, nullptr, 's'},
            {"genotype", required_argument, nullptr, 'g'},
            {"sample-ploidy", required_argument, nullptr, SamplePloidy},
            {"number", required_argument, nullptr, 'n'},
            {"granularity", required_argument, nullptr, Granularity},
            {"mode", required_argument, nullptr, 'm'},
            {"min-support", required_argument, nullptr, MinSupport},
            {"sv", no_argument, nullptr, StructuralVariants},
            {"use-depth", no_argument, nullptr, UseDepth}
        };
    }
    if (gpu_enabled) {
        long_options.push_back({"gpu", required_argument, nullptr, Gpu});
        long_options.push_back({"buffer", required_argument, nullptr, Buffer});
    }
    long_options.push_back({"threads", required_argument, nullptr, 't'});
    long_options.push_back({"debug", no_argument, nullptr, 'D'});
    long_options.push_back({"help", no_argument, nullptr, 'h'});
    long_options.push_back({"advanced", no_argument, nullptr, 'H'});
    long_options.push_back({nullptr, 0, nullptr, 0});

    const char* short_options = result.command == Command::Construct ? "r:v:k:t:DhH" : "s:g:n:m:t:DhH";
    auto& config = result.config.common;
    optind = 0;
    int option_code;
    while ((option_code = getopt_long(argc, argv, short_options, long_options.data(), nullptr)) != -1) {
        switch (option_code) {
        case 'r': config.refFileName = optarg; break;
        case 'v': config.vcfFileName = optarg; break;
        case SaveGraph: config.outputGraphFileName = optarg; break;
        case LoadGraph: config.inputGraphFileName = optarg; break;
        case 's': config.samplesConfigFileName = optarg; break;
        case 'g': config.sampleType = optarg; break;
        case VcfPloidy: config.vcfPloidy = std::max(std::stoi(optarg), 2); break;
        case SamplePloidy: config.samplePloidy = std::max(std::stoi(optarg), 2); break;
        case 'k': config.kmerLen = std::max(std::stoi(optarg), 5); break;
        case Fast: config.fastMode = true; break;
        case UseUniqueKmers: config.useUniqueKmers = true; break;
        case 'n': config.haploidNum = std::stoull(optarg); break;
        case Granularity: {
            const double bases = std::stof(optarg) * 1e6;
            if (!std::isfinite(bases) || bases < 0.0 ||
                bases >= static_cast<double>(std::numeric_limits<std::uint32_t>::max()) + 1.0) {
                throw std::out_of_range("--granularity must be finite and representable as an unsigned 32-bit base count.");
            }
            config.chrLenThread = static_cast<std::uint32_t>(bases);
            break;
        }
        case 'm': config.transitionProType = optarg; break;
        case MinSupport: config.minSupportingGQ = std::stof(optarg); break;
        case StructuralVariants: config.svGenotypeBool = true; break;
        case UseDepth: config.useDepth = true; break;
        case Gpu: result.config.gpu.gpu = std::stoi(optarg); break;
        case Buffer: result.config.gpu.buffer = std::stoi(optarg); break;
        case 't': config.threads = std::max(std::stoi(optarg), 1); break;
        case 'D': config.debug = true; break;
        case 'h':
        case 'H':
            print_command_help(argv[0], result.command, gpu_enabled, option_code == 'H');
            return result;
        case '?':
            print_command_help(argv[0], result.command, gpu_enabled);
            result.status = 1;
            return result;
        }
    }

    // ------ Validation ------
    if (argc < 3) {
        print_command_help(argv[0], result.command, gpu_enabled);
        result.status = 1;
        return result;
    }
    const std::string error = validate_options(result.command, result.config, gpu_enabled);
    if (!error.empty()) {
        error_stream() << "Parameter error: " << error << "\n\n";
        print_command_help(argv[0], result.command, gpu_enabled);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}
