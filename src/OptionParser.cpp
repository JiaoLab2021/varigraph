#include "../include/OptionParser.hpp"

#include <algorithm>
#include <cmath>
#include <getopt.h>
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
    help.line("-H, --advanced", "", "show all help");
    help.line("--version", "", "show the program version");
    help.print(diagnostic_stream());
}

// Print command options using the configuration defaults.
void print_command_help(const char* executable, Command command, bool advanced = false) {
    const AppConfig defaults;
    const auto& common = defaults.common;
    HelpPrinter help;

    if (command == Command::Construct) {
        diagnostic_stream() << "Usage: " << executable << " construct -r FILE -v FILE [options]\n\n"
                  << "Build a genome graph index from a reference FASTA and population VCF.\n\n";
        help.section("Input/output");
        help.line("-r, --reference", "FILE", "input reference FASTA");
        help.line("-v, --vcf", "FILE", "population VCF for index building");
        help.line("--save-graph", "FILE", "save the graph index [" + common.outputGraphFileName + "]");
        help.blank();
        help.section("Ploidy");
        help.line("--vcf-ploidy", "INT", "ploidy of genotypes in the VCF, 2-8 [" + std::to_string(common.vcfPloidy) + "]");
        help.blank();
        help.section("Indexing");
        help.line("-k, --kmer", "INT", "k-mer size, values below 5 use 5; maximum 28 [" + std::to_string(common.kmerLen) + "]");
        if (advanced) {
            help.line("--fast", "", "build faster with a possible reduction in genotyping accuracy");
            help.line("--use-unique-kmers", "", "use only unique k-mers for indexing");
        }
    } else {
        diagnostic_stream() << "Usage: " << executable << " genotype --load-graph FILE -s FILE [options]\n\n"
                  << "Genotype samples from their reads using a saved graph index.\n\n";
        help.section("Input/output");
        help.line("--load-graph", "FILE", "load a graph index made by construct [" + common.inputGraphFileName + "]");
        help.line("-s, --sample", "FILE", "sample configuration; each line: SAMPLE READS ...");
        help.blank();
        help.section("Sample");
        help.line("-g, --genotype", "STR", "sample genotype: hom (homozygous) or het (heterozygous) [" + common.sampleType + "]");
        help.line("--sample-ploidy", "INT", "sample ploidy, 2-8 [" + std::to_string(common.samplePloidy) + "]");
        help.blank();
        help.section("Genotyping");
        if (advanced) {
            std::ostringstream granularity;
            granularity << common.chrLenThread / 1e6;
            std::ostringstream min_support;
            min_support << common.minSupportingGQ;
            help.line("-n, --number", "INT", "number of haplotypes used for genotyping [" + std::to_string(common.haploidNum) + "]");
            help.line("--granularity", "FLOAT", "chromosome length processed per thread, in Mb [" + granularity.str() + "]");
            help.line("-m, --mode", "STR", "transition probability: fre (haplotype frequency) or rec (recombination rate) [" + common.transitionProType + "]");
            help.line("--min-support", "FLOAT", "minimum genotype quality (GQ) [" + min_support.str() + "]");
        }
        help.line("--sv", "", "genotype structural variants only");
        help.line("--use-depth", "", "use sequencing depth as the depth of homozygous k-mers");
    }

    if (gpu_build) {
        help.blank();
        help.section("GPU");
        help.line("--gpu", "INT", "GPU device ID [" + std::to_string(defaults.gpu.gpu) + "]");
        help.line("--buffer", "INT", "GPU buffer size in MB, maximum " + std::to_string(max_gpu_buffer_mb) + " [" + std::to_string(defaults.gpu.buffer) + "]");
    }

    help.blank();
    help.section("General options");
    help.line("-t, --threads", "INT", "number of compute threads [" + std::to_string(common.threads) + "]");
    if (advanced) help.line("-D, --debug", "", "enable debug output");
    help.line("-h, --help", "", "show common help");
    help.line("-H, --advanced", "", "show all help");
    help.print(diagnostic_stream());
}

std::string option_name(const std::vector<option>& options, int code) {
    if (code < 256) return std::string("-") + static_cast<char>(code);
    for (const auto& entry : options) {
        if (entry.name != nullptr && entry.val == code) return std::string("--") + entry.name;
    }
    return {};
}

std::string validate_options(Command command, const AppConfig& config) {
    const auto& common = config.common;

    if (command == Command::Construct) {
        if (common.refFileName.empty()) return "-r. The reference genome file cannot be empty.";
        if (common.vcfFileName.empty()) return "-v. The VCF file cannot be empty.";
        if (common.outputGraphFileName.empty()) return "--save-graph. The genome graph file cannot be empty.";
        if (common.vcfPloidy > 8) {
            return "--vcf-ploidy. The provided value must be between 2 and 8 (inclusive).";
        }
        if (common.kmerLen > 28) {
            return "-k. The k-mer size must not exceed 28; values below 5 use 5.";
        }
        if (gpu_build && config.gpu.gpu < 0) {
            return "--gpu. The provided value must be greater than or equal to 0.";
        }
    } else {
        if (common.inputGraphFileName.empty()) return "--load-graph. The genome graph file cannot be empty.";
        if (common.samplesConfigFileName.empty()) return "-s. The sample configuration file cannot be empty.";
        if (common.sampleType != "hom" && common.sampleType != "het") {
            return "-g. The provided value must be either 'hom' or 'het'.";
        }
        if (common.samplePloidy > 8) {
            return "--sample-ploidy. The provided value must be between 2 and 8 (inclusive).";
        }
        if (common.haploidNum == 0) return "-n. The provided value must be greater than 0.";
        if (common.chrLenThread < 1) {
            return "--granularity. The chromosome granularity must be at least 1 bp (" + std::to_string(common.chrLenThread) + ").";
        }
        if (common.transitionProType != "fre" && common.transitionProType != "rec") {
            return "-m. The transition probability type must be either 'fre' or 'rec'.";
        }
    }

    if (gpu_build) {
        if (config.gpu.buffer <= 0) return "--buffer. The provided value must be greater than 0.";
        if (config.gpu.buffer > max_gpu_buffer_mb) {
            return "--buffer. The provided value must not exceed " + std::to_string(max_gpu_buffer_mb) + " MB.";
        }
    }
    return {};
}

}

ParseResult parse_command(int argc, char** argv) {
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
        error_stream() << "Command not found: " << command << "\n\n";
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
    if (gpu_build) {
        long_options.push_back({"gpu", required_argument, nullptr, Gpu});
        long_options.push_back({"buffer", required_argument, nullptr, Buffer});
    }
    long_options.push_back({"threads", required_argument, nullptr, 't'});
    long_options.push_back({"debug", no_argument, nullptr, 'D'});
    long_options.push_back({"help", no_argument, nullptr, 'h'});
    long_options.push_back({"advanced", no_argument, nullptr, 'H'});
    long_options.push_back({nullptr, 0, nullptr, 0});

    const char* short_options = result.command == Command::Construct ? "r:v:k:t:DhH" : "s:g:n:m:t:DhH";
    auto& common = result.config.common;
    optind = 0;
    int option_code;
    while ((option_code = getopt_long(argc, argv, short_options, long_options.data(), nullptr)) != -1) {
        try {
            switch (option_code) {
            case 'r': common.refFileName = optarg; break;
            case 'v': common.vcfFileName = optarg; break;
            case SaveGraph: common.outputGraphFileName = optarg; break;
            case LoadGraph: common.inputGraphFileName = optarg; break;
            case 's': common.samplesConfigFileName = optarg; break;
            case 'g': common.sampleType = optarg; break;
            case VcfPloidy: common.vcfPloidy = std::max(std::stoi(optarg), 2); break;
            case SamplePloidy: common.samplePloidy = std::max(std::stoi(optarg), 2); break;
            case 'k': common.kmerLen = std::max(std::stoi(optarg), 5); break;
            case Fast: common.fastMode = true; break;
            case UseUniqueKmers: common.useUniqueKmers = true; break;
            case 'n': common.haploidNum = std::stoull(optarg); break;
            case Granularity: {
                const double bases = std::stof(optarg) * 1e6;
                if (!std::isfinite(bases) || bases < 0.0 ||
                    bases >= static_cast<double>(std::numeric_limits<std::uint32_t>::max()) + 1.0) {
                    throw std::out_of_range("granularity");
                }
                common.chrLenThread = static_cast<std::uint32_t>(bases);
                break;
            }
            case 'm': common.transitionProType = optarg; break;
            case MinSupport: common.minSupportingGQ = std::stof(optarg); break;
            case StructuralVariants: common.svGenotypeBool = true; break;
            case UseDepth: common.useDepth = true; break;
            case Gpu: result.config.gpu.gpu = std::stoi(optarg); break;
            case Buffer: result.config.gpu.buffer = std::stoi(optarg); break;
            case 't': common.threads = std::max(std::stoi(optarg), 1); break;
            case 'D': common.debug = true; break;
            case 'h':
            case 'H':
                print_command_help(argv[0], result.command, option_code == 'H');
                return result;
            case '?':
                print_command_help(argv[0], result.command);
                result.status = 1;
                return result;
            }
        } catch (const std::logic_error&) {  // invalid or out-of-range numeric value
            error_stream() << "Parameter error: " << option_name(long_options, option_code)
                           << ". Invalid value '" << optarg << "'.\n\n";
            print_command_help(argv[0], result.command);
            result.status = 1;
            return result;
        }
    }

    // ------ Validation ------
    if (argc < 3) {
        print_command_help(argv[0], result.command);
        result.status = 1;
        return result;
    }
    const std::string error = validate_options(result.command, result.config);
    if (!error.empty()) {
        error_stream() << "Parameter error: " << error << "\n\n";
        print_command_help(argv[0], result.command);
        result.status = 1;
        return result;
    }

    result.run = true;
    return result;
}
