#include "../include/CommandRunner.hpp"

#include "../include/logger.hpp"

#ifndef VARIGRAPH_GPU
#include "../include/varigraph.hpp"
#endif

// Print the effective CLI configuration before executing a command.
void print_config(Command command, const AppConfig& config, bool gpu_enabled) {
    const auto& options = config.common;
    stage_stream() << "Configuration ...\n";
    log_stream() << "  - Threads: " << options.threads << '\n';
    if (command == Command::Construct) {
        log_stream() << "  - Reference: " << options.refFileName << '\n';
        log_stream() << "  - Variants: " << options.vcfFileName << '\n';
        log_stream() << "  - Output graph: " << options.outputGraphFileName << '\n';
        log_stream() << "  - k-mer size: " << options.kmerLen << '\n';
        log_stream() << "  - VCF ploidy: " << options.vcfPloidy << '\n';
        log_stream() << "  - Fast mode: " << (options.fastMode ? "Enabled" : "Disabled") << '\n';
        log_stream() << "  - Unique k-mers: " << (options.useUniqueKmers ? "Enabled" : "Disabled") << '\n';
    } else {
        log_stream() << "  - Input graph: " << options.inputGraphFileName << '\n';
        log_stream() << "  - Sample configuration: " << options.samplesConfigFileName << '\n';
        log_stream() << "  - Sample genotype: " << options.sampleType << '\n';
        log_stream() << "  - Sample ploidy: " << options.samplePloidy << '\n';
        log_stream() << "  - Haplotypes: " << options.haploidNum << '\n';
        log_stream() << "  - Chromosome granularity: " << options.chrLenThread << " bp\n";
        log_stream() << "  - Transition probability: " << options.transitionProType << '\n';
        log_stream() << "  - SV only: " << (options.svGenotypeBool ? "Enabled" : "Disabled") << '\n';
        log_stream() << "  - Minimum GQ: " << options.minSupportingGQ << '\n';
        log_stream() << "  - Use sequencing depth: " << (options.useDepth ? "Enabled" : "Disabled") << '\n';
    }
    if (gpu_enabled) {
        log_stream() << "  - GPU ID: " << config.gpu.gpu << '\n';
        log_stream() << "  - GPU buffer: " << config.gpu.buffer << " MB\n";
    }
    log_stream() << "  - Debug: " << (options.debug ? "Enabled" : "Disabled") << '\n';

    if (command == Command::Genotype) {
        if (options.haploidNum < 10) {
            warning_stream() << "Parameter warning: -n. The number of haploids for genotyping is relatively low, which may affect the accuracy of genotyping.\n";
        }
        if (options.chrLenThread < 1000) {
            warning_stream() << "Parameter warning: --granularity. The chromosome granularity is less than 1000 bp ("
                             << options.chrLenThread << " bp).\n";
        }
    }
    if (gpu_enabled && config.gpu.buffer > 1000) {
        warning_stream() << "Parameter warning: --buffer. The buffer size is relatively large, which may consume a large amount of GPU memory.\n";
    }
}

#ifndef VARIGRAPH_GPU
// Run the CPU workflow while the caller owns its configuration.
int run_command(Command command, const AppConfig& config) {
    {
        Varigraph app(config.common);
        if (command == Command::Construct) {
            app.construct();
        } else {
            app.parse_sample_config();
            app.load();
            app.fastq_genotype();
        }
        stage_stream() << "Releasing graph memory ...\n";
    }
    return 0;
}
#endif
