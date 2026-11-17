#include "../include/CommandRunner.hpp"

#include "../include/logger.hpp"

#ifndef VARIGRAPH_GPU
#include "../include/varigraph.hpp"
#endif

void print_config(Command command, const AppConfig& config) {
    const auto& common = config.common;
    stage_stream() << "Configuration ...\n";
    log_stream() << "  - Threads: " << common.threads << '\n';
    if (command == Command::Construct) {
        log_stream() << "  - Reference: " << common.refFileName << '\n';
        log_stream() << "  - Variants: " << common.vcfFileName << '\n';
        log_stream() << "  - Output graph: " << common.outputGraphFileName << '\n';
        log_stream() << "  - k-mer size: " << common.kmerLen << '\n';
        log_stream() << "  - VCF ploidy: " << common.vcfPloidy << '\n';
        log_stream() << "  - Fast mode: " << (common.fastMode ? "Enabled" : "Disabled") << '\n';
        log_stream() << "  - Unique k-mers: " << (common.useUniqueKmers ? "Enabled" : "Disabled") << '\n';
    } else {
        log_stream() << "  - Input graph: " << common.inputGraphFileName << '\n';
        log_stream() << "  - Sample configuration: " << common.samplesConfigFileName << '\n';
        log_stream() << "  - Sample genotype: " << common.sampleType << '\n';
        log_stream() << "  - Sample ploidy: " << common.samplePloidy << '\n';
        log_stream() << "  - Haplotypes: " << common.haploidNum << '\n';
        log_stream() << "  - Chromosome granularity: " << common.chrLenThread << " bp\n";
        log_stream() << "  - Transition probability: " << common.transitionProType << '\n';
        log_stream() << "  - SV only: " << (common.svGenotypeBool ? "Enabled" : "Disabled") << '\n';
        log_stream() << "  - Minimum GQ: " << common.minSupportingGQ << '\n';
        log_stream() << "  - Use sequencing depth: " << (common.useDepth ? "Enabled" : "Disabled") << '\n';
    }
    if (gpu_build) {
        log_stream() << "  - GPU ID: " << config.gpu.gpu << '\n';
        log_stream() << "  - GPU buffer: " << config.gpu.buffer << " MB\n";
    }
    log_stream() << "  - Debug: " << (common.debug ? "Enabled" : "Disabled") << '\n';

    if (command == Command::Genotype) {
        if (common.haploidNum < 10) {
            warning_stream() << "Parameter warning: -n. The number of haploids for genotyping is relatively low, which may affect the accuracy of genotyping.\n";
        }
        if (common.chrLenThread < 1000) {
            warning_stream() << "Parameter warning: --granularity. The chromosome granularity is less than 1000 bp (" << common.chrLenThread << " bp).\n";
        }
    }
    if (gpu_build && config.gpu.buffer > 1000) {
        warning_stream() << "Parameter warning: --buffer. The buffer size is relatively large, which may consume a large amount of GPU memory.\n";
    }
}

#ifndef VARIGRAPH_GPU
void run_command(Command command, const AppConfig& config) {
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
}
#endif
