#include "../include/CommandRunner.hpp"
#include "../include/varigraph.cuh"

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include "../include/logger.hpp"

// Select the device and run the GPU workflow with the caller's configuration.
int run_command(Command command, const AppConfig& config) {
    stage_stream() << "Selecting GPU ...\n";
    const auto selection_status = cudaSetDevice(config.gpu.gpu);
    if (selection_status != cudaSuccess) {
        throw std::runtime_error("Cannot select GPU " + std::to_string(config.gpu.gpu) + ": " + cudaGetErrorString(selection_status));
    }
    cudaDeviceProp device{};
    const auto properties_status = cudaGetDeviceProperties(&device, config.gpu.gpu);
    if (properties_status != cudaSuccess) {
        throw std::runtime_error("Cannot read GPU " + std::to_string(config.gpu.gpu) + " properties: " + cudaGetErrorString(properties_status));
    }
    log_stream() << "  - Selected GPU: " << device.name << '\n';
    log_stream() << "  - Compute capability: " << device.major << '.' << device.minor << '\n';
    log_stream() << "  - GPU memory: " << device.totalGlobalMem / 1024 / 1024 / 1024 << " GB\n";

    {
        VarigraphKernel app(config.common, config.gpu);
        if (command == Command::Construct) {
            app.construct_kernel();
        } else {
            app.parse_sample_config();
            app.load();
            app.fastq_genotype_kernel();
        }
        stage_stream() << "Releasing graph memory ...\n";
    }
    return 0;
}
