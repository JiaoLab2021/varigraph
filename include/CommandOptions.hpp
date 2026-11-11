#pragma once

#include <cstdint>
#include <string>

enum class Command { Construct, Genotype };

struct VarigraphConfig {
    std::string refFileName;
    std::string vcfFileName;
    std::string samplesConfigFileName;
    std::string inputGraphFileName = "graph.bin";
    std::string outputGraphFileName = "graph.bin";
    bool fastMode = false;
    std::uint32_t kmerLen = 27;
    std::string sampleType = "het";
    std::uint32_t samplePloidy = 2;
    std::uint32_t vcfPloidy = 2;
    std::uint32_t haploidNum = 15;
    std::uint32_t chrLenThread = 1'000'000;
    std::string transitionProType = "rec";
    bool svGenotypeBool = false;
    bool debug = false;
    std::uint32_t threads = 10;
    float minSupportingGQ = 0.0f;
    bool useUniqueKmers = false;
    bool useDepth = false;
};

struct GpuOptions {
    int gpu = 0;
    int buffer = 100;
};

struct AppConfig {
    VarigraphConfig common;
    GpuOptions gpu;
};
