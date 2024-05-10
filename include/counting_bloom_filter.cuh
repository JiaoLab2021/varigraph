#ifndef COUNTING_BLOOM_FILTER_CUH
#define COUNTING_BLOOM_FILTER_CUH

#include "counting_bloom_filter.hpp"
#include "get_time.hpp"
#include "cuda_error_handling.hpp"

#include "MurmurHash3.cuh"
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/iterator/constant_iterator.h>

using namespace std;


// Define the capped_add struct
struct capped_add {
    __device__ uint8_t operator()(const uint8_t& a, const uint8_t& b) const {
        uint16_t sum = static_cast<uint16_t>(a) + static_cast<uint16_t>(b);
        return sum > 255 ? 255 : static_cast<uint8_t>(sum);
    }
};


// This function performs an atomic addition operation on a uint8_t value, ensuring that the result does not exceed 255.
__device__ uint8_t atomicAddUint8(uint8_t* address, uint8_t val);

// This kernel function performs counting for a bloom filter. It takes an array of unique hashes, their counts, and a filter array.
// Each thread handles one unique hash. It calculates the position in the filter array by taking the modulo of the hash value with the size of the filter.
// Then it atomically adds the count of the hash to the value at the calculated position in the filter array.
__global__ void add_filter_kernel(uint64_t* uniqueHashsD, uint8_t* uniqueHashCountsD, uint32_t numUniqueHashs, uint8_t* filter, uint64_t size);

// Define the BloomFilterKernel class
class BloomFilterKernel : public BloomFilter {
public:
    // Constructor
    BloomFilterKernel() : BloomFilter() {}
    BloomFilterKernel(uint64_t size, double errorRate) : BloomFilter(size, errorRate) {
        gpuErrchk(cudaMallocManaged(&_filterD, sizeof(uint8_t) * _size));  // Allocate unified memory to store the Counting Bloom Filter
    }
    // Destructor
    ~BloomFilterKernel() {
        if (_filterD != nullptr) {
            gpuErrchk(cudaFree(_filterD));  // Free the memory occupied by the Counting Bloom Filter
            _filterD = nullptr;
        }
    }

    /*
     * Variable: kmersD
     * -----------------
     * This variable is declared on the device (GPU) and is only accessible 
     * from the device code (i.e., CUDA kernels). It cannot be directly accessed 
     * from the host (CPU) code.
    */
    void add_kernel(uint64_t * kmersD, uint32_t numKmers);
    void find_kernel(uint64_t * kmersD, uint32_t numKmers, bool * findBoolsD);
    void count_kernel(uint64_t * kmersD, uint32_t numKmers, uint8_t * countsD);

    /*
     * Function: copyFilterDToHost
     * ------------------------------
     * This function copies the content of the Counting Bloom Filter from the device (GPU) to the host (CPU).
     *
     * It first ensures that all previous CUDA operations are completed by calling cudaDeviceSynchronize.
     * Then, it uses cudaMemcpy to copy the content of _filterD (on the device) to _filter (on the host).
     * The cudaMemcpy function is called with cudaMemcpyDeviceToHost, indicating that the copy is from the device to the host.
    */
    void copyFilterDToHost() {
        // Ensure that all previous CUDA operations are completed
        gpuErrchk(cudaDeviceSynchronize());

        // Copy the content of _filterD from the device to the host
        gpuErrchk(cudaMemcpy(_filter, _filterD, _size * sizeof(uint8_t), cudaMemcpyDeviceToHost));

        // Print a message indicating that the copy operation is complete
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom Filter copied from device to host ...\n";

        // Free the memory occupied by the Counting Bloom Filter on the device
        gpuErrchk(cudaFree(_filterD));
        _filterD = nullptr;
    }
private:
    // Counting Bloom Filter
    uint8_t* _filterD = nullptr;
};

#endif