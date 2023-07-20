#ifndef COUNTING_BLOOM_FILTER_HPP
#define COUNTING_BLOOM_FILTER_HPP

#include <iostream>
#include <bitset>
#include <random>
#include <vector>
#include <climits>
#include <algorithm>
#include <fstream>

#include "MurmurHash3.h"
#include "get_time.hpp"

using namespace std;

/*
'Less Hashing, Same Performance: Building a Better Bloom Filter'
https://hur.st/bloomfilter
m = ceil(n * log(p) / log(1 / (pow(2, log(2)))));
k = round(m * log(2) / n);

Among them, n represents the number of elements to be stored, and p represents the allowable false positive rate. 
These two parameters are used to determine the size of the Bloom filter and the number of hash functions.
*/

// Define the BloomFilter class
class BloomFilter
{
public:
    /**
     * @author zezhen du
     * @date 2023/06/23
     * @version v1.0
     * @brief counting bloom _filter
     * 
     * @param n          the estimated number of elements that the Bloom _filter will need to store.
     * @param p          the estimated error rate of the Bloom _filter
     * 
     * @return 0
    **/

    BloomFilter() {}
    BloomFilter(uint64_t n, double p) : _size(_calculate_size(n, p)), _numHashes(_calculate_num_hashes(n, _size)), _filter(_size, 0) {
        _init_seeds(_numHashes);
    }

    ~BloomFilter() {
        _filter.clear();
        vector<uint8_t>().swap(_filter);
    }

    // overloaded assignment operator
    BloomFilter& operator=(const BloomFilter& other);

    // Calculate the hash value and set the counter at the corresponding position plus one
    void add(uint64_t kmer);

    // Determine whether the k-mer exists in the bloom filter
    bool find(uint64_t kmer);

    // Get the frequency of kmer
    uint8_t count(uint64_t kmer);

    // get_state
    uint64_t get_size();
    uint32_t get_num();
    double get_cap();

    // clear
    void clear();

    // save or load _filter
    void save(const string& filename) const;
    void load(const string& filename);

private:
    uint64_t _size;  // Bloom filter size
    uint32_t _numHashes;  // number of hash functions
    vector<uint64_t> _seeds;  // Stores the seed of the hash function
    vector<uint8_t> _filter;  // Storing Bloom Filters

    // Calculate Bloom _filter size
    static uint64_t _calculate_size(uint64_t n, double p);

    // Calculate the number of hash functions
    static uint32_t _calculate_num_hashes(uint64_t n, uint64_t m);

    // Initialize the seed for the hash function
    void _init_seeds(uint32_t _numHashes);

    // MurmurHash3 hash function
    uint64_t _murmur_hash(const void* key, int len, unsigned int seed);
};

#endif