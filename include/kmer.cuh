#ifndef KMER_CUH
#define KMER_CUH

#include <map>
#include <unordered_set>

#include "counting_bloom_filter.cuh"
#include "seq_nt4_table.hpp"
#include "hash64.hpp"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/iterator/constant_iterator.h>

// k-mer (GPU version)
__global__ void kmer_sketch_kernel(const char* str, int len, uint32_t k, uint64_t* result);

/**
 * Find symmetric (k)-kmer on a DNA sequence (graph construct)
 * Skip k-mer with a frequency greater than 1 in the reference genome
 *
 * @param str                       DNA sequence
 * @param k                         k-mer size
 * @param freKmerHashSetMap         output kmer hash table - map<coverage, vector<kmerHash> >
 *                                  a[i].x = kMer<<8 | kmerSpan
 * @param bf                        Frequency of kmers in the reference genome, Counting Bloom Filter
**/
void kmer_sketch_construct_kernel(string& str, uint32_t k, map<uint8_t, unordered_set<uint64_t> >& freKmerHashSetMap, BloomFilterKernel* bf);

#endif // KMER_CUH