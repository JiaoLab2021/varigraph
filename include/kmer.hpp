#ifndef KMER_HPP
#define KMER_HPP

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include <tuple>

#include "counting_bloom_filter.hpp"
#include "construct_index.hpp"

using namespace std;

namespace kmerBit{
	uint64_t hash64(uint64_t key, uint64_t mask);


	/**
	 * Find symmetric (k)-kmer on a DNA sequence (Counting Bloom Filter)
	 *
	 * @param str            DNA sequence
	 * @param k              k-mer size
	 * @param bf             Bloom filter of reference genome
	 * 
	 * @return 0
	**/
	int kmer_sketch_bf(const string& str, const uint32_t& k, BloomFilter* bf);


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
	void kmer_sketch_construct(string str, uint32_t k, map<uint8_t, unordered_set<uint64_t> >& freKmerHashSetMap, BloomFilter* bf);
	

	/**
	 * Find symmetric (k)-kmer on a DNA sequence (read)
	 *
	 * @param str                       DNA sequence
	 * @param k                         k-mer size
	 * @param GraphKmerHashHapStrMap    Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVec>
	 * 
	 * @return hashVec                  output kmer hash table - vector<kmerHash>
	 *                                  a[i].x = kMer<<8 | kmerSpan
	**/
	vector<uint64_t> kmer_sketch_fastq(
		const string& str, 
		const uint32_t& k, 
		const unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap
	);


	/**
	 * Find symmetric (k)-kmer on a DNA sequence (genotype)
	 *
	 * @param str                       DNA sequence
	 * @param k                         k-mer size
	 * 
	 * @return hashSet                  output kmer hash table - set<kmerHash>
	 *                                  a[i].x = kMer<<8 | kmerSpan
	**/
	unordered_set<uint64_t> kmer_sketch_genotype(
		const string& str, 
		const uint32_t& k
	);


	/**
     * @author               zezhen du
     * @date                 2022/12/13
     * @version              v1.0
	 * @brief                uncode hash function
	 * @param hashNum        hash number of kmer by kmer_sketch

	 * @return               make_tuple(chrId, kmerEnd, kmerStrand)
	 *                       chrId
	 *                       kmerEnd
	 * 						 kmerStrand
	**/
	tuple<uint64_t, uint32_t, int> hash_uncode(const uint64_t& hashNum);
	
}

#endif // KMER_HPP