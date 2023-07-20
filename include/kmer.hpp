#ifndef KMER_HPP
#define KMER_HPP

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <assert.h>
#include <tuple>

#include "counting_bloom_filter.hpp"
#include "construct_index.hpp"

using namespace std;

namespace kmerBit{
	uint64_t hash64(uint64_t key, uint64_t mask);


	/**
	 * Find symmetric (k)-kmer on a DNA sequence
	 *
	 * @param chromosomeId   reference ID
	 * @param str            DNA sequence
	 * @param k              k-mer size
	 * @param outMap         output kmer hash table - map<kmerHash, vec<chrId<<32 | end | strand>>
	 *                       a[i].x = kMer<<8 | kmerSpan
	 *                       a[i].y = rid<<32 | lastPos<<1 | strand
	 *                       where lastPos is the position of the last base of the i-th minimizer,
	 *                       and strand indicates whether the minimizer comes from the top or the bottom strand.
	 *                       Callers may want to set "p->n = 0"; otherwise results are appended to p
	**/
	void kmer_sketch(const uint32_t chromosomeId, const string & str, uint32_t k, unordered_map<uint64_t, vector<uint64_t> > & outMap);


	/**
	 * Find symmetric (k)-kmer on a DNA sequence
	 *
	 * @param chromosomeId   reference ID
	 * @param str            DNA sequence
	 * @param k              k-mer size
	 * @param outMap         output kmer hash table - map<kmerHash, map<chrId, kmerNum>>
	 *                       a[i].x = kMer<<8 | kmerSpan
	 *                       p->a[i].y = chromosomeId<<32 | lastPos<<1 | strand
	 * 						 outMap map<kmerScore, map<readId, kmerNum>>
	 *                       where lastPos is the position of the last base of the i-th minimizer,
	 *                       and strand indicates whether the minimizer comes from the top or the bottom strand.
	 *                       Callers may want to set "p->n = 0"; otherwise results are appended to p
	**/
	void kmer_sketch(const uint32_t chromosomeId, const string & str, uint32_t k, unordered_map<uint64_t, unordered_map<uint64_t, uint32_t> > & outMap);


	/**
	 * Find symmetric (k)-kmer on a DNA sequence (reference genome)
	 *
	 * @param str            DNA sequence
	 * @param k              k-mer size
	 * @param bf             Bloom filter of reference genome
	 * 
	 * @return 0
	**/
	int kmer_sketch_fasta(const string& str, const uint32_t& k, BloomFilter* bf);


	/**
	 * Find symmetric (k)-kmer on a DNA sequence (graph)
	 * 只添加在参考基因组中频率为1的kmer
	 *
	 * @param str            DNA sequence
	 * @param k              k-mer size
	 * @param outMap         output kmer hash table - map<kmerHash, coverage>
	 *                       a[i].x = kMer<<8 | kmerSpan
	 * @param bf             参考基因组中kmer的频率, Counting Bloom Filter
	**/
	void kmer_sketch(string str, uint32_t k, unordered_map<uint64_t, uint8_t> & outMap, BloomFilter* bf);
	

	/**
	 * Find symmetric (k)-kmer on a DNA sequence (read)
	 *
	 * @param str                   DNA sequence
	 * @param k                     k-mer size
	 * @param GraphKmerCovFreMap    Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFre>
	 * 
	 * @return hashVec              output kmer hash table - vector<kmerHash>
	 *                              a[i].x = kMer<<8 | kmerSpan
	**/
	vector<uint64_t> kmer_sketch_fastq(
		const string& str, 
		const uint32_t& k, 
		const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap
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
	tuple<uint64_t, uint32_t, int> hash_uncode(const uint64_t & hashNum);
	
}

#endif // KMER_HPP