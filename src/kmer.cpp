// g++ -c kmer.cpp -std=c++17 -O3 -march=native

#include "../include/kmer.hpp"

using namespace std;

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;


/**
 * Find symmetric (k)-kmer on a DNA sequence (Counting Bloom Filter)
 *
 * @param str            DNA sequence
 * @param k              k-mer size
 * @param bf             Bloom filter of reference genome
 * 
 * @return 0
**/
int kmerBit::kmer_sketch_bf(const string& str, const uint32_t& k, BloomFilter* bf) {
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	unsigned int len = str.length();

	assert(len > 0 && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	for (i = l = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) // not an ambiguous base
		{
			int z;
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;

				bf->add(info.x);
			}
		} 
		else l = 0, kmer_span = 0;
	}

	return 0;
}


/**
 * Find symmetric (k)-kmer on a DNA sequence (graph construct)
 * Skip k-mer with a frequency greater than 1 in the reference genome
 *
 * @param str                       DNA sequence
 * @param k                         k-mer size
 * @param freKmerHashSetMap         output kmer hash table - map<frequence, unordered_set<kmerHash> >
 *                                  a[i].x = kMer<<8 | kmerSpan
 * @param bf                        Frequency of kmers in the reference genome, Counting Bloom Filter
**/
void kmerBit::kmer_sketch_construct(string& str, uint32_t k, map<uint8_t, unordered_set<uint64_t> >& freKmerHashSetMap, BloomFilter* bf) {
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	unsigned int len = str.length();

	assert(len > 0 && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	for (i = l = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) // not an ambiguous base
		{
			int z;
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;

				// k-mers with their frequence in bf
				auto& emplacedValue = freKmerHashSetMap.emplace(bf->count(info.x), unordered_set<uint64_t>()).first->second;
				emplacedValue.insert(info.x);
			}
		}
		else l = 0, kmer_span = 0;
	}
}


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
vector<uint64_t> kmerBit::kmer_sketch_fastq(
	const string& str, 
	const uint32_t& k, 
	const unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap
) {
	vector<uint64_t> hashVec;  // save k-mers
	hashVec.reserve(str.size());

	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	unsigned int len = str.length();

	assert(len > 0 && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	for (i = l = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) {  // not an ambiguous base
			int z;
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;

				if (GraphKmerHashHapStrMap.find(info.x) != GraphKmerHashHapStrMap.end()) {  // If there is a corresponding k-mer in the graph, add it
					hashVec.push_back(info.x);
				}
			}
		} 
		else l = 0, kmer_span = 0;
	}

	return hashVec;
}


/**
 * Find symmetric (k)-kmer on a DNA sequence (genotype)
 *
 * @param str                       DNA sequence
 * @param k                         k-mer size
 * 
 * @return hashSet                  output kmer hash table - set<kmerHash>
 *                                  a[i].x = kMer<<8 | kmerSpan
**/
unordered_set<uint64_t> kmerBit::kmer_sketch_genotype(
	const string& str, 
	const uint32_t& k
) {
	unordered_set<uint64_t> hashSet;
	hashSet.reserve(str.size());
	
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	unsigned int len = str.length();

	assert(len > 0 && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	for (i = l = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) // not an ambiguous base
		{
			int z;
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;

				// Records
				hashSet.insert(info.x);
			}
		}
		else l = 0, kmer_span = 0;
	}

	return hashSet;
}


/**
 * @author               zezhen du
 * @date                 2022/12/13
 * @version              v1.0
 * @brief                uncode hash function
 * @param hashNum        hash number of kmer by kmer_sketch

 * @return               make_tuple(chrId, kmerEnd, kmerStrand)
 *                       chrId
 *                       kmerEnd
 * 				   	     kmerStrand
**/
tuple<uint64_t, uint32_t, int> kmerBit::hash_uncode(const uint64_t & hashNum)
{
	uint64_t chrId = hashNum>>32;
	uint32_t kmerEnd = (hashNum - (chrId<<32))>>1;
	int kmerStrand = hashNum - (chrId<<32) - (kmerEnd<<1);

	return make_tuple(chrId, kmerEnd, kmerStrand);
}