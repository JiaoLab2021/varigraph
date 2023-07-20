// g++ -c kmer.cpp -std=c++17 -O3 -march=native

#include "../include/kmer.hpp"

using namespace std;

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint64_t kmerBit::hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}


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
void kmerBit::kmer_sketch(const uint32_t chromosomeId, const string & str, uint32_t k, unordered_map<uint64_t, vector<uint64_t> > & outMap)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

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
			if (l >= k && kmer_span < 256)
			{
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				if (z == 0) // 正向
				{
					info.y = (uint64_t)chromosomeId<<32 | (uint32_t)(i+1)<<1 | z;
				}
				else // 反向
				{
					info.y = (uint64_t)chromosomeId<<32 | (uint32_t)(len-i-1+k)<<1 | z;
				}
				outMap[info.x].push_back(info.y);
			}
		} 
		else l = 0, tq.count = tq.front = 0, kmer_span = 0;
	}
}


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
void kmerBit::kmer_sketch(const uint32_t chromosomeId, const string & str, uint32_t k, unordered_map<uint64_t, unordered_map<uint64_t, uint32_t> > & outMap)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

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
			if (l >= k && kmer_span < 256)
			{
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				// 不要终止位置和kmer方向信息，所以注释掉
				if (z == 0) // 正向
				{
					info.y = (uint64_t)chromosomeId<<32 | (uint32_t)(i+1)<<1 | z;
				}
				else // 反向
				{
					info.y = (uint64_t)chromosomeId<<32 | (uint32_t)(len-i-1+k)<<1 | z;
				}
				outMap[info.x][chromosomeId]++;
			}
		} 
		else l = 0, tq.count = tq.front = 0, kmer_span = 0;
	}
}


/**
 * Find symmetric (k)-kmer on a DNA sequence (reference genome)
 *
 * @param str            DNA sequence
 * @param k              k-mer size
 * @param bf             Bloom filter of reference genome
 * 
 * @return 0
**/
int kmerBit::kmer_sketch_fasta(const string& str, const uint32_t& k, BloomFilter* bf)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

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
		else l = 0, tq.count = tq.front = 0, kmer_span = 0;
	}

	return 0;
}


/**
 * Find symmetric (k)-kmer on a DNA sequence (graph)
 * 跳过在参考基因组中频率大于1的kmer
 *
 * @param str            DNA sequence
 * @param k              k-mer size
 * @param outMap         output kmer hash table - map<kmerHash, coverage>
 *                       a[i].x = kMer<<8 | kmerSpan
 * @param bf             参考基因组中kmer的频率, Counting Bloom Filter
**/
void kmerBit::kmer_sketch(string str, uint32_t k, unordered_map<uint64_t, uint8_t> & outMap, BloomFilter* bf)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

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
			if (l >= k && kmer_span < 256)
			{
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;

				if (bf->count(info.x) > 1) continue; // 在参考基因组中频率大于1的kmer跳过

				if (outMap[info.x] < UINT8_MAX)  // 防止变量越界
				{
					outMap[info.x]++;
				}
			}
		}
		else l = 0, tq.count = tq.front = 0, kmer_span = 0;
	}
}


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
vector<uint64_t> kmerBit::kmer_sketch_fastq(
	const string& str, 
	const uint32_t& k, 
	const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap
)
{
	vector<uint64_t> hashVec;  // save k-mers

	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

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

				if (GraphKmerCovFreMap.find(info.x) != GraphKmerCovFreMap.end()) {  // If there is a corresponding k-mer in the figure, add it
					hashVec.push_back(info.x);
				}
			}
		} 
		else l = 0, tq.count = tq.front = 0, kmer_span = 0;
	}

	return hashVec;
}


/**
 * @author               zezhen du
 * @date                 2022/12/13
 * @version              v1.0
 * @brief                uncode hash function
 * @param hashNum        hash number of kmer by kmer_sketch

 * @return            make_tuple(chrId, kmerEnd, kmerStrand)
 *                    chrId
 *                    kmerEnd
 * 					 kmerStrand
**/
tuple<uint64_t, uint32_t, int> kmerBit::hash_uncode(const uint64_t & hashNum)
{
	uint64_t chrId = hashNum>>32;
	uint32_t kmerEnd = (hashNum - (chrId<<32))>>1;
	int kmerStrand = hashNum - (chrId<<32) - (kmerEnd<<1);

	return make_tuple(chrId, kmerEnd, kmerStrand);
}