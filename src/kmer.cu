#include "../include/kmer.cuh"


// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;


__device__ unsigned char seq_nt4_table_kernel[256] = {
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

__device__ uint64_t hash64_kernel(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

// k-mer (GPU version)
__global__ void kmer_sketch_kernel(const char* str, int len, uint32_t k, uint64_t* result) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	// prevent out-of-bounds
    if (tid >= len - k + 1) return;

    uint32_t start = tid, end = tid + k;

    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
    uint32_t l = 0;
    int kmer_span = 0;

    for (uint32_t i = start; i < end; i++) {
        int c = seq_nt4_table_kernel[(uint8_t)str[i]];
        if (c < 4) {
            int z;
            kmer_span = l + 1 < k ? l + 1 : k;
            kmer[0] = (kmer[0] << 2 | c) & mask;
            kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1;
            if (kmer[0] == kmer[1]) continue;
            z = kmer[0] < kmer[1] ? 0 : 1;
            ++l;
            if (l >= k && kmer_span < 256) {
                result[tid] = hash64_kernel(kmer[z], mask) << 8 | kmer_span;
            }
        } else {
            l = 0;
            kmer_span = 0;
        }
    }
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
void kmer_sketch_construct_kernel(string& str, uint32_t k, map<uint8_t, unordered_set<uint64_t> >& freKmerHashSetMap, BloomFilterKernel* bf) {
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, kmer_span = 0;

	unsigned int len = str.length();

	assert(len > 0 && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	for (i = l = 0; i < len; ++i) {
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