//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

// Note - The x86 and x64 versions do _not_ produce the same results, as the
// algorithms are optimized for their respective platforms. You can still
// compile and run any of them on any platform, but your performance with the
// non-native version will be less than optimal.

#include "../include/MurmurHash3.cuh"

//-----------------------------------------------------------------------------
#define FORCE_INLINE inline __attribute__((always_inline))

__device__ uint64_t rotl64_kernel(uint64_t x, int8_t r) {
    return (x << r) | (x >> (64 - r));
}

#define BIG_CONSTANT_kernel(x) (x##LLU)
#define ROTL64_kernel(x,y)     rotl64_kernel(x,y)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

__device__ uint64_t getblock64_kernel(const uint64_t* p, int i) {
    return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

__device__ uint64_t fmix64_kernel ( uint64_t k ) {
  k ^= k >> 33;
  k *= BIG_CONSTANT_kernel(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT_kernel(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}

//-----------------------------------------------------------------------------

__global__ void MurmurHash3_x64_128_kernel (
  uint64_t * keys, 
  uint32_t num_keys, 
  int len,
  uint32_t seed,
  uint64_t * results
) {
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  // prevent out-of-bounds
  if (tid >= num_keys) return;

  const uint8_t *data = (const uint8_t*)(keys + tid);
  const int nblocks = len / 16;

  uint64_t h1 = seed;
  uint64_t h2 = seed;

  const uint64_t c1 = BIG_CONSTANT_kernel(0x87c37b91114253d5);
  const uint64_t c2 = BIG_CONSTANT_kernel(0x4cf5ad432745937f);

  //----------
  // body

  const uint64_t * blocks = (const uint64_t *)(data);

   for(int i = 0; i < nblocks; i++) {
    uint64_t k1 = getblock64_kernel(blocks,i*2+0);
    uint64_t k2 = getblock64_kernel(blocks,i*2+1);

    k1 *= c1; k1  = ROTL64_kernel(k1,31); k1 *= c2; h1 ^= k1;

    h1 = ROTL64_kernel(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

    k2 *= c2; k2  = ROTL64_kernel(k2,33); k2 *= c1; h2 ^= k2;

    h2 = ROTL64_kernel(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint64_t k1 = 0;
  uint64_t k2 = 0;

  switch(len & 15)
  {
  case 15: k2 ^= ((uint64_t)tail[14]) << 48;
  case 14: k2 ^= ((uint64_t)tail[13]) << 40;
  case 13: k2 ^= ((uint64_t)tail[12]) << 32;
  case 12: k2 ^= ((uint64_t)tail[11]) << 24;
  case 11: k2 ^= ((uint64_t)tail[10]) << 16;
  case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
  case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
           k2 *= c2; k2  = ROTL64_kernel(k2,33); k2 *= c1; h2 ^= k2;

  case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
  case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
  case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
  case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
  case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
  case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
  case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
  case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
           k1 *= c1; k1  = ROTL64_kernel(k1,31); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len; h2 ^= len;

  h1 += h2;
  h2 += h1;

  h1 = fmix64_kernel(h1);
  h2 = fmix64_kernel(h2);

  h1 += h2;
  h2 += h1;

  results[tid] = h1 + h2;
}


//-----------------------------------------------------------------------------