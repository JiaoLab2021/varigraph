//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

#ifndef _MURMURHASH3_CUH
#define _MURMURHASH3_CUH

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

// Microsoft Visual Studio

#if defined(_MSC_VER) && (_MSC_VER < 1600)

typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;

// Other compilers

#else	// defined(_MSC_VER)

#include <cuda_runtime.h>
#include <thrust/device_vector.h>

#endif // !defined(_MSC_VER)

//-----------------------------------------------------------------------------

__global__ void MurmurHash3_x64_128_kernel(uint64_t * keys, uint32_t num_keys, int len, uint32_t seed, uint64_t * results);

//-----------------------------------------------------------------------------

#endif // _MURMURHASH3_H
