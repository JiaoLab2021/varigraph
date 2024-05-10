#ifndef CUDA_ERROR_HANDLING_H
#define CUDA_ERROR_HANDLING_H

#include <iostream>
#include <cuda_runtime.h>

#include "get_time.hpp"


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "GPUassert: " << cudaGetErrorString(code) << endl;
        if (abort) exit(code);
   }
}

#endif // CUDA_ERROR_HANDLING_H