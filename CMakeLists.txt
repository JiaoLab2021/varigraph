cmake_minimum_required(VERSION 3.12)

project(varigraph CUDA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add source files for CPU version
set(SRC_FILES_CPU
    src/construct_index.cpp
    src/counting_bloom_filter.cpp
    src/fastq_kmer.cpp
    src/genotype.cpp
    src/get_time.cpp
    src/haplotype_select.cpp
    src/kmer.cpp
    src/MurmurHash3.cpp
    src/save.cpp
    src/strip_split_join.cpp
    src/varigraph.cpp
)

# Add source files for GPU version
set(SRC_FILES_GPU
    src/construct_index.cu
    src/counting_bloom_filter.cu
    src/fastq_kmer.cu
    src/kmer.cu
    src/MurmurHash3.cu
    src/varigraph.cu
)

# Add executable for CPU version
add_executable(varigraph main.cpp ${SRC_FILES_CPU})
target_compile_options(varigraph PRIVATE -O3)

# Link libraries for CPU version
target_include_directories(varigraph PRIVATE include)
find_package(Threads REQUIRED)
target_link_libraries(varigraph Threads::Threads)
find_package(ZLIB REQUIRED)
target_link_libraries(varigraph ZLIB::ZLIB)
target_link_libraries(varigraph stdc++fs)

# Check for GPU, CUDA, and NVCC
cmake_policy(SET CMP0104 NEW)
find_package(CUDA)
if(CUDA_FOUND)
    message(STATUS "Found CUDA, building GPU version.")
    # Set compiler options for CUDA
    set(CUDA_NVCC_FLAGS -O3 -std=c++17 --extended-lambda)
    
    # Add executable for GPU version
    cuda_add_executable(varigraph-gpu main.cu ${SRC_FILES_CPU} ${SRC_FILES_GPU})

    # Set CUDA_ARCHITECTURES
    set_property(TARGET varigraph-gpu PROPERTY CUDA_ARCHITECTURES "sm_50")

    # Include directories and link libraries
    target_include_directories(varigraph-gpu PRIVATE include)
    target_link_libraries(varigraph-gpu Threads::Threads)
    target_link_libraries(varigraph-gpu ZLIB::ZLIB)
    target_link_libraries(varigraph-gpu stdc++fs)

    # Install varigraph-gpu
    install(TARGETS varigraph-gpu DESTINATION bin)
else()
    message(STATUS "CUDA not found, building CPU version only.")
endif()

# Install function for CPU version
install(TARGETS varigraph DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
file(GLOB EXECUTABLES "bin/varigraph*")
foreach(EXECUTABLE ${EXECUTABLES})
    install(PROGRAMS ${EXECUTABLE} DESTINATION bin)
endforeach()
