// g++ -c counting_bloom_filter.cpp -std=c++17 -O3 -march=native

#include "../include/counting_bloom_filter.hpp"

// overloaded assignment operator
BloomFilter& BloomFilter::operator=(const BloomFilter& other) {
    if (this != &other) {
        // copy member variable
        _size = other._size;
        _numHashes = other._numHashes;
        _seeds = other._seeds;

        // delete old memory
        if (_filter != nullptr) {
            delete[] _filter;
            _filter = nullptr;
        }

        // allocate new memory
        _filter = new uint8_t[_size];

        // copy data
        std::copy(other._filter, other._filter + _size, _filter);
    }
    return *this;
}

void BloomFilter::add(uint64_t kmer) {
    for (int i = 0; i < _numHashes; ++i) {
        uint64_t hashValue = _murmur_hash(static_cast<void*>(reinterpret_cast<uint64_t*>(&kmer)), sizeof(kmer), _seeds[i]);
        uint64_t position = hashValue % _size;
        if (_filter[position] < 255) { // Determine whether the counter has reached the maximum value
            _filter[position]++;
        }
    }
}


// Determine whether the k-mer exists in the bloom filter
bool BloomFilter::find(uint64_t kmer) {
    for (int i = 0; i < _numHashes; ++i) {
        uint64_t hashValue = _murmur_hash(static_cast<void*>(reinterpret_cast<uint64_t*>(&kmer)), sizeof(kmer), _seeds[i]);
        if (_filter[hashValue % _size] == 0) {
            return false;
        }
    }
    return true;
}

// Get the frequency of kmer
uint8_t BloomFilter::count(uint64_t kmer) {
    uint8_t minFreq = UINT8_MAX;
    int count = 0;
    for (int i = 0; i < _numHashes; ++i) {
        uint64_t hashValue = _murmur_hash(static_cast<void*>(reinterpret_cast<uint64_t*>(&kmer)), sizeof(kmer), _seeds[i]);
        uint64_t position = hashValue % _size;
        if (_filter[position] < minFreq) {
            minFreq = _filter[position];
            count = 1;
        } else if (_filter[position] == minFreq) {
            count++;
        }
    }

    // return (count == _numHashes) ? minFreq : 0;
    return minFreq;
}

// Calculate Bloom _filter size
uint64_t BloomFilter::_calculate_size(uint64_t n, double p) {
    return ceil((n * log(p)) / log(1.0 / pow(2.0, log(2.0))));
}

// Calculate the number of hash functions
uint32_t BloomFilter::_calculate_num_hashes(uint64_t n, uint64_t m) {
    return round(m * log(2.0) / n);
}

// Initialize the seed for the hash function
void BloomFilter::_init_seeds(uint32_t _numHashes) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<size_t> dis(1, numeric_limits<size_t>::max());
    for (int i = 0; i < _numHashes; ++i) {
        _seeds.push_back(dis(gen));
    }
}

// MurmurHash3 hash function
uint64_t BloomFilter::_murmur_hash(const void* key, int len, unsigned int seed) {
    // uint64_t hashValue;
    uint64_t hashValue[2];
    
    MurmurHash3_x64_128(key, len, seed, &hashValue);

    // return hashValue;
    return hashValue[0] + hashValue[1];
}

double BloomFilter::get_cap() {
    if (_size == 0) {
        // Return a special value indicating an invalid or undefined result
        return std::numeric_limits<double>::quiet_NaN();  // or return std::numeric_limits<double>::infinity();
    }

    // uint64_t totalNum = std::count_if(_filter.begin(), _filter.end(), [](const uint8_t& it) { return it > 0; });
    uint64_t totalNum = 0;
    for (uint64_t i = 0; i < _size; ++i) {
        if (_filter[i] > 0) {
            ++totalNum;
        }
    }

    return static_cast<double>(totalNum) / _size;
}

// clear
void BloomFilter::clear() {
    if (_filter != nullptr) {
        delete[] _filter;
        _filter = nullptr;
    }
    _filter = new uint8_t[_size]();
}

void BloomFilter::save(const string& filename) const {
    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom filter saved to file: " << filename << endl;

    std::ofstream file(filename, std::ios::binary);
    if (file) {
        // Holds the size of the bloom filter and the number of hash functions
        file.write(reinterpret_cast<const char*>(&_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&_numHashes), sizeof(uint32_t));

        // Save the seeds
        for (const auto& seed : _seeds) {
            file.write(reinterpret_cast<const char*>(&seed), sizeof(uint64_t));
        }
        
        // Save the content of the bloom filter
        for (uint64_t i = 0; i < _size; ++i) {
            file.write(reinterpret_cast<const char*>(&_filter[i]), sizeof(uint8_t));
        }
        
        file.close();
    } else {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << filename << "': No such file or directory." << endl;
        exit(1);
    }
}

void BloomFilter::load(const string& filename) {
    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom filter loaded from file: " << filename << endl;

    std::ifstream file(filename, std::ios::binary);
    if (file) {
        // Load size of bloom filter and number of hash functions
        file.read(reinterpret_cast<char*>(&_size), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&_numHashes), sizeof(uint32_t));
        
        // Empty the original Bloom filter
        if (_filter != nullptr) {
            delete[] _filter;
            _filter = nullptr;
        }
        _filter = new uint8_t[_size]();
        _seeds.clear();
        _seeds.resize(_numHashes);

        // Load the seeds
        for (auto& seed : _seeds) {
            file.read(reinterpret_cast<char*>(&seed), sizeof(uint64_t));
        }
        
        // Load the content of the bloom filter
        for (uint64_t i = 0; i < _size; ++i) {
            file.read(reinterpret_cast<char*>(&_filter[i]), sizeof(uint8_t));
        }
        
        file.close();
    } else {
        // log
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << filename << "': No such file or directory." << endl;
        exit(1);
    }
}
