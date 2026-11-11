#ifndef SAVE_HPP
#define SAVE_HPP

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>

#include "strip_split_join.hpp"
#include "logger.hpp"

using namespace std;

/**
 * @brief save result
 *
 * @param outputFileName
 *
 * @return
**/
class SAVE
{
private:
    // output stream
    string outputFileName_;

    // file handle
    ofstream fpO;
    gzFile gzfpO = nullptr;
public:
    SAVE() = default;
    explicit SAVE(string outputFileName);
    ~SAVE();

    SAVE(const SAVE&) = delete;
    SAVE& operator=(const SAVE&) = delete;

    int save(string & outTxt);
    // Finish the output and report buffered write errors before claiming success.
    void close();
};

#endif
