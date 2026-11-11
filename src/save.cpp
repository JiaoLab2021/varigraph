#include "../include/save.hpp"
#include "../include/logger.hpp"

#include <climits>
#include <stdexcept>
#include <utility>

using namespace std;

SAVE::SAVE(string outputFileName) : outputFileName_(std::move(outputFileName)) {
    if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos) {
        gzfpO = gzopen(outputFileName_.c_str(), "wb");
        if (!gzfpO) throw runtime_error("Cannot open output '" + outputFileName_ + "'.");
    } else if (!outputFileName_.empty()) {
        fpO.open(outputFileName_, ios::out);
        if (!fpO) throw runtime_error("Cannot open output '" + outputFileName_ + "'.");
    }
}

// Normal callers close explicitly; unwinding still releases an open gzip handle.
SAVE::~SAVE() {
    if (gzfpO) gzclose(gzfpO);
}

void SAVE::close() {
    if (gzfpO) {
        const int status = gzclose(std::exchange(gzfpO, nullptr));
        if (status != Z_OK) throw runtime_error("Cannot finish gzip output '" + outputFileName_ + "'.");
    } else if (fpO.is_open()) {
        fpO.close();
        if (!fpO) throw runtime_error("Cannot finish output '" + outputFileName_ + "'.");
    } else if (outputFileName_.empty()) {
        output_stream().flush();
        if (!output_stream()) throw runtime_error("Cannot finish standard output.");
    }
}

/**
 * @brief save result to file
 *
 * @return int
**/
int SAVE::save(
    string & outTxt
)
{
    outTxt = strip(outTxt, '\n');  // remove '\n'

    if (outTxt.empty()) {
        return 0;
    }

    if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos) {
        outTxt += "\n";
        if (!gzfpO) throw runtime_error("Output is already closed: '" + outputFileName_ + "'.");
        for (size_t offset = 0; offset < outTxt.size();) {
            const auto count = static_cast<unsigned int>(std::min(outTxt.size() - offset, static_cast<size_t>(INT_MAX)));
            if (gzwrite(gzfpO, outTxt.data() + offset, count) != static_cast<int>(count)) {
                throw runtime_error("Cannot write gzip output '" + outputFileName_ + "'.");
            }
            offset += count;
        }
    } else if (outputFileName_.size() > 0) {
        fpO << outTxt << endl;
        if (!fpO) throw runtime_error("Cannot write output '" + outputFileName_ + "'.");
    } else {
        output_stream() << outTxt << endl;
        if (!output_stream()) throw runtime_error("Cannot write standard output.");
    }

    return 0;
}
