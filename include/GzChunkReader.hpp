#ifndef GZCHUNKREADER_HPP
#define GZCHUNKREADER_HPP

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string.h>
#include <zlib.h>
#include "logger.hpp"

// Keep zlib errors out of both line readers and kseq's unsigned buffer lengths.
inline int read_gzip_checked(gzFile file, void* buffer, unsigned int size) {
    const int bytes_read = gzread(file, buffer, size);
    int error_code = Z_OK;
    const char* error_message = gzerror(file, &error_code);
    if (bytes_read < 0 || (error_code != Z_OK && error_code != Z_STREAM_END)) {
        throw std::runtime_error(std::string("Cannot read compressed input: ") + error_message);
    }
    return bytes_read;
}

class GzChunkReader {
public:
    explicit GzChunkReader(const std::string& file_path, const size_t chunk_size = 1024 * 1024)  // Default 1Mb cache
        : file_path_(file_path), chunk_size_(chunk_size),
          chunk_buffer_(new char[chunk_size_]), buffer_pos_(0), buffer_end_(0)
    {
        // File handle is of type gzFile
        file_ = gzopen(file_path_.c_str(), "rb");

        if (!file_) {
            error_stream() << "'" << file_path_ << "': No such file or directory." << std::endl;
            exit(1);
        }
    }

    ~GzChunkReader() {
        if (file_) {
            gzclose(file_);
        }
    }

    bool read_line(std::string& line) {
        line.clear();

        while (true) {
            // Fill the buffer if it's empty
            if (buffer_pos_ == buffer_end_) {
                if (!read_chunk()) {
                    // If the end of the file is reached but there is no newline, return the last portion of content. Otherwise, one line will be missed.
                    if (!line.empty()) {
                        return true;
                    }
                    return false;
                }
            }

            // Extract a line from the buffer
            const char* line_begin = buffer_ + buffer_pos_;
            const char* line_end = static_cast<const char*>(memchr(line_begin, '\n', buffer_end_ - buffer_pos_));
            if (line_end) {
                line.append(line_begin, line_end - line_begin);
                buffer_pos_ += static_cast<size_t>(line_end - line_begin) + 1;
                return true;
            } else {
                line.append(line_begin, buffer_end_ - buffer_pos_);
                buffer_pos_ = buffer_end_;
            }
        }

        return false;
    }

private:
    bool read_chunk() {
        buffer_pos_ = 0;
        const int bytes_read = read_gzip_checked(file_, chunk_buffer_.get(), chunk_size_);
        buffer_end_ = static_cast<size_t>(bytes_read);
        if (buffer_end_ == 0) return false;
        buffer_ = chunk_buffer_.get();
        return true;
    }

private:
    std::string file_path_;
    size_t chunk_size_;
    std::unique_ptr<char[]> chunk_buffer_;
    char* buffer_;
    size_t buffer_pos_;
    size_t buffer_end_;
    gzFile file_;  // File handle is of type gzFile
};

#endif
