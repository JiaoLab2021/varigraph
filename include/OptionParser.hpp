#pragma once

#include "CommandOptions.hpp"

struct ParseResult {
    Command command = Command::Construct;
    AppConfig config;
    bool run = false;
    int status = 0;
};

ParseResult parse_command(int argc, char** argv, bool gpu_enabled);
