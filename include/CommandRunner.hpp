#pragma once

#include "CommandOptions.hpp"

void print_config(Command command, const AppConfig& config, bool gpu_enabled);
int run_command(Command command, const AppConfig& config);
