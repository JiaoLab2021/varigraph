#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "include/CommandRunner.hpp"
#include "include/OptionParser.hpp"
#include "include/ProgramMetadata.hpp"
#include "include/get_time.hpp"
#include "include/logger.hpp"
#include "include/sys.hpp"

int main(int argc, char** argv) {
#ifdef VARIGRAPH_GPU
    constexpr bool gpu_enabled = true;
#else
    constexpr bool gpu_enabled = false;
#endif

    try {
        const std::string command_line = program::cmdline(argc, argv);
        const auto parsed = parse_command(argc, argv, gpu_enabled);
        if (!parsed.run) return parsed.status;

        set_debug(parsed.config.common.debug);
        const double start_time = realtime();
        log_stream() << "Version: " << program::version << '\n';
        log_stream() << "Date: " << getTime() << '\n';
        log_stream() << "Command: " << command_line << '\n';

        print_config(parsed.command, parsed.config, gpu_enabled);

        const int status = run_command(parsed.command, parsed.config);
        if (status != 0) return status;

        stage_stream()
            << "Real time: " << std::fixed << std::setprecision(3)
            << (realtime() - start_time) << " sec; CPU: " << cputime()
            << " sec; Peak RSS: " << peakrss() / 1024.0 / 1024.0 / 1024.0 << " GB\n";
        return 0;
    } catch (const std::exception& error) {
        error_stream() << error.what() << '\n';
        return 1;
    }
}
