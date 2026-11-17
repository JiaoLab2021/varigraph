#include <exception>
#include <iomanip>
#include <iostream>
#include <locale>
#include <string>

#include "include/CommandRunner.hpp"
#include "include/OptionParser.hpp"
#include "include/ProgramMetadata.hpp"
#include "include/get_time.hpp"
#include "include/logger.hpp"
#include "include/sys.hpp"

int main(int argc, char** argv) {
    try {
        const std::string command_line = program::cmdline(argc, argv);
        const auto parsed = parse_command(argc, argv);
        if (!parsed.run) return parsed.status;

        set_debug(parsed.config.common.debug);
        diagnostic_stream().imbue(std::locale(""));  // thousands separators in log numbers
        const double start_time = realtime();
        log_stream() << "Version: " << program::version << '\n';
        log_stream() << "Date: " << getTime() << '\n';
        log_stream() << "Command: " << command_line << '\n';

        print_config(parsed.command, parsed.config);
        run_command(parsed.command, parsed.config);

        diagnostic_stream().imbue(std::locale::classic());
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
