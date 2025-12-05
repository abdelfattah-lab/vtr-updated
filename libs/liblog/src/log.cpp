/**
 * Lightweight logging tool.  Automatically prepend messages with prefixes and store in log file.
 *
 * Author: Jason Luu
 * Date: Sept 5, 2014
 */

#include <stdio.h>
#include <stdarg.h> /* Allows for variable arguments, necessary for wrapping printf */
#include <string>
#include <vector>
#include <fstream>
#include "log.h"

#define LOG_DEFAULT_FILE_NAME "output.log"
#define MAX_LOG_SIZE_MB 100
#define LINES_TO_KEEP 10000

static int log_warning = 0;
static int log_error = 0;
FILE* log_stream = nullptr;
static std::string current_log_filename;
static int log_call_count = 0;

static void check_init();
static void check_and_truncate_log();

/* Set the output file of logger.
 * If different than current log file, close current log file and reopen to new log file
 * Modified to use append mode to prevent rewriting entire log files
 */
void log_set_output_file(const char* filename) {
    if (log_stream != nullptr) {
        fclose(log_stream);
    }

    if (filename == nullptr) {
        log_stream = nullptr;
        current_log_filename.clear();
    } else {
        current_log_filename = filename;
        log_stream = fopen(filename, "a");
        if (log_stream == nullptr) {
            printf("Error writing to file %s\n\n", filename);
        }
    }
}

void log_print_direct(const char* message, ...) {
    va_list args;
    va_start(args, message);
    vprintf(message, args);
    va_end(args);
}

void log_print_info(const char* message, ...) {
    check_init(); /* Check if output log file setup, if not, then this function also sets it up */
    check_and_truncate_log(); /* Check if log file is too large and truncate if needed */

    va_list args;
    va_start(args, message);
    vprintf(message, args);
    va_end(args);

    if (log_stream) {
        va_start(args, message); /* Must reset variable arguments so that they can be read again */
        vfprintf(log_stream, message, args);
        va_end(args);

        fflush(log_stream);
    }
}

void log_print_warning(const char* /*filename*/, unsigned int /*line_num*/, const char* message, ...) {
    check_init(); /* Check if output log file setup, if not, then this function also sets it up */
    check_and_truncate_log(); /* Check if log file is too large and truncate if needed */

    va_list args;
    va_start(args, message);
    log_warning++;

    printf("Warning %d: ", log_warning);
    vprintf(message, args);
    va_end(args);

    if (log_stream) {
        va_start(args, message); /* Must reset variable arguments so that they can be read again */
        fprintf(log_stream, "Warning %d: ", log_warning);
        vfprintf(log_stream, message, args);

        va_end(args);
        fflush(log_stream);
    }
}

void log_print_error(const char* /*filename*/, unsigned int /*line_num*/, const char* message, ...) {
    check_init(); /* Check if output log file setup, if not, then this function also sets it up */
    check_and_truncate_log(); /* Check if log file is too large and truncate if needed */

    va_list args;
    va_start(args, message);
    log_error++;

    fprintf(stderr, "Error %d: ", log_error);
    vfprintf(stderr, message, args);
    va_end(args);

    if (log_stream) {
        va_start(args, message); /* Must reset variable arguments so that they can be read again */
        fprintf(log_stream, "Error %d: ", log_error);
        vfprintf(log_stream, message, args);

        va_end(args);

        fflush(log_stream);
    }
}

/**
 * Check if output log file setup, if not, then this function also sets it up
 */
static void check_init() {
    //We now allow a nullptr log_stream (i.e. no log file) so nothing to do here
}

/**
 * Check log file size and truncate to last N lines if too large
 * Only checks every 1000 calls to avoid overhead
 */
static void check_and_truncate_log() {
    if (!log_stream || current_log_filename.empty()) {
        return;
    }

    // Only check every 1000 log calls to reduce overhead
    log_call_count++;
    if (log_call_count % 1000 != 0) {
        return;
    }

    // Check file size
    long current_pos = ftell(log_stream);
    fseek(log_stream, 0, SEEK_END);
    long file_size = ftell(log_stream);
    long max_size = MAX_LOG_SIZE_MB * 1024L * 1024L;

    if (file_size <= max_size) {
        // File is not too large, restore original position
        fseek(log_stream, current_pos, SEEK_SET);
        return;
    }

    // File is too large, truncate to last N lines
    fclose(log_stream);
    log_stream = nullptr;

    // Read all lines from the file
    std::ifstream infile(current_log_filename);
    std::vector<std::string> lines;
    std::string line;

    while (std::getline(infile, line)) {
        lines.push_back(line);
    }
    infile.close();

    // Keep only the last LINES_TO_KEEP lines
    size_t start_idx = 0;
    if (lines.size() > LINES_TO_KEEP) {
        start_idx = lines.size() - LINES_TO_KEEP;
    }

    // Rewrite the file with only the last N lines
    std::ofstream outfile(current_log_filename, std::ios::trunc);
    for (size_t i = start_idx; i < lines.size(); i++) {
        outfile << lines[i] << '\n';
    }
    outfile.close();

    // Reopen the file in append mode
    log_stream = fopen(current_log_filename.c_str(), "a");
    if (log_stream == nullptr) {
        printf("Error reopening log file %s after truncation\n", current_log_filename.c_str());
    }
}

void log_close() {
    if (log_stream) {
        fclose(log_stream);
    }
}
