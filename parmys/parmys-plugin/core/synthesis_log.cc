/*
 * Synthesis logging implementation.
 */

#include "synthesis_log.h"
#include "odin_types.h"
#include <iomanip>
#include <ctime>

namespace synthesis_log {

bool enabled = false;
std::ofstream log_file;

// Counters for summary
static int total_mults = 0;
static int hard_mults = 0;
static int const_mults_ternary_dp = 0;
static int const_mults_binary_dp = 0;
static int const_mults_compressor = 0;
static int soft_mults = 0;
static int total_adders_created = 0;
static int chained_adders = 0;
static int standalone_adders = 0;

// Counters for standalone $add cells (from Yosys, not from multipliers)
static int total_standalone_adds = 0;
static int standalone_adds_to_hard = 0;
static int standalone_adds_soft = 0;

void init(const std::string& filename) {
    log_file.open(filename);
    if (log_file.is_open()) {
        enabled = true;

        // Header
        log_file << "================================================================================\n";
        log_file << "                    PARMYS SYNTHESIS LOG\n";
        log_file << "================================================================================\n";

        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        log_file << "Generated: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
        log_file << "================================================================================\n\n";

        // Reset counters
        total_mults = hard_mults = const_mults_ternary_dp = 0;
        const_mults_binary_dp = const_mults_compressor = soft_mults = 0;
        total_adders_created = chained_adders = standalone_adders = 0;
        total_standalone_adds = standalone_adds_to_hard = standalone_adds_soft = 0;
    }
}

void close() {
    if (enabled && log_file.is_open()) {
        log_summary();
        log_file.close();
        enabled = false;
    }
}

std::string get_node_info(nnode_t* node) {
    std::stringstream ss;
    if (node) {
        ss << "name=\"" << (node->name ? node->name : "unnamed") << "\"";
        if (node->loc.line >= 0) {
            ss << " line=" << node->loc.line;
        }
    }
    return ss.str();
}

void log_input_operation(nnode_t* node, const std::string& source_info) {
    if (!enabled) return;

    log_file << "[INPUT] ";
    if (node) {
        // Get operation type string
        const char* type_str = "UNKNOWN";
        switch (node->type) {
            case MULTIPLY: type_str = "MULTIPLY"; break;
            case ADD: type_str = "ADD"; break;
            case MINUS: type_str = "SUBTRACT"; break;
            default: break;
        }
        log_file << type_str << " " << get_node_info(node);

        if (node->num_input_port_sizes >= 2) {
            log_file << " ports=[" << node->input_port_sizes[0]
                     << "," << node->input_port_sizes[1] << "]";
        }
        if (node->num_output_port_sizes >= 1) {
            log_file << " out=" << node->output_port_sizes[0];
        }
    }
    if (!source_info.empty()) {
        log_file << " src=\"" << source_info << "\"";
    }
    log_file << "\n";
}

void log_mult_decision(nnode_t* node, const std::string& route,
                       int input_a_width, int input_b_width,
                       const std::string& constant_value) {
    if (!enabled) return;

    total_mults++;

    log_file << "\n[MULT] " << get_node_info(node) << "\n";
    log_file << "       Input widths: A=" << input_a_width << ", B=" << input_b_width << "\n";
    log_file << "       Route: " << route << "\n";

    if (!constant_value.empty()) {
        log_file << "       Constant: " << constant_value << "\n";
    }

    // Update counters
    if (route == "HARD_MULTIPLIER") {
        hard_mults++;
    } else if (route == "CONST_MULT_TERNARY_DP") {
        const_mults_ternary_dp++;
    } else if (route == "CONST_MULT_BINARY_DP") {
        const_mults_binary_dp++;
    } else if (route == "CONST_MULT_COMPRESSOR") {
        const_mults_compressor++;
    } else if (route == "SOFT_MULTIPLIER") {
        soft_mults++;
    }
}

void log_adder_decision(nnode_t* node, const std::string& route,
                        int width, bool is_chained) {
    if (!enabled) return;

    total_adders_created++;
    if (is_chained) {
        chained_adders++;
    } else {
        standalone_adders++;
    }

    log_file << "  [ADDER] " << get_node_info(node)
             << " width=" << width
             << " route=" << route
             << " chained=" << (is_chained ? "YES" : "NO") << "\n";
}

void log_standalone_adder(nnode_t* node, int input_a_width, int input_b_width,
                          const std::string& route) {
    if (!enabled) return;

    total_standalone_adds++;

    log_file << "\n[ADD] " << get_node_info(node) << "\n";
    log_file << "      Input widths: A=" << input_a_width << ", B=" << input_b_width << "\n";
    log_file << "      Route: " << route << "\n";

    // Update counters based on route
    // FROM_YOSYS means this is an original $add cell, counted but route not yet determined
    if (route == "HARD_ADDER" || route == "SPLIT_TO_HARD") {
        standalone_adds_to_hard++;
    } else if (route == "BELOW_THRESHOLD" || route == "SOFT") {
        standalone_adds_soft++;
    }
    // FROM_YOSYS: counted in total but route determined later by iterate_adders()
}

void log_adder_chain(const std::string& source_mult_name,
                     const std::vector<std::string>& adder_names,
                     const std::vector<int>& adder_widths,
                     const std::string& chain_type) {
    if (!enabled) return;

    log_file << "  [CHAIN] from \"" << source_mult_name << "\" type=" << chain_type << "\n";
    log_file << "          Structure (sumout -> input):\n";

    for (size_t i = 0; i < adder_names.size(); i++) {
        log_file << "            ";
        if (i > 0) {
            log_file << "|\n            v\n            ";
        }
        log_file << "[" << adder_names[i] << " w=" << adder_widths[i] << "]";
        if (i < adder_names.size() - 1) {
            log_file << " --sumout-->";
        }
        log_file << "\n";
    }
}

void log_dp_solution(const std::string& mult_name,
                     int num_rows,
                     int num_adders,
                     int total_bit_width,
                     const std::string& pairing_info) {
    if (!enabled) return;

    log_file << "  [DP] Solution for \"" << mult_name << "\":\n";
    log_file << "       Rows (partial products): " << num_rows << "\n";
    log_file << "       Adders created: " << num_adders << "\n";
    log_file << "       Total bit-width: " << total_bit_width << "\n";
    if (!pairing_info.empty()) {
        log_file << "       Pairings: " << pairing_info << "\n";
    }
}

void log_compressor_tree(const std::string& mult_name,
                         const std::string& tree_type,
                         int input_height,
                         int output_rows,
                         int num_fa,
                         int num_ha) {
    if (!enabled) return;

    log_file << "  [COMPRESSOR] for \"" << mult_name << "\":\n";
    log_file << "       Tree type: " << tree_type << "\n";
    log_file << "       Input height: " << input_height << " -> Output rows: " << output_rows << "\n";
    log_file << "       Full adders (FA): " << num_fa << ", Half adders (HA): " << num_ha << "\n";
}

void log_summary() {
    if (!enabled) return;

    log_file << "\n================================================================================\n";
    log_file << "                           SYNTHESIS SUMMARY\n";
    log_file << "================================================================================\n\n";

    log_file << "MULTIPLIERS:\n";
    log_file << "  Total:                    " << total_mults << "\n";
    log_file << "  Hard multipliers:         " << hard_mults << "\n";
    log_file << "  Const mult (ternary DP):  " << const_mults_ternary_dp << "\n";
    log_file << "  Const mult (binary DP):   " << const_mults_binary_dp << "\n";
    log_file << "  Const mult (compressor):  " << const_mults_compressor << "\n";
    log_file << "  Soft multipliers:         " << soft_mults << "\n\n";

    log_file << "STANDALONE ADDERS (from $add cells, NOT from multipliers):\n";
    log_file << "  Total:                    " << total_standalone_adds << "\n";
    log_file << "  Mapped to hard adders:    " << standalone_adds_to_hard << "\n";
    log_file << "  Soft/below threshold:     " << standalone_adds_soft << "\n";
    log_file << "  NOTE: These become SIMPLE carry chains (no ternary optimization)\n\n";

    log_file << "ADDERS CREATED (from all sources):\n";
    log_file << "  Total created:            " << total_adders_created << "\n";
    log_file << "  In chains (sumout->in):   " << chained_adders << "\n";
    log_file << "  Standalone:               " << standalone_adders << "\n\n";

    if (total_adders_created > 0) {
        float chain_ratio = 100.0f * chained_adders / total_adders_created;
        log_file << "  Chain ratio:              " << std::fixed << std::setprecision(1)
                 << chain_ratio << "%\n";
        log_file << "  (Higher is better for DCC3 packing)\n";
    }

    log_file << "\n================================================================================\n";
}

} // namespace synthesis_log
