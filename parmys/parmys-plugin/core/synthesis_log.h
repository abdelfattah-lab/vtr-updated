/*
 * Synthesis logging for debugging and visualization.
 * Traces operations from Verilog through to hardware structures.
 */

#ifndef SYNTHESIS_LOG_H
#define SYNTHESIS_LOG_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

// Forward declarations
struct nnode_t;
struct netlist_t;

namespace synthesis_log {

// Enable/disable logging globally
extern bool enabled;
extern std::ofstream log_file;

void init(const std::string& filename);
void close();

// Log an original operation from Verilog
void log_input_operation(nnode_t* node, const std::string& source_info);

// Log synthesis decision for a multiplier
void log_mult_decision(nnode_t* node, const std::string& route,
                       int input_a_width, int input_b_width,
                       const std::string& constant_value = "");

// Log synthesis decision for an adder
void log_adder_decision(nnode_t* node, const std::string& route,
                        int width, bool is_chained);

// Log a standalone $add cell (from Yosys, not from multiplier)
void log_standalone_adder(nnode_t* node, int input_a_width, int input_b_width,
                          const std::string& route);

// Log when an adder chain is created (sumout -> input connection)
void log_adder_chain(const std::string& source_mult_name,
                     const std::vector<std::string>& adder_names,
                     const std::vector<int>& adder_widths,
                     const std::string& chain_type);

// Log DP solution for constant multiplication
void log_dp_solution(const std::string& mult_name,
                     int num_rows,
                     int num_adders,
                     int total_bit_width,
                     const std::string& pairing_info);

// Log compressor tree structure
void log_compressor_tree(const std::string& mult_name,
                         const std::string& tree_type,
                         int input_height,
                         int output_rows,
                         int num_fa,
                         int num_ha);

// Summary statistics
void log_summary();

// Helper to get node info string
std::string get_node_info(nnode_t* node);

} // namespace synthesis_log

#endif // SYNTHESIS_LOG_H
