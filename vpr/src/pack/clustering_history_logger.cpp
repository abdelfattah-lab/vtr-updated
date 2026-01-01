/**
 * @file
 * @brief   Implementation of CLB creation history logging.
 */

#include "clustering_history_logger.h"

#include <algorithm>
#include <functional>
#include <iomanip>
#include <queue>
#include <set>
#include <sstream>

#include "atom_netlist.h"
#include "echo_files.h"
#include "globals.h"
#include "lb_type_rr_graph.h"
#include "pack_types.h"
#include "physical_types.h"
#include "vpr_types.h"
#include "vtr_util.h"

// Global instance
ClusteringHistoryLogger* g_clustering_history_logger = nullptr;

ClusteringHistoryLogger::ClusteringHistoryLogger() {
    // Check if echo is enabled for clustering history
    if (isEchoFileEnabled(E_ECHO_CLUSTERING_HISTORY)) {
        const char* filename = getEchoFileName(E_ECHO_CLUSTERING_HISTORY);
        if (filename) {
            file_.open(filename);
            if (file_.is_open()) {
                file_ << "================================================================================\n";
                file_ << "                     CLUSTERING HISTORY LOG\n";
                file_ << "================================================================================\n";
                file_ << "\n";
                file_ << "This file logs the creation of CLBs during clustering, including:\n";
                file_ << "  - Route type: SKIP_INTRA_LB_ROUTE (quick) vs FULL (detailed routing)\n";
                file_ << "  - Success/Failure status of each molecule packing attempt\n";
                file_ << "  - Molecules and their placements within the CLB\n";
                file_ << "  - Congestion details on routing failures\n";
                file_ << "  - Timing for each CLB creation\n";
                file_ << "\n";
                file_ << "================================================================================\n\n";
            }
        }
    }

    // Check if echo is enabled for clustering profile (finalized CLB summary)
    if (isEchoFileEnabled(E_ECHO_CLUSTERING_PROFILE)) {
        const char* profile_filename = getEchoFileName(E_ECHO_CLUSTERING_PROFILE);
        if (profile_filename) {
            profile_file_.open(profile_filename);
        }
    }
}

ClusteringHistoryLogger::~ClusteringHistoryLogger() {
    if (file_.is_open()) {
        file_ << "\n================================================================================\n";
        file_ << "                     END OF CLUSTERING HISTORY LOG\n";
        file_ << "================================================================================\n";
        file_.close();
    }
    if (profile_file_.is_open()) {
        profile_file_.close();
    }
}

std::string ClusteringHistoryLogger::status_to_string(e_block_pack_status status) {
    switch (status) {
        case e_block_pack_status::BLK_PASSED:
            return "PASSED";
        case e_block_pack_status::BLK_FAILED_FEASIBLE:
            return "FAILED_FEASIBLE";
        case e_block_pack_status::BLK_FAILED_ROUTE:
            return "FAILED_ROUTE";
        case e_block_pack_status::BLK_FAILED_FLOORPLANNING:
            return "FAILED_FLOORPLANNING";
        case e_block_pack_status::BLK_FAILED_NOC_GROUP:
            return "FAILED_NOC_GROUP";
        case e_block_pack_status::BLK_STATUS_UNDEFINED:
        default:
            return "UNDEFINED";
    }
}

std::string ClusteringHistoryLogger::strategy_to_string(ClusterLegalizationStrategy strategy) {
    switch (strategy) {
        case ClusterLegalizationStrategy::FULL:
            return "FULL (with intra-LB routing)";
        case ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE:
            return "SKIP_INTRA_LB_ROUTE (quick)";
        default:
            return "UNKNOWN";
    }
}

std::string ClusteringHistoryLogger::get_atom_name(const t_pack_molecule* molecule, int index) const {
    if (!molecule || index < 0 || index >= static_cast<int>(molecule->atom_block_ids.size())) {
        return "<invalid>";
    }
    AtomBlockId atom_id = molecule->atom_block_ids[index];
    if (!atom_id.is_valid()) {
        return "<empty>";
    }
    try {
        const auto& atom_ctx = g_vpr_ctx.atom();
        return atom_ctx.nlist.block_name(atom_id);
    } catch (...) {
        return "<error>";
    }
}

std::string ClusteringHistoryLogger::get_placement_description(t_pb_graph_node* primitive) const {
    if (!primitive) {
        return "<no placement>";
    }

    std::string desc;
    const t_pb_graph_node* curr = primitive;

    // Build path from root to primitive
    std::vector<std::string> path;
    while (curr != nullptr) {
        if (curr->pb_type && curr->pb_type->name) {
            std::string node_name = curr->pb_type->name;
            node_name += "[" + std::to_string(curr->placement_index) + "]";
            path.push_back(node_name);
        } else {
            path.push_back("<unknown>[" + std::to_string(curr->placement_index) + "]");
        }
        curr = curr->parent_pb_graph_node;
    }

    // Reverse to get root-to-leaf order
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        if (!desc.empty()) {
            desc += "/";
        }
        desc += *it;
    }

    return desc;
}

std::string ClusteringHistoryLogger::get_placement_description_with_mode(const t_pb* atom_pb) const {
    if (!atom_pb || !atom_pb->pb_graph_node) {
        return "<no placement>";
    }

    std::string desc;

    // Build path from root to primitive, including mode information from t_pb
    std::vector<std::string> path;
    const t_pb* curr_pb = atom_pb;

    while (curr_pb != nullptr) {
        const t_pb_graph_node* gnode = curr_pb->pb_graph_node;
        if (gnode && gnode->pb_type && gnode->pb_type->name) {
            std::string node_name = gnode->pb_type->name;
            node_name += "[" + std::to_string(gnode->placement_index) + "]";

            // Add mode name if this node has multiple modes and a parent selected a mode
            // The mode is stored on the pb and indicates which child mode is selected
            if (gnode->pb_type->num_modes > 1 && curr_pb->mode >= 0 &&
                curr_pb->mode < gnode->pb_type->num_modes) {
                const char* mode_name = gnode->pb_type->modes[curr_pb->mode].name;
                if (mode_name) {
                    node_name += "[" + std::string(mode_name) + "]";
                }
            }

            path.push_back(node_name);
        } else if (gnode) {
            path.push_back("<unknown>[" + std::to_string(gnode->placement_index) + "]");
        }
        curr_pb = curr_pb->parent_pb;
    }

    // Reverse to get root-to-leaf order
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        if (!desc.empty()) {
            desc += "/";
        }
        desc += *it;
    }

    return desc;
}

std::string ClusteringHistoryLogger::describe_congestion(const t_lb_router_data* router_data) const {
    if (!router_data || !router_data->lb_type_graph || !router_data->lb_rr_node_stats) {
        return "  <No routing data available>\n";
    }

    std::ostringstream oss;
    const auto& lb_type_graph = *router_data->lb_type_graph;
    const auto* lb_rr_node_stats = router_data->lb_rr_node_stats;
    const auto& atom_ctx = g_vpr_ctx.atom();
    t_logical_block_type_ptr lb_type = router_data->lb_type;

    // Find congested nodes
    std::vector<int> congested_nodes;
    for (size_t inode = 0; inode < lb_type_graph.size(); ++inode) {
        if (lb_rr_node_stats[inode].occ > lb_type_graph[inode].capacity) {
            congested_nodes.push_back(inode);
        }
    }

    if (congested_nodes.empty()) {
        oss << "  <No congested routing nodes found>\n";
        return oss.str();
    }

    // Build a map from congested nodes to nets using them
    std::multimap<int, AtomNetId> congested_node_to_nets;
    if (router_data->intra_lb_nets) {
        const auto& lb_nets = *router_data->intra_lb_nets;
        for (size_t inet = 0; inet < lb_nets.size(); inet++) {
            const auto& lb_net = lb_nets[inet];
            if (!lb_net.rt_tree) continue;

            // BFS through the route tree with safety limit
            std::queue<const t_lb_trace*> q;
            q.push(lb_net.rt_tree);
            size_t max_iterations = lb_type_graph.size() * 10; // Safety limit
            size_t iterations = 0;
            while (!q.empty() && iterations < max_iterations) {
                iterations++;
                const t_lb_trace* curr = q.front();
                q.pop();
                if (!curr) continue;

                int inode = curr->current_node;
                if (inode >= 0 && static_cast<size_t>(inode) < lb_type_graph.size() &&
                    lb_rr_node_stats[inode].occ > lb_type_graph[inode].capacity) {
                    congested_node_to_nets.insert({inode, lb_net.atom_net_id});
                }

                for (const auto& next : curr->next_nodes) {
                    q.push(&next);
                }
            }
        }
    }

    oss << "  CONGESTED INPUT PINS:\n";
    for (int inode : congested_nodes) {
        const t_lb_type_rr_node& rr_node = lb_type_graph[inode];
        const t_lb_rr_node_stats& rr_stats = lb_rr_node_stats[inode];

        // Describe the node
        std::string node_desc;
        if (rr_node.pb_graph_pin) {
            node_desc = rr_node.pb_graph_pin->to_string(false);
        } else if (lb_type && inode == get_lb_type_rr_graph_ext_source_index(lb_type)) {
            node_desc = "cluster-external source";
        } else if (lb_type && inode == get_lb_type_rr_graph_ext_sink_index(lb_type)) {
            node_desc = "cluster-external sink";
        } else {
            switch (rr_node.type) {
                case LB_SOURCE:
                    node_desc = "internal source";
                    break;
                case LB_SINK:
                    node_desc = "internal sink";
                    break;
                case LB_INTERMEDIATE:
                    node_desc = "intermediate node";
                    break;
                default:
                    node_desc = "unknown node";
            }
        }

        oss << "    RR Node " << inode << " (" << node_desc << ")\n";
        oss << "      Occupancy: " << rr_stats.occ << " > Capacity: " << rr_node.capacity << "\n";

        // List the nets using this node
        auto range = congested_node_to_nets.equal_range(inode);
        if (range.first != range.second) {
            oss << "      Competing nets:\n";
            for (auto itr = range.first; itr != range.second; ++itr) {
                AtomNetId net_id = itr->second;
                if (net_id.is_valid()) {
                    oss << "        - " << atom_ctx.nlist.net_name(net_id) << "\n";
                }
            }
        }
    }

    return oss.str();
}

void ClusteringHistoryLogger::log_clb_start(LegalizationClusterId cluster_id,
                                             const std::string& cluster_type_name,
                                             const t_pack_molecule* seed_molecule) {
    // Reset candidate failure stats for the new CLB
    reset_candidate_stats();

    if (!file_.is_open()) return;

    clb_count_++;
    clb_start_time_ = std::chrono::high_resolution_clock::now();

    file_ << "--------------------------------------------------------------------------------\n";
    file_ << "CLB #" << clb_count_ << " (ID: " << (size_t)cluster_id << ")\n";
    file_ << "--------------------------------------------------------------------------------\n";
    file_ << "  Type: " << cluster_type_name << "\n";

    if (seed_molecule) {
        file_ << "  Seed molecule: " << get_atom_name(seed_molecule, seed_molecule->root);
        if (seed_molecule->pack_pattern && seed_molecule->pack_pattern->name) {
            file_ << " (pattern: " << seed_molecule->pack_pattern->name << ")";
        }
        file_ << "\n";
        file_ << "  Molecule size: " << seed_molecule->num_blocks << " atom(s)\n";
    }
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_molecule_attempt(const t_pack_molecule* molecule,
                                                    ClusterLegalizationStrategy strategy,
                                                    e_block_pack_status status,
                                                    t_pb_graph_node** primitives_list,
                                                    int molecule_size,
                                                    double elapsed_us) {
    if (!file_.is_open() || !molecule) return;

    file_ << "  MOLECULE ATTEMPT:\n";
    file_ << "    Root atom: " << get_atom_name(molecule, molecule->root) << "\n";
    file_ << "    Route type: " << strategy_to_string(strategy) << "\n";
    file_ << "    Status: " << status_to_string(status) << "\n";
    file_ << "    Time: " << std::fixed << std::setprecision(3) << (elapsed_us / 1000.0) << " ms\n";

    // Only log placements if we have valid data and the status indicates we got far enough
    if (primitives_list && molecule_size > 0) {
        file_ << "    Atom placements:\n";
        int num_atoms = static_cast<int>(molecule->atom_block_ids.size());
        int limit = std::min(molecule_size, num_atoms);
        for (int i = 0; i < limit; i++) {
            if (i < 0 || i >= num_atoms) continue;
            if (!molecule->atom_block_ids[i].is_valid()) continue;

            file_ << "      [" << i << "] " << get_atom_name(molecule, i) << "\n";
            if (primitives_list[i]) {
                file_ << "          -> " << get_placement_description(primitives_list[i]) << "\n";
            } else {
                file_ << "          -> <no placement>\n";
            }
        }
    }
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_routing_failure(const t_pack_molecule* molecule,
                                                   const t_lb_router_data* router_data,
                                                   t_pb_graph_node** primitives_list,
                                                   int molecule_size) {
    if (!file_.is_open()) return;

    file_ << "  ROUTING FAILURE DETAILS:\n";

    if (molecule) {
        file_ << "    Failed molecule: " << get_atom_name(molecule, molecule->root) << "\n";

        // Log attempted placements if available
        if (primitives_list && molecule_size > 0) {
            file_ << "    Attempted placements:\n";
            int num_atoms = static_cast<int>(molecule->atom_block_ids.size());
            int limit = std::min(molecule_size, num_atoms);
            for (int i = 0; i < limit; i++) {
                if (i < 0 || i >= num_atoms) continue;
                if (!molecule->atom_block_ids[i].is_valid()) continue;
                file_ << "      [" << i << "] " << get_atom_name(molecule, i);
                if (primitives_list[i]) {
                    file_ << " -> " << get_placement_description(primitives_list[i]);
                }
                file_ << "\n";
            }
        }
    }

    // Log congestion details
    file_ << "\n    CONGESTION ANALYSIS:\n";
    file_ << describe_congestion(router_data);

    // Note: Detailed routing paths are not available here because trees are freed
    // after routing failure. The failure_description in router_data contains
    // congestion info captured before cleanup.

    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_routing_paths(const t_lb_router_data* router_data) {
    if (!file_.is_open()) return;

    file_ << "\n  ======== DETAILED ROUTING PATHS ========\n";

    if (!router_data || !router_data->lb_type_graph || !router_data->intra_lb_nets) {
        file_ << "    <No routing data available>\n";
        return;
    }

    const auto& atom_ctx = g_vpr_ctx.atom();
    const auto& lb_type_graph = *router_data->lb_type_graph;
    const auto& lb_nets = *router_data->intra_lb_nets;
    const auto* lb_rr_node_stats = router_data->lb_rr_node_stats;
    t_logical_block_type_ptr lb_type = router_data->lb_type;

    // Find congested nodes for highlighting
    std::set<int> congested_nodes;
    if (lb_rr_node_stats) {
        for (size_t inode = 0; inode < lb_type_graph.size(); ++inode) {
            if (lb_rr_node_stats[inode].occ > lb_type_graph[inode].capacity) {
                congested_nodes.insert(inode);
            }
        }
    }

    // Helper to get pin name from RR node
    auto get_node_name = [&](int inode) -> std::string {
        if (inode < 0 || static_cast<size_t>(inode) >= lb_type_graph.size()) {
            return "<invalid>";
        }
        const t_lb_type_rr_node& rr_node = lb_type_graph[inode];
        if (rr_node.pb_graph_pin) {
            return rr_node.pb_graph_pin->to_string(false);
        } else if (lb_type && inode == get_lb_type_rr_graph_ext_source_index(lb_type)) {
            return "EXT_SOURCE";
        } else if (lb_type && inode == get_lb_type_rr_graph_ext_sink_index(lb_type)) {
            return "EXT_SINK";
        } else {
            switch (rr_node.type) {
                case LB_SOURCE: return "SOURCE";
                case LB_SINK: return "SINK";
                case LB_INTERMEDIATE: return "INTERMEDIATE";
                default: return "UNKNOWN";
            }
        }
    };

    // Helper to check if a node is congested
    auto is_congested = [&](int inode) -> bool {
        return congested_nodes.count(inode) > 0;
    };

    // Recursive function to trace and print route tree
    std::function<void(const t_lb_trace*, int, std::vector<std::string>&)> trace_route;
    trace_route = [&](const t_lb_trace* trace, int depth, std::vector<std::string>& path) {
        if (!trace) return;

        int inode = trace->current_node;
        std::string node_name = get_node_name(inode);

        // Mark congested nodes with ***
        if (is_congested(inode)) {
            node_name = "***" + node_name + "*** (CONGESTED)";
        }

        path.push_back(node_name);

        if (trace->next_nodes.empty()) {
            // This is a sink - print the complete path
            file_ << "      PATH: ";
            for (size_t i = 0; i < path.size(); ++i) {
                if (i > 0) file_ << " -> ";
                file_ << path[i];
            }
            file_ << "\n";
        } else {
            // Continue tracing to children
            for (const auto& next : trace->next_nodes) {
                trace_route(&next, depth + 1, path);
            }
        }

        path.pop_back();
    };

    file_ << "    (Paths marked with *** indicate congested nodes)\n\n";

    // Process each net
    size_t net_count = 0;
    for (size_t inet = 0; inet < lb_nets.size(); inet++) {
        const auto& lb_net = lb_nets[inet];
        if (!lb_net.rt_tree) continue;

        // Get net name
        std::string net_name;
        if (lb_net.atom_net_id.is_valid()) {
            net_name = atom_ctx.nlist.net_name(lb_net.atom_net_id);
        } else {
            net_name = "<unknown_net_" + std::to_string(inet) + ">";
        }

        // Get source and sink pin info
        std::string source_info = "<unknown_source>";
        std::vector<std::string> sink_infos;

        for (size_t term = 0; term < lb_net.terminals.size(); term++) {
            int term_rr_node = lb_net.terminals[term];
            std::string term_pin_name = get_node_name(term_rr_node);

            // Try to get atom pin info
            std::string atom_pin_info;
            if (term < lb_net.atom_pins.size() && lb_net.atom_pins[term].is_valid()) {
                AtomPinId pin_id = lb_net.atom_pins[term];
                AtomBlockId blk_id = atom_ctx.nlist.pin_block(pin_id);
                AtomPortId port_id = atom_ctx.nlist.pin_port(pin_id);
                int pin_index = atom_ctx.nlist.pin_port_bit(pin_id);
                std::string blk_name = atom_ctx.nlist.block_name(blk_id);
                std::string port_name = atom_ctx.nlist.port_name(port_id);
                atom_pin_info = blk_name + "." + port_name + "[" + std::to_string(pin_index) + "]";
            }

            if (term == 0) {
                // Source
                source_info = atom_pin_info.empty() ? term_pin_name : atom_pin_info + " @ " + term_pin_name;
            } else {
                // Sink
                std::string sink = atom_pin_info.empty() ? term_pin_name : atom_pin_info + " @ " + term_pin_name;
                sink_infos.push_back(sink);
            }
        }

        file_ << "    NET [" << net_count++ << "]: " << net_name << "\n";
        file_ << "      Source: " << source_info << "\n";
        for (size_t i = 0; i < sink_infos.size(); i++) {
            file_ << "      Sink " << i << ": " << sink_infos[i] << "\n";
        }

        // Trace and print routing paths
        std::vector<std::string> path;
        trace_route(lb_net.rt_tree, 0, path);

        file_ << "\n";
    }

    if (net_count == 0) {
        file_ << "      <No routed nets found>\n";
    }

    file_.flush();
}

void ClusteringHistoryLogger::log_clb_success(LegalizationClusterId cluster_id,
                                               const t_pb* cluster_pb,
                                               const std::vector<t_pack_molecule*>& molecules,
                                               ClusterLegalizationStrategy strategy) {
    if (!file_.is_open()) return;

    file_ << "  CLB CREATION: SUCCESS\n";
    file_ << "    Route type: " << strategy_to_string(strategy) << "\n";
    file_ << "    Total molecules packed: " << molecules.size() << "\n";

    if (cluster_pb && cluster_pb->name) {
        file_ << "    CLB name: " << cluster_pb->name << "\n";
    }

    // List all molecules in the CLB
    file_ << "    Packed molecules:\n";
    const auto& atom_ctx = g_vpr_ctx.atom();
    for (size_t i = 0; i < molecules.size(); i++) {
        const t_pack_molecule* mol = molecules[i];
        if (!mol) continue;

        file_ << "      [" << i << "] Root: " << get_atom_name(mol, mol->root);
        if (mol->pack_pattern && mol->pack_pattern->name) {
            file_ << " (pattern: " << mol->pack_pattern->name << ")";
        }
        file_ << "\n";

        // List atoms and their placements (with mode info)
        int num_atoms = static_cast<int>(mol->atom_block_ids.size());
        int limit = std::min(mol->num_blocks, num_atoms);
        for (int j = 0; j < limit; j++) {
            AtomBlockId atom_id = mol->atom_block_ids[j];
            if (!atom_id.is_valid()) continue;

            const t_pb* atom_pb = atom_ctx.lookup.atom_pb(atom_id);
            file_ << "          Atom: " << atom_ctx.nlist.block_name(atom_id);
            if (atom_pb) {
                file_ << " @ " << get_placement_description_with_mode(atom_pb);
            }
            file_ << "\n";
        }
    }
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_clb_failure(LegalizationClusterId cluster_id,
                                               const std::string& reason,
                                               ClusterLegalizationStrategy strategy,
                                               const t_pb* cluster_pb,
                                               const std::vector<t_pack_molecule*>* molecules) {
    if (!file_.is_open()) return;

    file_ << "  CLB CREATION: FAILED\n";
    file_ << "    Route type: " << strategy_to_string(strategy) << "\n";
    file_ << "    Reason: " << reason << "\n";
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_clb_timing() {
    if (!file_.is_open()) return;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - clb_start_time_);
    double elapsed_ms = duration.count() / 1000.0;

    file_ << "  Time elapsed: " << std::fixed << std::setprecision(3) << elapsed_ms << " ms\n";
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_feasibility_failure(const t_pack_molecule* molecule,
                                                       int num_placements_tried,
                                                       const std::string& last_failure_reason) {
    if (!file_.is_open()) return;

    const auto& atom_ctx = g_vpr_ctx.atom();

    file_ << "  FEASIBILITY FAILURE DETAILS:\n";
    file_ << "    Primitive placements tried: " << num_placements_tried << "\n";
    file_ << "    Last failure reason: " << last_failure_reason << "\n";

    if (molecule) {
        file_ << "    Molecule info:\n";
        file_ << "      Type: ";
        if (molecule->type == MOLECULE_SINGLE_ATOM) {
            file_ << "single_atom";
        } else if (molecule->pack_pattern && molecule->pack_pattern->name) {
            file_ << molecule->pack_pattern->name;
        } else {
            file_ << "forced_pack";
        }
        file_ << "\n";
        file_ << "      Num blocks: " << molecule->num_blocks << "\n";

        if (molecule->is_chain() && molecule->chain_info) {
            file_ << "      Chain info:\n";
            file_ << "        is_long_chain: " << molecule->chain_info->is_long_chain << "\n";
            file_ << "        chain_id: " << molecule->chain_info->chain_id << "\n";
            file_ << "        required_entry_chain_id: " << molecule->required_entry_chain_id << "\n";
        }

        // List all atoms in the molecule
        file_ << "      Atoms in molecule:\n";
        for (int i = 0; i < molecule->num_blocks; i++) {
            AtomBlockId atom_id = molecule->atom_block_ids[i];
            if (!atom_id.is_valid()) {
                file_ << "        [" << i << "] <empty>\n";
                continue;
            }
            const t_model* model = atom_ctx.nlist.block_model(atom_id);
            file_ << "        [" << i << "] " << atom_ctx.nlist.block_name(atom_id);
            if (model) {
                file_ << " (primitive: " << model->name << ")";
            }
            if (molecule->root == i) {
                file_ << " [ROOT]";
            }
            file_ << "\n";
        }
    }

    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_feasibility_failure(const t_pack_molecule* molecule,
                                                       int num_placements_tried,
                                                       const std::string& last_failure_reason,
                                                       const std::vector<PlacementAttemptInfo>& placement_attempts) {
    if (!file_.is_open()) return;

    const auto& atom_ctx = g_vpr_ctx.atom();

    file_ << "  FEASIBILITY FAILURE DETAILS:\n";
    file_ << "    Primitive placements tried: " << num_placements_tried << "\n";
    file_ << "    Last failure reason: " << last_failure_reason << "\n";

    // Log each placement attempt with its failure reason
    if (!placement_attempts.empty()) {
        file_ << "    Placement attempts:\n";
        for (size_t i = 0; i < placement_attempts.size(); ++i) {
            file_ << "      [" << (i + 1) << "] " << placement_attempts[i].primitive_path << "\n";
            if (!placement_attempts[i].failure_reason.empty()) {
                file_ << "          Failure: " << placement_attempts[i].failure_reason << "\n";
            }
        }
    }

    if (molecule) {
        file_ << "    Molecule info:\n";
        file_ << "      Type: ";
        if (molecule->type == MOLECULE_SINGLE_ATOM) {
            file_ << "single_atom";
        } else if (molecule->pack_pattern && molecule->pack_pattern->name) {
            file_ << molecule->pack_pattern->name;
        } else {
            file_ << "forced_pack";
        }
        file_ << "\n";
        file_ << "      Num blocks: " << molecule->num_blocks << "\n";

        if (molecule->is_chain() && molecule->chain_info) {
            file_ << "      Chain info:\n";
            file_ << "        is_long_chain: " << molecule->chain_info->is_long_chain << "\n";
            file_ << "        chain_id: " << molecule->chain_info->chain_id << "\n";
            file_ << "        required_entry_chain_id: " << molecule->required_entry_chain_id << "\n";
        }

        // List all atoms in the molecule
        file_ << "      Atoms in molecule:\n";
        for (int i = 0; i < molecule->num_blocks; i++) {
            AtomBlockId atom_id = molecule->atom_block_ids[i];
            if (!atom_id.is_valid()) {
                file_ << "        [" << i << "] <empty>\n";
                continue;
            }
            const t_model* model = atom_ctx.nlist.block_model(atom_id);
            file_ << "        [" << i << "] " << atom_ctx.nlist.block_name(atom_id);
            if (model) {
                file_ << " (primitive: " << model->name << ")";
            }
            if (molecule->root == i) {
                file_ << " [ROOT]";
            }
            file_ << "\n";
        }
    }

    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_iteration_start(int iteration, ClusterLegalizationStrategy strategy) {
    if (!file_.is_open()) return;

    file_ << "\n";
    file_ << "  >>> SKIP_INTRA_LB_ROUTE failed - RETRYING with FULL routing <<<\n";
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::record_candidate_failure(e_block_pack_status status,
                                                        const std::string& primitive_type) {
    // Always record, even if file is not open (stats might be used elsewhere)
    candidate_failure_stats_[status][primitive_type]++;
}

void ClusteringHistoryLogger::reset_candidate_stats() {
    candidate_failure_stats_.clear();
    routing_stats_ = RoutingStats();  // Reset routing stats too
}

void ClusteringHistoryLogger::log_candidate_failure_stats() {
    if (!file_.is_open()) return;
    if (candidate_failure_stats_.empty()) return;

    // Calculate total failures
    int total_failures = 0;
    for (const auto& reason_pair : candidate_failure_stats_) {
        for (const auto& prim_pair : reason_pair.second) {
            total_failures += prim_pair.second;
        }
    }

    file_ << "  CANDIDATE PACKING FAILURES:\n";
    file_ << "    Total failures: " << total_failures << "\n";

    // Print breakdown by reason
    for (const auto& reason_pair : candidate_failure_stats_) {
        e_block_pack_status status = reason_pair.first;
        const auto& prim_counts = reason_pair.second;

        // Calculate total for this reason
        int reason_total = 0;
        for (const auto& prim_pair : prim_counts) {
            reason_total += prim_pair.second;
        }

        file_ << "    " << status_to_string(status) << ": " << reason_total << "\n";

        // Print breakdown by primitive type
        for (const auto& prim_pair : prim_counts) {
            file_ << "      " << prim_pair.first << ": " << prim_pair.second << "\n";
        }
    }
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::record_routing_attempt(size_t lb_rr_graph_size,
                                                      size_t num_nets,
                                                      int num_iterations,
                                                      double elapsed_us,
                                                      bool success,
                                                      bool is_impossible) {
    // Set RR graph size on first attempt
    if (routing_stats_.lb_rr_graph_size == 0) {
        routing_stats_.lb_rr_graph_size = lb_rr_graph_size;
    }

    routing_stats_.total_attempts++;
    if (success) {
        routing_stats_.successful_attempts++;
    } else {
        routing_stats_.failed_attempts++;
        // Track failure by primitive type
        std::string type_key = current_molecule_type_.empty() ? "<unknown>" : current_molecule_type_;
        routing_stats_.failed_by_type[type_key]++;
        if (is_impossible) {
            routing_stats_.impossible_attempts++;
            routing_stats_.impossible_by_type[type_key]++;
        }
    }

    // Track source using the current routing source
    if (current_routing_source_ == RoutingSource::MOLECULE) {
        routing_stats_.attempts_from_molecule++;
    } else {
        routing_stats_.attempts_from_legality++;
    }

    routing_stats_.total_nets_routed += num_nets;
    routing_stats_.total_iterations += num_iterations;
    routing_stats_.max_iterations = std::max(routing_stats_.max_iterations, num_iterations);
    routing_stats_.total_time_us += elapsed_us;
    routing_stats_.max_time_us = std::max(routing_stats_.max_time_us, elapsed_us);
}

void ClusteringHistoryLogger::record_mode_retry() {
    routing_stats_.mode_retry_count++;
}

void ClusteringHistoryLogger::log_routing_stats() {
    if (!file_.is_open()) return;
    if (routing_stats_.total_attempts == 0) return;

    file_ << "  INTRA-LB ROUTING STATISTICS:\n";
    file_ << "    LB RR graph size: " << routing_stats_.lb_rr_graph_size << " nodes\n";
    file_ << "    Total routing attempts: " << routing_stats_.total_attempts << "\n";
    file_ << "      Successful: " << routing_stats_.successful_attempts << "\n";
    file_ << "      Failed: " << routing_stats_.failed_attempts << "\n";
    // Breakdown of failures by primitive type
    if (!routing_stats_.failed_by_type.empty()) {
        file_ << "        By molecule type:\n";
        for (const auto& pair : routing_stats_.failed_by_type) {
            file_ << "          " << pair.first << ": " << pair.second << "\n";
        }
    }
    if (routing_stats_.impossible_attempts > 0) {
        file_ << "        Impossible (immediate dead-end): " << routing_stats_.impossible_attempts << "\n";
        // Breakdown of impossible by primitive type
        if (!routing_stats_.impossible_by_type.empty()) {
            file_ << "          By molecule type:\n";
            for (const auto& pair : routing_stats_.impossible_by_type) {
                file_ << "            " << pair.first << ": " << pair.second << "\n";
            }
        }
    }

    // Source breakdown
    file_ << "    Attempt sources:\n";
    file_ << "      From molecule packing: " << routing_stats_.attempts_from_molecule << "\n";
    file_ << "      From legality check: " << routing_stats_.attempts_from_legality << "\n";
    if (routing_stats_.mode_retry_count > 0) {
        file_ << "      Mode conflict retries: " << routing_stats_.mode_retry_count << "\n";
    }

    file_ << "    Total nets routed: " << routing_stats_.total_nets_routed << "\n";

    if (routing_stats_.total_attempts > 0) {
        double avg_nets = static_cast<double>(routing_stats_.total_nets_routed) / routing_stats_.total_attempts;
        file_ << "    Avg nets per attempt: " << std::fixed << std::setprecision(1) << avg_nets << "\n";
    }

    file_ << "    Pathfinder iterations:\n";
    file_ << "      Total: " << routing_stats_.total_iterations << "\n";
    file_ << "      Max per attempt: " << routing_stats_.max_iterations << "\n";

    if (routing_stats_.total_attempts > 0) {
        double avg_iters = static_cast<double>(routing_stats_.total_iterations) / routing_stats_.total_attempts;
        file_ << "      Avg per attempt: " << std::fixed << std::setprecision(1) << avg_iters << "\n";
    }

    file_ << "    Routing time:\n";
    file_ << "      Total: " << std::fixed << std::setprecision(3) << (routing_stats_.total_time_us / 1000.0) << " ms\n";
    file_ << "      Max per attempt: " << std::fixed << std::setprecision(3) << (routing_stats_.max_time_us / 1000.0) << " ms\n";

    if (routing_stats_.total_attempts > 0) {
        double avg_time = routing_stats_.total_time_us / routing_stats_.total_attempts;
        file_ << "      Avg per attempt: " << std::fixed << std::setprecision(3) << (avg_time / 1000.0) << " ms\n";
    }

    file_ << "\n";
    file_.flush();
}

// Helper function to collect atoms placed within a BLE5
static void collect_ble5_atoms(const t_pb* pb,
                               ClusteringHistoryLogger::Ble5Utilization& ble5_util) {
    if (!pb || !pb->pb_graph_node) return;

    // If this is a primitive with a name, it's a placed atom
    if (pb->pb_graph_node->is_primitive() && pb->name) {
        ble5_util.atoms.push_back(pb->name);
        return;
    }

    // Recurse into children
    if (pb->child_pbs) {
        const t_pb_type* pb_type = pb->pb_graph_node->pb_type;
        const t_mode* mode = &pb_type->modes[pb->mode];
        for (int child_type = 0; child_type < mode->num_pb_type_children; child_type++) {
            int num_children = mode->pb_type_children[child_type].num_pb;
            for (int child_inst = 0; child_inst < num_children; child_inst++) {
                const t_pb* child = &pb->child_pbs[child_type][child_inst];
                if (child->name) {
                    collect_ble5_atoms(child, ble5_util);
                }
            }
        }
    }
}

// Helper function to collect input pins for a specific pb node
static void collect_pb_input_pins(const t_pb* pb, const t_pb* root_pb,
                                  ClusteringHistoryLogger::PbNodeInfo& node_info) {
    if (!pb || !pb->pb_graph_node || !root_pb) return;

    const t_pb_graph_node* gnode = pb->pb_graph_node;
    const auto& atom_ctx = g_vpr_ctx.atom();

    // For primitives, check which input pins are used via atom netlist
    if (gnode->is_primitive() && pb->name) {
        AtomBlockId atom_id = atom_ctx.lookup.pb_atom(pb);
        if (!atom_id.is_valid()) return;

        for (int port = 0; port < gnode->num_input_ports; port++) {
            const char* port_name = gnode->input_pins[port][0].port->name;
            AtomPortId atom_port = atom_ctx.nlist.find_atom_port(
                atom_id, gnode->input_pins[port][0].port->model_port);

            for (int pin = 0; pin < gnode->num_input_pins[port]; pin++) {
                std::string pin_name = std::string(port_name) + "[" + std::to_string(pin) + "]";
                bool pin_used = false;
                if (atom_port.is_valid()) {
                    AtomPinId atom_pin = atom_ctx.nlist.port_pin(atom_port, pin);
                    if (atom_pin.is_valid()) {
                        AtomNetId net = atom_ctx.nlist.pin_net(atom_pin);
                        pin_used = net.is_valid();
                    }
                }
                if (pin_used) {
                    node_info.input_pins[pin_name].push_back("used");
                } else {
                    node_info.input_pins[pin_name] = std::vector<std::string>();
                }
            }
        }
    } else {
        // For non-primitives, check pb_route for which input pins have routing
        for (int port = 0; port < gnode->num_input_ports; port++) {
            const char* port_name = gnode->input_pins[port][0].port->name;
            for (int pin = 0; pin < gnode->num_input_pins[port]; pin++) {
                const t_pb_graph_pin* gpin = &gnode->input_pins[port][pin];
                int pin_id = gpin->pin_count_in_cluster;
                std::string pin_name = std::string(port_name) + "[" + std::to_string(pin) + "]";

                bool pin_used = false;
                if (root_pb->pb_route.count(pin_id)) {
                    const auto& route = root_pb->pb_route.at(pin_id);
                    pin_used = route.atom_net_id.is_valid() || (route.driver_pb_pin_id != OPEN);
                }
                if (pin_used) {
                    node_info.input_pins[pin_name].push_back("routed");
                } else {
                    node_info.input_pins[pin_name] = std::vector<std::string>();
                }
            }
        }
    }
}

// Recursively collect the full pb hierarchy within a BLE5
static ClusteringHistoryLogger::PbNodeInfo collect_pb_hierarchy(const t_pb* pb, const t_pb* root_pb) {
    ClusteringHistoryLogger::PbNodeInfo node_info;
    if (!pb || !pb->pb_graph_node) return node_info;

    const t_pb_graph_node* gnode = pb->pb_graph_node;

    // Set basic info
    node_info.pb_type_name = gnode->pb_type->name;
    node_info.pb_index = gnode->placement_index;
    node_info.is_primitive = gnode->is_primitive();

    if (node_info.is_primitive) {
        node_info.atom_name = pb->name ? pb->name : "";
        node_info.mode = "";  // Primitives don't have modes
    } else {
        node_info.atom_name = "";
        if (pb->mode < gnode->pb_type->num_modes) {
            node_info.mode = gnode->pb_type->modes[pb->mode].name;
        } else {
            node_info.mode = "<unknown>";
        }
    }

    // Collect input pins for this level
    collect_pb_input_pins(pb, root_pb, node_info);

    // Recurse into children (skip for primitives)
    if (!node_info.is_primitive && pb->child_pbs) {
        const t_pb_type* pb_type = gnode->pb_type;
        const t_mode* mode = &pb_type->modes[pb->mode];
        for (int child_type = 0; child_type < mode->num_pb_type_children; child_type++) {
            int num_children = mode->pb_type_children[child_type].num_pb;
            for (int child_inst = 0; child_inst < num_children; child_inst++) {
                const t_pb* child = &pb->child_pbs[child_type][child_inst];
                if (child->name) {
                    node_info.children.push_back(collect_pb_hierarchy(child, root_pb));
                }
            }
        }
    }

    return node_info;
}

// Forward declaration for recursive tracing
static void trace_primitive_to_ble5_inputs(const t_pb* pb, const t_pb* root_pb,
                                           const std::set<int>& ble5_input_pin_ids,
                                           const std::map<int, std::string>& pin_id_to_name,
                                           ClusteringHistoryLogger::Ble5Utilization& ble5_util,
                                           std::ofstream* debug_file);

// Helper function to collect BLE5-level input pin usage
// Traces from used primitive pins backward through pb_route to find which BLE5 inputs are used
static void collect_ble5_input_pins(const t_pb* ble5_pb, const t_pb* root_pb,
                                    ClusteringHistoryLogger::Ble5Utilization& ble5_util,
                                    std::ofstream* debug_file = nullptr) {
    if (!ble5_pb || !ble5_pb->pb_graph_node || !root_pb) return;

    const t_pb_graph_node* ble5_gnode = ble5_pb->pb_graph_node;

    // Build a set of BLE5 input pin IDs and a map to pin names
    std::set<int> ble5_input_pin_ids;
    std::map<int, std::string> pin_id_to_name;

    for (int port = 0; port < ble5_gnode->num_input_ports; port++) {
        const char* port_name = ble5_gnode->input_pins[port][0].port->name;
        for (int pin = 0; pin < ble5_gnode->num_input_pins[port]; pin++) {
            const t_pb_graph_pin* gpin = &ble5_gnode->input_pins[port][pin];
            int pin_id = gpin->pin_count_in_cluster;
            ble5_input_pin_ids.insert(pin_id);
            std::string pin_name = std::string(port_name) + "[" + std::to_string(pin) + "]";
            pin_id_to_name[pin_id] = pin_name;
            // Initialize all with empty vector (unused)
            ble5_util.input_pins[pin_name] = std::vector<std::string>();
        }
    }

    // Trace from each used primitive pin backward to find which BLE5 inputs are used
    trace_primitive_to_ble5_inputs(ble5_pb, root_pb, ble5_input_pin_ids, pin_id_to_name, ble5_util, debug_file);
}

// Recursively find primitives and trace their used pins back to BLE5 inputs
// debug_file is optional - if provided, writes debug info
static void trace_primitive_to_ble5_inputs(const t_pb* pb, const t_pb* root_pb,
                                           const std::set<int>& ble5_input_pin_ids,
                                           const std::map<int, std::string>& pin_id_to_name,
                                           ClusteringHistoryLogger::Ble5Utilization& ble5_util,
                                           std::ofstream* debug_file = nullptr) {
    if (!pb || !pb->pb_graph_node) return;

    // If this is a primitive, check its used input pins and trace back
    if (pb->pb_graph_node->is_primitive() && pb->name) {
        const auto& atom_ctx = g_vpr_ctx.atom();
        AtomBlockId atom_id = atom_ctx.lookup.pb_atom(pb);
        if (!atom_id.is_valid()) return;

        const t_pb_graph_node* prim_gnode = pb->pb_graph_node;

        if (debug_file && debug_file->is_open()) {
            *debug_file << "      [DEBUG] Primitive: " << pb->name
                        << " (type: " << prim_gnode->pb_type->name << ")\n";
            *debug_file << "        pb_route size on root: " << root_pb->pb_route.size() << "\n";
            *debug_file << "        BLE5 input pin IDs: ";
            for (int id : ble5_input_pin_ids) *debug_file << id << " ";
            *debug_file << "\n";
        }

        // Check each input port
        for (int port = 0; port < prim_gnode->num_input_ports; port++) {
            AtomPortId atom_port = atom_ctx.nlist.find_atom_port(
                atom_id, prim_gnode->input_pins[port][0].port->model_port);

            for (int pin = 0; pin < prim_gnode->num_input_pins[port]; pin++) {
                bool pin_is_used = false;
                if (atom_port.is_valid()) {
                    AtomPinId atom_pin = atom_ctx.nlist.port_pin(atom_port, pin);
                    if (atom_pin.is_valid()) {
                        AtomNetId net = atom_ctx.nlist.pin_net(atom_pin);
                        pin_is_used = net.is_valid();
                    }
                }

                if (pin_is_used) {
                    // This primitive pin is used - trace back through pb_route to find BLE5 input
                    const t_pb_graph_pin* gpin = &prim_gnode->input_pins[port][pin];
                    int current_pin_id = gpin->pin_count_in_cluster;

                    if (debug_file && debug_file->is_open()) {
                        *debug_file << "        Used pin: " << gpin->port->name << "[" << pin << "]"
                                    << " pin_count_in_cluster=" << current_pin_id << "\n";
                        *debug_file << "          Tracing: ";
                    }

                    // Trace backward through driver chain
                    int max_iterations = 100;  // Safety limit
                    bool found_ble5_input = false;
                    for (int i = 0; i < max_iterations; i++) {
                        if (debug_file && debug_file->is_open()) {
                            *debug_file << current_pin_id;
                        }

                        // Check if current pin is a BLE5 input
                        if (ble5_input_pin_ids.count(current_pin_id)) {
                            // Record the atom name and which port/pin it drives
                            std::string atom_with_pin = std::string(pb->name) + "." +
                                std::string(gpin->port->name) + "[" + std::to_string(pin) + "]";
                            ble5_util.input_pins[pin_id_to_name.at(current_pin_id)].push_back(atom_with_pin);
                            found_ble5_input = true;
                            if (debug_file && debug_file->is_open()) {
                                *debug_file << " -> FOUND BLE5 input: " << pin_id_to_name.at(current_pin_id);
                            }
                            break;
                        }

                        // Look up driver in pb_route
                        if (!root_pb->pb_route.count(current_pin_id)) {
                            if (debug_file && debug_file->is_open()) {
                                *debug_file << " -> NOT IN pb_route";
                            }
                            break;  // No routing info for this pin
                        }

                        const auto& route = root_pb->pb_route.at(current_pin_id);
                        if (route.driver_pb_pin_id == OPEN) {
                            if (debug_file && debug_file->is_open()) {
                                *debug_file << " -> driver=OPEN (source)";
                            }
                            break;  // No driver (reached source)
                        }

                        if (debug_file && debug_file->is_open()) {
                            *debug_file << " -> ";
                        }
                        current_pin_id = route.driver_pb_pin_id;
                    }

                    if (debug_file && debug_file->is_open()) {
                        if (!found_ble5_input) {
                            *debug_file << " [NOT FOUND]";
                        }
                        *debug_file << "\n";
                    }
                }
            }
        }
        return;
    }

    // Recurse into children
    if (pb->child_pbs) {
        const t_pb_type* pb_type = pb->pb_graph_node->pb_type;
        const t_mode* mode = &pb_type->modes[pb->mode];
        for (int child_type = 0; child_type < mode->num_pb_type_children; child_type++) {
            int num_children = mode->pb_type_children[child_type].num_pb;
            for (int child_inst = 0; child_inst < num_children; child_inst++) {
                const t_pb* child = &pb->child_pbs[child_type][child_inst];
                if (child->name) {
                    trace_primitive_to_ble5_inputs(child, root_pb, ble5_input_pin_ids,
                                                   pin_id_to_name, ble5_util, debug_file);
                }
            }
        }
    }
}

// Helper function to recursively find BLE5 nodes and their pin usage
static void collect_ble5_pin_usage(const t_pb* pb, const t_pb* root_pb,
                                   std::map<int, ClusteringHistoryLogger::FleUtilization>& fle_map,
                                   int current_fle_idx, const std::string& current_fle_mode,
                                   std::ofstream* debug_file = nullptr) {
    if (!pb || !pb->pb_graph_node) return;

    const std::string pb_type_name = pb->pb_graph_node->pb_type->name;

    // Check if this is a BLE5 node
    if (pb_type_name.find("ble5") != std::string::npos || pb_type_name == "ble5") {
        int ble5_idx = pb->pb_graph_node->placement_index;

        // Get the mode of this BLE5
        std::string ble5_mode = "<unknown>";
        if (pb->mode < pb->pb_graph_node->pb_type->num_modes) {
            ble5_mode = pb->pb_graph_node->pb_type->modes[pb->mode].name;
        }

        if (debug_file && debug_file->is_open()) {
            *debug_file << "    [DEBUG] BLE5[" << ble5_idx << "] mode=" << ble5_mode << "\n";
        }

        // Create BLE5 utilization entry
        ClusteringHistoryLogger::Ble5Utilization ble5_util;
        ble5_util.ble5_index = ble5_idx;
        ble5_util.ble5_mode = ble5_mode;

        // Collect atoms placed within this BLE5
        collect_ble5_atoms(pb, ble5_util);

        // Collect BLE5-level input pin usage (not primitive pins)
        collect_ble5_input_pins(pb, root_pb, ble5_util, debug_file);

        // Collect full pb hierarchy within this BLE5
        // We start from the children of BLE5, not BLE5 itself (since BLE5 info is already captured above)
        if (pb->child_pbs) {
            const t_pb_type* pb_type = pb->pb_graph_node->pb_type;
            const t_mode* mode = &pb_type->modes[pb->mode];
            for (int child_type = 0; child_type < mode->num_pb_type_children; child_type++) {
                int num_children = mode->pb_type_children[child_type].num_pb;
                for (int child_inst = 0; child_inst < num_children; child_inst++) {
                    const t_pb* child = &pb->child_pbs[child_type][child_inst];
                    if (child->name) {
                        ble5_util.hierarchy.push_back(collect_pb_hierarchy(child, root_pb));
                    }
                }
            }
        }

        // Add to FLE map
        if (fle_map.find(current_fle_idx) != fle_map.end()) {
            fle_map[current_fle_idx].ble5_usage[ble5_idx] = ble5_util;
        }
        return;  // Don't recurse further into BLE5 children
    }

    // Check if this is a FLE node
    int fle_idx = current_fle_idx;
    std::string fle_mode = current_fle_mode;
    if (pb_type_name.find("fle") != std::string::npos || pb_type_name == "fle") {
        fle_idx = pb->pb_graph_node->placement_index;
        if (pb->mode < pb->pb_graph_node->pb_type->num_modes) {
            fle_mode = pb->pb_graph_node->pb_type->modes[pb->mode].name;
        }
        // Create FLE entry if not exists
        if (fle_map.find(fle_idx) == fle_map.end()) {
            ClusteringHistoryLogger::FleUtilization fle_util;
            fle_util.fle_index = fle_idx;
            fle_util.fle_mode = fle_mode;
            fle_map[fle_idx] = fle_util;
        }
    }

    // Recurse into children
    if (pb->child_pbs) {
        const t_pb_type* pb_type = pb->pb_graph_node->pb_type;
        const t_mode* mode = &pb_type->modes[pb->mode];
        for (int child_type = 0; child_type < mode->num_pb_type_children; child_type++) {
            int num_children = mode->pb_type_children[child_type].num_pb;
            for (int child_inst = 0; child_inst < num_children; child_inst++) {
                const t_pb* child = &pb->child_pbs[child_type][child_inst];
                if (child->name) {  // Only process if child is used
                    collect_ble5_pin_usage(child, root_pb, fle_map, fle_idx, fle_mode, debug_file);
                }
            }
        }
    }
}

void ClusteringHistoryLogger::record_finalized_clb(LegalizationClusterId cluster_id,
                                                    const t_pb* cluster_pb,
                                                    const std::vector<t_pack_molecule*>& molecules) {
    FinalizedClbInfo info;
    info.cluster_id = cluster_id;

    if (cluster_pb) {
        info.cluster_name = cluster_pb->name ? cluster_pb->name : "<unnamed>";
        if (cluster_pb->pb_graph_node && cluster_pb->pb_graph_node->pb_type) {
            info.cluster_type = cluster_pb->pb_graph_node->pb_type->name;
        } else {
            info.cluster_type = "<unknown>";
        }
    }

    const auto& atom_ctx = g_vpr_ctx.atom();

    // Record molecules with their atoms grouped (similar to log_clb_success format)
    for (const auto* mol : molecules) {
        if (!mol) continue;

        MoleculeInfo mol_info;
        mol_info.root_atom_name = get_atom_name(mol, mol->root);
        mol_info.num_blocks = mol->num_blocks;
        if (mol->pack_pattern && mol->pack_pattern->name) {
            mol_info.pattern_name = mol->pack_pattern->name;
        }

        // Record atom placements within this molecule
        int num_atoms = static_cast<int>(mol->atom_block_ids.size());
        int limit = std::min(mol->num_blocks, num_atoms);
        for (int i = 0; i < limit; i++) {
            AtomBlockId atom_id = mol->atom_block_ids[i];
            if (!atom_id.is_valid()) continue;

            std::string atom_name = atom_ctx.nlist.block_name(atom_id);
            const t_pb* atom_pb = atom_ctx.lookup.atom_pb(atom_id);
            std::string placement = get_placement_description_with_mode(atom_pb);
            mol_info.atom_placements.emplace_back(atom_name, placement);
        }

        info.molecules.push_back(std::move(mol_info));
    }

    // Traverse pb hierarchy to collect FLE/BLE5 utilization with pin usage
    if (cluster_pb) {
        collect_ble5_pin_usage(cluster_pb, cluster_pb, info.fle_utilization, -1, "", nullptr);
    }

    // Count total and used FLEs
    info.total_fles = 10;  // Could be extracted from architecture
    info.used_fles = static_cast<int>(info.fle_utilization.size());

    finalized_clbs_.push_back(std::move(info));
}

// Helper function to recursively print pb hierarchy
static void print_pb_hierarchy(std::ofstream& out,
                               const std::vector<ClusteringHistoryLogger::PbNodeInfo>& nodes,
                               int indent_level) {
    std::string indent(indent_level * 2, ' ');

    for (const auto& node : nodes) {
        // Print node header
        if (node.is_primitive) {
            out << indent << node.pb_type_name << "[" << node.pb_index << "] (primitive)";
            if (!node.atom_name.empty()) {
                out << " atom=" << node.atom_name;
            }
            out << "\n";
        } else {
            out << indent << node.pb_type_name << "[" << node.pb_index << "]";
            if (!node.mode.empty()) {
                out << " mode=" << node.mode;
            }
            out << "\n";
        }

        // Print input pins for this node
        if (!node.input_pins.empty()) {
            out << indent << "  Pins:\n";
            for (const auto& [pin_name, targets] : node.input_pins) {
                if (targets.empty()) {
                    out << indent << "    " << pin_name << " = 0\n";
                } else {
                    out << indent << "    " << pin_name << " = 1\n";
                }
            }
        }

        // Recurse into children
        if (!node.children.empty()) {
            print_pb_hierarchy(out, node.children, indent_level + 1);
        }
    }
}

// Convert FleActiveMode to string for output
static std::string mode_to_string(FleActiveMode mode) {
    switch (mode) {
        case FleActiveMode::LUT5:
            return "LUT5";
        case FleActiveMode::SIMPLE_CHAIN:
            return "simple_chain";
        case FleActiveMode::CHAIN:
            return "chain";
        default:
            return "unknown";
    }
}

// Convert a multiset of active modes to a string like "{LUT5, LUT5, chain}"
static std::string mode_set_to_string(const std::multiset<FleActiveMode>& modes) {
    if (modes.empty()) {
        return "{}";
    }

    std::string result = "{";
    bool first = true;
    // Output in priority order: LUT5, chain, simple_chain (with duplicates)
    for (FleActiveMode m : {FleActiveMode::LUT5, FleActiveMode::CHAIN, FleActiveMode::SIMPLE_CHAIN}) {
        size_t count = modes.count(m);
        for (size_t i = 0; i < count; i++) {
            if (!first) result += ", ";
            result += mode_to_string(m);
            first = false;
        }
    }
    result += "}";
    return result;
}

// Collect all active modes for a BLE5 into a multiset (preserves duplicates across BLE5s)
static void collect_ble5_modes(const ClusteringHistoryLogger::Ble5Utilization& ble5,
                                const std::map<std::string, std::string>& atom_to_pattern,
                                std::multiset<FleActiveMode>& modes) {
    if (ble5.atoms.empty()) {
        return;
    }

    // Check BLE5 mode for LUT5
    std::string mode_lower = ble5.ble5_mode;
    std::transform(mode_lower.begin(), mode_lower.end(), mode_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (mode_lower.find("lut") != std::string::npos ||
        mode_lower.find("blut") != std::string::npos) {
        modes.insert(FleActiveMode::LUT5);
    }

    // Check atom patterns from molecules - use a local set to avoid duplicates within same BLE5
    bool has_simple_chain = false;
    bool has_chain = false;
    for (const auto& atom : ble5.atoms) {
        auto it = atom_to_pattern.find(atom);
        if (it != atom_to_pattern.end()) {
            const std::string& pattern = it->second;
            if (pattern == "simple_chain") {
                has_simple_chain = true;
            } else if (pattern.find("chain") != std::string::npos) {
                has_chain = true;
            }
        }
    }
    if (has_simple_chain) {
        modes.insert(FleActiveMode::SIMPLE_CHAIN);
    }
    if (has_chain) {
        modes.insert(FleActiveMode::CHAIN);
    }
}

void ClusteringHistoryLogger::write_summary() {
    if (!profile_file_.is_open()) return;

    profile_file_ << "================================================================================\n";
    profile_file_ << "                     FINALIZED CLB SUMMARY\n";
    profile_file_ << "================================================================================\n\n";

    profile_file_ << "Total CLBs created: " << finalized_clbs_.size() << "\n\n";

    for (const auto& clb : finalized_clbs_) {
        profile_file_ << "--------------------------------------------------------------------------------\n";
        profile_file_ << "CLB ID: " << size_t(clb.cluster_id) << " | Name: " << clb.cluster_name
              << " | Type: " << clb.cluster_type << "\n";
        profile_file_ << "--------------------------------------------------------------------------------\n";

        // FLE utilization
        double fle_util_percent = clb.total_fles > 0
            ? (100.0 * clb.used_fles / clb.total_fles)
            : 0.0;
        profile_file_ << "  FLE UTILIZATION: " << clb.used_fles << " / " << clb.total_fles
              << " (" << std::fixed << std::setprecision(1) << fle_util_percent << "%)\n\n";

        // FLE breakdown with BLE5 details and pin usage
        profile_file_ << "  FLE DETAILS:\n";
        for (const auto& [fle_idx, fle_util] : clb.fle_utilization) {
            profile_file_ << "    FLE[" << fle_idx << "] mode=" << fle_util.fle_mode << "\n";

            for (const auto& [ble5_idx, ble5_util] : fle_util.ble5_usage) {
                profile_file_ << "      BLE5[" << ble5_idx << "] mode=" << ble5_util.ble5_mode << "\n";

                // Show atoms in this BLE5 (one per line)
                if (!ble5_util.atoms.empty()) {
                    profile_file_ << "        Atoms:\n";
                    for (const auto& atom : ble5_util.atoms) {
                        profile_file_ << "          - " << atom << "\n";
                    }
                }

                // Show BLE5-level input pin usage
                if (!ble5_util.input_pins.empty()) {
                    profile_file_ << "        BLE5 Pins:\n";
                    for (const auto& [pin_name, targets] : ble5_util.input_pins) {
                        if (targets.empty()) {
                            profile_file_ << "          " << pin_name << " = 0\n";
                        } else {
                            profile_file_ << "          " << pin_name << " = 1 (";
                            for (size_t i = 0; i < targets.size(); i++) {
                                if (i > 0) profile_file_ << ", ";
                                profile_file_ << targets[i];
                            }
                            profile_file_ << ")\n";
                        }
                    }
                }

                // Show full hierarchy within BLE5
                if (!ble5_util.hierarchy.empty()) {
                    profile_file_ << "        Hierarchy:\n";
                    print_pb_hierarchy(profile_file_, ble5_util.hierarchy, 5);
                }
            }
        }

        // Molecules with grouped atom placements (similar to clustering_history format)
        profile_file_ << "\n  PACKED MOLECULES (" << clb.molecules.size() << "):\n";
        for (size_t i = 0; i < clb.molecules.size(); i++) {
            const auto& mol = clb.molecules[i];
            profile_file_ << "    [" << i << "] Root: " << mol.root_atom_name;
            if (!mol.pattern_name.empty()) {
                profile_file_ << " (pattern: " << mol.pattern_name << ")";
            }
            profile_file_ << "\n";

            // List atoms and their placements within this molecule
            for (const auto& [atom, placement] : mol.atom_placements) {
                profile_file_ << "        Atom: " << atom << " @ " << placement << "\n";
            }
        }

        profile_file_ << "\n";
    }

    // ==================== FLE ACTIVITY SUMMARY ====================
    // Build atom -> pattern map from all molecules across all CLBs
    std::map<std::string, std::string> atom_to_pattern;
    for (const auto& clb : finalized_clbs_) {
        for (const auto& mol : clb.molecules) {
            for (const auto& [atom_name, placement] : mol.atom_placements) {
                atom_to_pattern[atom_name] = mol.pattern_name;
            }
        }
    }

    // Count FLE mode combinations across all CLBs
    // Each FLE gets a multiset of active modes (preserves duplicates, e.g., {LUT5, LUT5})
    std::map<std::multiset<FleActiveMode>, int> fle_counts;
    for (const auto& clb : finalized_clbs_) {
        for (const auto& [fle_idx, fle_util] : clb.fle_utilization) {
            // Collect all active modes across all BLE5s in this FLE
            std::multiset<FleActiveMode> fle_modes;
            for (const auto& [ble5_idx, ble5_util] : fle_util.ble5_usage) {
                collect_ble5_modes(ble5_util, atom_to_pattern, fle_modes);
            }

            // Only count FLEs that have at least one active mode
            if (!fle_modes.empty()) {
                fle_counts[fle_modes]++;
            }
        }
    }

    // Output FLE activity summary
    profile_file_ << "================================================================================\n";
    profile_file_ << "                     FLE ACTIVITY SUMMARY\n";
    profile_file_ << "================================================================================\n\n";

    // Sort by count (descending) for readability
    std::vector<std::pair<std::multiset<FleActiveMode>, int>> sorted_counts(
        fle_counts.begin(), fle_counts.end());
    std::sort(sorted_counts.begin(), sorted_counts.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    int total_fles = 0;
    for (const auto& [modes, count] : sorted_counts) {
        profile_file_ << "  " << mode_set_to_string(modes) << ": " << count << "\n";
        total_fles += count;
    }
    profile_file_ << "\n  Total FLEs: " << total_fles << "\n";

    profile_file_.flush();
}
