/**
 * @file
 * @brief   Implementation of CLB creation history logging.
 */

#include "clustering_history_logger.h"

#include <algorithm>
#include <iomanip>
#include <queue>
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
}

ClusteringHistoryLogger::~ClusteringHistoryLogger() {
    if (file_.is_open()) {
        file_ << "\n================================================================================\n";
        file_ << "                     END OF CLUSTERING HISTORY LOG\n";
        file_ << "================================================================================\n";
        file_.close();
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
    file_ << "\n";
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

        // List atoms and their placements
        int num_atoms = static_cast<int>(mol->atom_block_ids.size());
        int limit = std::min(mol->num_blocks, num_atoms);
        for (int j = 0; j < limit; j++) {
            AtomBlockId atom_id = mol->atom_block_ids[j];
            if (!atom_id.is_valid()) continue;

            const t_pb* atom_pb = atom_ctx.lookup.atom_pb(atom_id);
            file_ << "          Atom: " << atom_ctx.nlist.block_name(atom_id);
            if (atom_pb && atom_pb->pb_graph_node) {
                file_ << " @ " << get_placement_description(atom_pb->pb_graph_node);
            }
            file_ << "\n";
        }
    }
    file_ << "\n";
    file_.flush();
}

void ClusteringHistoryLogger::log_clb_failure(LegalizationClusterId cluster_id,
                                               const std::string& reason,
                                               ClusterLegalizationStrategy strategy) {
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
