/**
 * @file
 * @author  Vaughn Betz (first revision - VPack),
 *          Alexander Marquardt (second revision - T-VPack),
 *          Jason Luu (third revision - AAPack),
 *          Alex Singer (fourth revision - APPack)
 * @date    June 8, 2011
 * @brief   Main clustering algorithm
 *
 * The clusterer uses several key data structures:
 *
 *      t_pb_type (and related types):
 *          Represent the architecture as described in the architecture file.
 *
 *      t_pb_graph_node (and related types):
 *          Represents a flattened version of the architecture with t_pb_types
 *          expanded (according to num_pb) into unique t_pb_graph_node instances,
 *          and the routing connectivity converted to a graph of t_pb_graph_pin (nodes)
 *          and t_pb_graph_edge.
 *
 *      t_pb:
 *          Represents a clustered instance of a t_pb_graph_node containing netlist primitives
 *
 *  t_pb_type and t_pb_graph_node (and related types) describe the targeted FPGA architecture, while t_pb represents
 *  the actual clustering of the user netlist.
 *
 *  For example:
 *      Consider an architecture where CLBs contain 4 BLEs, and each BLE is a LUT + FF pair.
 *      We wish to map a netlist of 400 LUTs and 400 FFs.
 *
 *      A BLE corresponds to one t_pb_type (which has num_pb = 4).
 *
 *      Each of the 4 BLE positions corresponds to a t_pb_graph_node (each of which references the BLE t_pb_type).
 *
 *      The output of clustering is 400 t_pb of type BLE which represent the clustered user netlist.
 *      Each of the 400 t_pb will reference one of the 4 BLE-type t_pb_graph_nodes.
 */

#include "greedy_clusterer.h"
#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include "SetupGrid.h"
#include "atom_netlist.h"
#include "attraction_groups.h"
#include "cluster_legalizer.h"
#include "cluster_util.h"
#include "clustering_history_logger.h"
#include "echo_files.h"
#include "constraints_report.h"
#include "greedy_seed_selector.h"
#include "pack_types.h"
#include "physical_types.h"
#include "prepack.h"
#include "vpr_context.h"
#include "vtr_math.h"
#include "vtr_vector.h"

// Set to true to enable verbose debug logging for multi-chain packing
constexpr bool MULTI_CHAIN_DEBUG = false;

namespace {

/**
 * @brief Struct to hold statistics on the progress of clustering.
 */
struct t_cluster_progress_stats {
    // The total number of molecules in the design.
    int num_molecules = 0;
    // The number of molecules which have been clustered.
    int num_molecules_processed = 0;
    // The number of molecules clustered since the last time the status was
    // logged.
    int mols_since_last_print = 0;
};

} // namespace

GreedyClusterer::GreedyClusterer(const t_packer_opts& packer_opts,
                                 const t_analysis_opts& analysis_opts,
                                 const AtomNetlist& atom_netlist,
                                 const t_arch& arch,
                                 const t_pack_high_fanout_thresholds& high_fanout_thresholds,
                                 const std::unordered_set<AtomNetId>& is_clock,
                                 const std::unordered_set<AtomNetId>& is_global)
    : packer_opts_(packer_opts)
    , analysis_opts_(analysis_opts)
    , atom_netlist_(atom_netlist)
    , arch_(arch)
    , high_fanout_thresholds_(high_fanout_thresholds)
    , is_clock_(is_clock)
    , is_global_(is_global)
    , primitive_candidate_block_types_(identify_primitive_candidate_block_types())
    , log_verbosity_(packer_opts.pack_verbosity)
    , net_output_feeds_driving_block_input_(identify_net_output_feeds_driving_block_input(atom_netlist)) {
}

std::map<t_logical_block_type_ptr, size_t>
GreedyClusterer::do_clustering(ClusterLegalizer& cluster_legalizer,
                               Prepacker& prepacker,
                               bool allow_unrelated_clustering,
                               bool balance_block_type_utilization,
                               AttractionInfo& attraction_groups,
                               DeviceContext& mutable_device_ctx) {
    // This routine returns a map that details the number of used block type
    // instances.
    std::map<t_logical_block_type_ptr, size_t> num_used_type_instances;

    /****************************************************************
     * Initialization
     *****************************************************************/

    // The clustering stats holds information used for logging the progress
    // of the clustering to the user.
    t_cluster_progress_stats clustering_stats;
    clustering_stats.num_molecules = prepacker.get_num_molecules();

    // TODO: Create a ClusteringTimingManager class.
    //       This code relies on the prepacker, once the prepacker is moved to
    //       the constructor, this code can also move to the constructor.
    std::shared_ptr<PreClusterDelayCalculator> clustering_delay_calc;
    std::shared_ptr<SetupTimingInfo> timing_info;
    // Default criticalities set to zero (e.g. if not timing driven)
    vtr::vector<AtomBlockId, float> atom_criticality(atom_netlist_.blocks().size(), 0.f);
    if (packer_opts_.timing_driven) {
        calc_init_packing_timing(packer_opts_, analysis_opts_, prepacker,
                                 clustering_delay_calc, timing_info, atom_criticality);
    }

    // Calculate the max molecule stats, which is used for gain calculation.
    const t_molecule_stats max_molecule_stats = prepacker.calc_max_molecule_stats(atom_netlist_);

    // Initialize the information for the greedy candidate selector.
    // TODO: Abstract into a candidate selector class.
    /* TODO: This is memory inefficient, fix if causes problems */
    /* Store stats on nets used by packed block, useful for determining transitively connected blocks
     * (eg. [A1, A2, ..]->[B1, B2, ..]->C implies cluster [A1, A2, ...] and C have a weak link) */
    vtr::vector<LegalizationClusterId, std::vector<AtomNetId>> clb_inter_blk_nets(atom_netlist_.blocks().size());
    // FIXME: This should be abstracted into a selector class. This is only used
    //        for gain calculation and selecting candidate molecules.
    t_clustering_data clustering_data;
    alloc_and_init_clustering(max_molecule_stats,
                              prepacker,
                              clustering_data,
                              clustering_stats.num_molecules);

    // Create the greedy seed selector.
    GreedySeedSelector seed_selector(atom_netlist_,
                                     prepacker,
                                     packer_opts_.cluster_seed_type,
                                     max_molecule_stats,
                                     atom_criticality);

    // Pick the first seed molecule.
    t_pack_molecule* seed_mol = seed_selector.get_next_seed(prepacker,
                                                            cluster_legalizer);

    /****************************************************************
     * Clustering
     *****************************************************************/

    print_pack_status_header();

    if (seed_mol) {
        AtomBlockId root_atom = seed_mol->atom_block_ids[seed_mol->root];
        if (root_atom.is_valid()) {
        }
    }

    // Continue clustering as long as a valid seed is returned from the seed
    // selector.
    while (seed_mol != nullptr) {
        // Check to ensure that this molecule is unclustered.
        VTR_ASSERT(!cluster_legalizer.is_mol_clustered(seed_mol));

        // The basic algorithm:
        // 1) Try to put all the molecules in that you can without doing the
        //    full intra-lb route. Then do full legalization at the end.
        // 2) If the legalization at the end fails, try again, but this time
        //    do full legalization for each molecule added to the cluster.

        // Check if we should use multi-chain packing for this seed
        bool use_multi_chain = packer_opts_.pack_multi_chain &&
                               seed_mol->is_chain() &&
                               seed_mol->chain_info &&
                               seed_mol->chain_info->is_long_chain;

        std::vector<LegalizationClusterId> new_cluster_ids;

        if (use_multi_chain) {
            // Use multi-chain packing to create multiple CLBs at once
            new_cluster_ids = try_grow_multi_chain_cluster(seed_mol,
                                                            ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE,
                                                            cluster_legalizer,
                                                            prepacker,
                                                            allow_unrelated_clustering,
                                                            balance_block_type_utilization,
                                                            *timing_info,
                                                            clb_inter_blk_nets,
                                                            clustering_data,
                                                            attraction_groups,
                                                            num_used_type_instances,
                                                            mutable_device_ctx);

            if (new_cluster_ids.empty()) {
                // If the previous strategy failed, try again with full legalization
                // for each molecule added to the cluster.

                // Log the retry with FULL strategy
                if (g_clustering_history_logger && g_clustering_history_logger->is_enabled()) {
                    g_clustering_history_logger->log_iteration_start(2, ClusterLegalizationStrategy::FULL);
                }

                new_cluster_ids = try_grow_multi_chain_cluster(seed_mol,
                                                                ClusterLegalizationStrategy::FULL,
                                                                cluster_legalizer,
                                                                prepacker,
                                                                allow_unrelated_clustering,
                                                                balance_block_type_utilization,
                                                                *timing_info,
                                                                clb_inter_blk_nets,
                                                                clustering_data,
                                                                attraction_groups,
                                                                num_used_type_instances,
                                                                mutable_device_ctx);
            }
        } else {
            // Try to grow a cluster from the seed molecule without doing intra-lb
            // route for each molecule (i.e. just use faster but not fully
            // conservative legality checks).
            LegalizationClusterId new_cluster_id = try_grow_cluster(seed_mol,
                                                                    ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE,
                                                                    cluster_legalizer,
                                                                    prepacker,
                                                                    allow_unrelated_clustering,
                                                                    balance_block_type_utilization,
                                                                    *timing_info,
                                                                    clb_inter_blk_nets,
                                                                    clustering_data,
                                                                    attraction_groups,
                                                                    num_used_type_instances,
                                                                    mutable_device_ctx);

            if (!new_cluster_id.is_valid()) {
                // If the previous strategy failed, try to grow the cluster again,
                // but this time perform full legalization for each molecule added
                // to the cluster.

                // Log the retry with FULL strategy
                if (g_clustering_history_logger && g_clustering_history_logger->is_enabled()) {
                    g_clustering_history_logger->log_iteration_start(2, ClusterLegalizationStrategy::FULL);
                }

                new_cluster_id = try_grow_cluster(seed_mol,
                                                  ClusterLegalizationStrategy::FULL,
                                                  cluster_legalizer,
                                                  prepacker,
                                                  allow_unrelated_clustering,
                                                  balance_block_type_utilization,
                                                  *timing_info,
                                                  clb_inter_blk_nets,
                                                  clustering_data,
                                                  attraction_groups,
                                                  num_used_type_instances,
                                                  mutable_device_ctx);
            }

            if (new_cluster_id.is_valid()) {
                new_cluster_ids.push_back(new_cluster_id);
            }
        }

        // Ensure that at the seed was packed successfully.
        VTR_ASSERT(!new_cluster_ids.empty());
        VTR_ASSERT(cluster_legalizer.is_mol_clustered(seed_mol));

        // Update the clustering progress stats for all created clusters.
        for (LegalizationClusterId new_cluster_id : new_cluster_ids) {
            size_t num_molecules_in_cluster = cluster_legalizer.get_num_molecules_in_cluster(new_cluster_id);
            clustering_stats.num_molecules_processed += num_molecules_in_cluster;
            clustering_stats.mols_since_last_print += num_molecules_in_cluster;
        }

        // Print the current progress of the packing after cluster(s) have been
        // successfully created.
        print_pack_status(clustering_stats.num_molecules,
                          clustering_stats.num_molecules_processed,
                          clustering_stats.mols_since_last_print,
                          mutable_device_ctx.grid.width(),
                          mutable_device_ctx.grid.height(),
                          attraction_groups,
                          cluster_legalizer);

        // Pick new seed.
        seed_mol = seed_selector.get_next_seed(prepacker,
                                               cluster_legalizer);
    }

    // If this architecture has LE physical block, report its usage.
    report_le_physical_block_usage(cluster_legalizer);

    // Write the summary of all finalized CLBs to clustering profile
    if (g_clustering_history_logger && g_clustering_history_logger->is_profile_enabled()) {
        g_clustering_history_logger->write_summary();
    }

    // Free the clustering data.
    // FIXME: This struct should use standard data structures so it does not
    //        have to be freed like this. This is also specific to the candidate
    //        gain calculation.
    free_clustering_data(clustering_data);

    return num_used_type_instances;
}

LegalizationClusterId GreedyClusterer::try_grow_cluster(
    t_pack_molecule* seed_mol,
    ClusterLegalizationStrategy strategy,
    ClusterLegalizer& cluster_legalizer,
    Prepacker& prepacker,
    bool allow_unrelated_clustering,
    bool balance_block_type_utilization,
    SetupTimingInfo& timing_info,
    vtr::vector<LegalizationClusterId, std::vector<AtomNetId>>& clb_inter_blk_nets,
    t_clustering_data& clustering_data,
    AttractionInfo& attraction_groups,
    std::map<t_logical_block_type_ptr, size_t>& num_used_type_instances,
    DeviceContext& mutable_device_ctx) {
    // Check to ensure that this molecule is unclustered.
    VTR_ASSERT(!cluster_legalizer.is_mol_clustered(seed_mol));

    // Set the legalization strategy of the cluster legalizer.
    cluster_legalizer.set_legalization_strategy(strategy);

    // Use the seed to start a new cluster.
    LegalizationClusterId legalization_cluster_id = start_new_cluster(seed_mol,
                                                                      cluster_legalizer,
                                                                      balance_block_type_utilization,
                                                                      num_used_type_instances,
                                                                      mutable_device_ctx);

    auto cluster_type = cluster_legalizer.get_cluster_type(legalization_cluster_id);
    // VTR_LOG("try_grow_cluster: cluster type name: %s\n", cluster_type->name.c_str());

    int high_fanout_threshold = high_fanout_thresholds_.get_threshold(cluster_type->name);

    update_cluster_stats(seed_mol,
                         cluster_legalizer,
                         is_clock_,  //Set of clock nets
                         is_global_, //Set of global nets (currently all clocks)
                         packer_opts_.global_clocks,
                         packer_opts_.alpha, packer_opts_.beta,
                         packer_opts_.timing_driven, packer_opts_.connection_driven,
                         high_fanout_threshold,
                         timing_info,
                         attraction_groups,
                         net_output_feeds_driving_block_input_);

    int num_unrelated_clustering_attempts = 0;
    t_pack_molecule* candidate_mol;
    candidate_mol = get_molecule_for_cluster(cluster_legalizer.get_cluster_pb(legalization_cluster_id),
                                             attraction_groups,
                                             allow_unrelated_clustering,
                                             packer_opts_.prioritize_transitive_connectivity,
                                             packer_opts_.transitive_fanout_threshold,
                                             packer_opts_.feasible_block_array_size,
                                             &num_unrelated_clustering_attempts,
                                             prepacker,
                                             cluster_legalizer,
                                             clb_inter_blk_nets,
                                             legalization_cluster_id,
                                             log_verbosity_,
                                             clustering_data.unclustered_list_head,
                                             clustering_data.unclustered_list_head_size,
                                             primitive_candidate_block_types_);

    /*
     * When attraction groups are created, the purpose is to pack more densely by adding more molecules
     * from the cluster's attraction group to the cluster. In a normal flow, (when attraction groups are
     * not on), the cluster keeps being packed until the get_molecule routines return either a repeated
     * molecule or a nullptr. When attraction groups are on, we want to keep exploring molecules for the
     * cluster until a nullptr is returned. So, the number of repeated molecules allowed is increased to a
     * large value.
     */
    int max_num_repeated_molecules = 1;
    if (attraction_groups.num_attraction_groups() > 0)
        max_num_repeated_molecules = attraction_groups_max_repeated_molecules_;

    // Continuously try to cluster candidate molecules into the cluster
    // until one of the following occurs:
    //  1) No candidate molecule is proposed.
    //  2) The same candidate was proposed multiple times.
    int num_repeated_molecules = 0;
    while (candidate_mol != nullptr && num_repeated_molecules < max_num_repeated_molecules) {
        // Try to cluster the candidate molecule into the cluster.
        bool success = try_add_candidate_mol_to_cluster(candidate_mol,
                                                        legalization_cluster_id,
                                                        cluster_legalizer);

        // If the candidate molecule was clustered successfully, update
        // the cluster stats.
        if (success) {
            update_cluster_stats(candidate_mol,
                                 cluster_legalizer,
                                 is_clock_,  //Set of all clocks
                                 is_global_, //Set of all global signals (currently clocks)
                                 packer_opts_.global_clocks,
                                 packer_opts_.alpha,
                                 packer_opts_.beta,
                                 packer_opts_.timing_driven,
                                 packer_opts_.connection_driven,
                                 high_fanout_threshold,
                                 timing_info,
                                 attraction_groups,
                                 net_output_feeds_driving_block_input_);
            num_unrelated_clustering_attempts = 0;
        }

        // Get the next candidate molecule.
        t_pack_molecule* prev_candidate_mol = candidate_mol;
        candidate_mol = get_molecule_for_cluster(cluster_legalizer.get_cluster_pb(legalization_cluster_id),
                                                 attraction_groups,
                                                 allow_unrelated_clustering,
                                                 packer_opts_.prioritize_transitive_connectivity,
                                                 packer_opts_.transitive_fanout_threshold,
                                                 packer_opts_.feasible_block_array_size,
                                                 &num_unrelated_clustering_attempts,
                                                 prepacker,
                                                 cluster_legalizer,
                                                 clb_inter_blk_nets,
                                                 legalization_cluster_id,
                                                 log_verbosity_,
                                                 clustering_data.unclustered_list_head,
                                                 clustering_data.unclustered_list_head_size,
                                                 primitive_candidate_block_types_);

        // If the next candidate molecule is the same as the previous
        // candidate molecule, increment the number of repreated
        // molecules counter.
        if (candidate_mol == prev_candidate_mol)
            num_repeated_molecules++;
    }

    // Ensure that the cluster is legal. When the cluster legalization
    // strategy is full, it must be legal.
    if (strategy != ClusterLegalizationStrategy::FULL) {
        // If the legalizer did not check everything for every molecule,
        // need to check that the full cluster is legal (need to perform
        // intra-lb routing).
        bool is_cluster_legal = cluster_legalizer.check_cluster_legality(legalization_cluster_id);

        if (!is_cluster_legal) {
            // Log CLB creation failure due to final legality check
            // Include the attempted molecules and placements for debugging
            if (g_clustering_history_logger && g_clustering_history_logger->is_enabled()) {
                g_clustering_history_logger->log_candidate_failure_stats();
                g_clustering_history_logger->log_routing_stats();
                // Get molecules and pb before the cluster is destroyed
                const auto& attempted_molecules = cluster_legalizer.get_cluster_molecules(legalization_cluster_id);
                const t_pb* attempted_pb = cluster_legalizer.get_cluster_pb(legalization_cluster_id);
                g_clustering_history_logger->log_clb_failure(
                    legalization_cluster_id,
                    "Final intra-LB routing check failed (congestion)",
                    strategy,
                    attempted_pb,
                    &attempted_molecules);
                g_clustering_history_logger->log_clb_timing();
            }

            // If the cluster is not legal, undo the cluster.
            // Update the used type instances.
            num_used_type_instances[cluster_legalizer.get_cluster_type(legalization_cluster_id)]--;
            // Destroy the illegal cluster.
            cluster_legalizer.destroy_cluster(legalization_cluster_id);
            cluster_legalizer.compress();
            // Cluster failed to grow.
            return LegalizationClusterId();
        }
    }

    VTR_ASSERT(legalization_cluster_id.is_valid());

    // Legal cluster was created. Store cluster info and clean cluster.

    // store info that will be used later in packing from pb_stats.
    // FIXME: If this is used for gain, it should be moved into the selector
    //        class. Perhaps a finalize_cluster_gain method.
    t_pb* cur_pb = cluster_legalizer.get_cluster_pb(legalization_cluster_id);
    t_pb_stats* pb_stats = cur_pb->pb_stats;
    for (const AtomNetId mnet_id : pb_stats->marked_nets) {
        int external_terminals = atom_netlist_.net_pins(mnet_id).size() - pb_stats->num_pins_of_net_in_pb[mnet_id];
        // Check if external terminals of net is within the fanout limit and
        // that there exists external terminals.
        if (external_terminals < packer_opts_.transitive_fanout_threshold && external_terminals > 0) {
            clb_inter_blk_nets[legalization_cluster_id].push_back(mnet_id);
        }
    }

    // Log CLB creation success to history file
    if (g_clustering_history_logger && g_clustering_history_logger->is_enabled()) {
        g_clustering_history_logger->log_candidate_failure_stats();
        g_clustering_history_logger->log_routing_stats();
        g_clustering_history_logger->log_clb_success(
            legalization_cluster_id,
            cluster_legalizer.get_cluster_pb(legalization_cluster_id),
            cluster_legalizer.get_cluster_molecules(legalization_cluster_id),
            strategy);
        g_clustering_history_logger->log_clb_timing();
    }

    // Since the cluster will no longer be added to beyond this point,
    // clean the cluster of any data not strictly necessary for
    // creating the clustered netlist.
    // NOTE: clean_cluster populates pb_route, so must be called before record_finalized_clb
    cluster_legalizer.clean_cluster(legalization_cluster_id);

    // Record finalized CLB for profile summary (separate file)
    // This must be called AFTER clean_cluster since that's when pb_route is populated
    if (g_clustering_history_logger && g_clustering_history_logger->is_profile_enabled()) {
        g_clustering_history_logger->record_finalized_clb(
            legalization_cluster_id,
            cluster_legalizer.get_cluster_pb(legalization_cluster_id),
            cluster_legalizer.get_cluster_molecules(legalization_cluster_id));
    }

    // Cluster has been grown successfully.
    return legalization_cluster_id;
}

std::vector<LegalizationClusterId> GreedyClusterer::try_grow_multi_chain_cluster(
    t_pack_molecule* seed_mol,
    ClusterLegalizationStrategy strategy,
    ClusterLegalizer& cluster_legalizer,
    Prepacker& prepacker,
    bool allow_unrelated_clustering,
    bool balance_block_type_utilization,
    SetupTimingInfo& timing_info,
    vtr::vector<LegalizationClusterId, std::vector<AtomNetId>>& clb_inter_blk_nets,
    t_clustering_data& clustering_data,
    AttractionInfo& attraction_groups,
    std::map<t_logical_block_type_ptr, size_t>& num_used_type_instances,
    DeviceContext& mutable_device_ctx) {

    std::vector<LegalizationClusterId> created_clusters;

    // Only use multi-chain packing for long chain seeds
    if (!seed_mol->is_chain() || !seed_mol->chain_info ||
        !seed_mol->chain_info->is_long_chain) {
        // Fall back to single cluster growth
        LegalizationClusterId cluster_id = try_grow_cluster(seed_mol,
                                                             strategy,
                                                             cluster_legalizer,
                                                             prepacker,
                                                             allow_unrelated_clustering,
                                                             balance_block_type_utilization,
                                                             timing_info,
                                                             clb_inter_blk_nets,
                                                             clustering_data,
                                                             attraction_groups,
                                                             num_used_type_instances,
                                                             mutable_device_ctx);
        if (cluster_id.is_valid()) {
            created_clusters.push_back(cluster_id);
        }
        return created_clusters;
    }

    // Get all molecules in the seed chain
    std::vector<t_pack_molecule*> seed_chain_mols = prepacker.get_chain_molecules(seed_mol);
    size_t chain_length = seed_chain_mols.size();

    if (chain_length <= 1) {
        // Single molecule chain, use regular clustering
        LegalizationClusterId cluster_id = try_grow_cluster(seed_mol,
                                                             strategy,
                                                             cluster_legalizer,
                                                             prepacker,
                                                             allow_unrelated_clustering,
                                                             balance_block_type_utilization,
                                                             timing_info,
                                                             clb_inter_blk_nets,
                                                             clustering_data,
                                                             attraction_groups,
                                                             num_used_type_instances,
                                                             mutable_device_ctx);
        if (cluster_id.is_valid()) {
            created_clusters.push_back(cluster_id);
        }
        return created_clusters;
    }

    if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
        VTR_LOG("Multi-chain packing: seed chain has %zu molecules\n", chain_length);
    }

    // Find suitable block type and mode for the seed
    AtomBlockId root_atom = seed_mol->atom_block_ids[seed_mol->root];
    const t_model* root_model = atom_netlist_.block_model(root_atom);
    auto itr = primitive_candidate_block_types_.find(root_model);
    VTR_ASSERT(itr != primitive_candidate_block_types_.end());
    std::vector<t_logical_block_type_ptr> candidate_types = itr->second;

    // Sort by utilization if requested
    if (balance_block_type_utilization) {
        std::stable_sort(candidate_types.begin(), candidate_types.end(),
                         [&](t_logical_block_type_ptr lhs, t_logical_block_type_ptr rhs) {
                             int lhs_num_instances = 0;
                             int rhs_num_instances = 0;
                             for (auto type : lhs->equivalent_tiles)
                                 lhs_num_instances += mutable_device_ctx.grid.num_instances(type, -1);
                             for (auto type : rhs->equivalent_tiles)
                                 rhs_num_instances += mutable_device_ctx.grid.num_instances(type, -1);
                             float lhs_util = vtr::safe_ratio<float>(num_used_type_instances[lhs], lhs_num_instances);
                             float rhs_util = vtr::safe_ratio<float>(num_used_type_instances[rhs], rhs_num_instances);
                             return lhs_util < rhs_util;
                         });
    }

    // Find a working type and mode
    // Temporarily suppress logging during test cluster creation
    ClusteringHistoryLogger* saved_logger = g_clustering_history_logger;
    g_clustering_history_logger = nullptr;

    t_logical_block_type_ptr chosen_type = nullptr;
    int chosen_mode = 0;
    for (auto type : candidate_types) {
        for (int mode = 0; mode < type->pb_graph_head->pb_type->num_modes; mode++) {
            // Quick check: try to start a test cluster with the seed
            auto [status, test_id] = cluster_legalizer.start_new_cluster(seed_mol, type, mode);
            if (status == e_block_pack_status::BLK_PASSED) {
                chosen_type = type;
                chosen_mode = mode;
                // Destroy the test cluster - we'll use batch creation instead
                cluster_legalizer.destroy_cluster(test_id);
                cluster_legalizer.compress();
                break;
            }
        }
        if (chosen_type != nullptr) break;
    }

    // Restore the logger
    g_clustering_history_logger = saved_logger;

    if (chosen_type == nullptr) {
        VPR_FATAL_ERROR(VPR_ERROR_PACK,
                        "Multi-chain packing: Cannot find suitable block type for seed molecule.\n");
    }

    // Set legalization strategy
    cluster_legalizer.set_legalization_strategy(strategy);

    // Create the batch of CLBs
    auto batch = cluster_legalizer.create_clb_batch(seed_mol, chosen_type, chosen_mode, chain_length);
    if (!batch || batch->clb_ids.size() != chain_length) {
        if (MULTI_CHAIN_DEBUG) {
            VTR_LOG("Multi-chain packing: Failed to create CLB batch\n");
        }
        return created_clusters;  // Empty
    }

    // Update used type instances count
    num_used_type_instances[chosen_type] += chain_length;

    // Expand FPGA if needed
    unsigned int num_instances = 0;
    for (auto equivalent_tile : chosen_type->equivalent_tiles) {
        num_instances += mutable_device_ctx.grid.num_instances(equivalent_tile, -1);
    }
    if (num_used_type_instances[chosen_type] > num_instances) {
        mutable_device_ctx.grid = create_device_grid(packer_opts_.device_layout,
                                                     arch_.grid_layouts,
                                                     num_used_type_instances,
                                                     packer_opts_.target_device_utilization);
    }

    // Pack seed chain molecules into the batch (one per CLB)
    // Get seed chain name for logging
    AtomBlockId seed_root_atom = seed_mol->atom_block_ids[seed_mol->root];
    std::string seed_chain_name = atom_netlist_.block_name(seed_root_atom);

    if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
        VTR_LOG("Multi-chain packing: Packing seed chain '%s' (length %zu) into %zu CLBs\n",
                seed_chain_name.c_str(), chain_length, batch->clb_ids.size());
    }

    bool seed_pack_success = true;
    for (size_t i = 0; i < chain_length && seed_pack_success; i++) {
        t_pack_molecule* mol = seed_chain_mols[i];
        LegalizationClusterId clb_id = batch->clb_ids[i];

        // Log CLB start for this batch CLB (using the chain molecule as the "seed")
        if (g_clustering_history_logger && g_clustering_history_logger->is_enabled()) {
            g_clustering_history_logger->log_clb_start(clb_id, chosen_type->name, mol);
        }

        e_block_pack_status status = cluster_legalizer.add_mol_to_cluster(mol, clb_id);
        if (status != e_block_pack_status::BLK_PASSED) {
            if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
                VTR_LOG("Multi-chain packing: Failed to pack seed chain molecule %zu into CLB %zu\n",
                        i, (size_t)clb_id);
            }
            seed_pack_success = false;
        } else {
            // Track molecule in batch for potential rollback
            batch->chain_molecules[seed_mol->chain_info.get()].push_back({i, mol});

            // Log seed placement details
            if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
                auto [num_tried, detail, attempts] = cluster_legalizer.get_last_molecule_failure_info();
                if (!attempts.empty()) {
                    VTR_LOG("Multi-chain packing:   Seed CLB[%zu] placed at: %s\n",
                            i, attempts[0].first.c_str());
                }
            }
        }
    }

    if (!seed_pack_success) {
        // Seed chain failed - destroy the entire batch and fall back
        num_used_type_instances[chosen_type] -= chain_length;
        cluster_legalizer.destroy_batch(*batch);
        cluster_legalizer.compress();

        // Fall back to regular single-cluster growth
        LegalizationClusterId cluster_id = try_grow_cluster(seed_mol,
                                                             strategy,
                                                             cluster_legalizer,
                                                             prepacker,
                                                             allow_unrelated_clustering,
                                                             balance_block_type_utilization,
                                                             timing_info,
                                                             clb_inter_blk_nets,
                                                             clustering_data,
                                                             attraction_groups,
                                                             num_used_type_instances,
                                                             mutable_device_ctx);
        if (cluster_id.is_valid()) {
            created_clusters.push_back(cluster_id);
        }
        return created_clusters;
    }

    // Find candidate chains that could potentially share these CLBs
    // They must:
    // 1. Be long chains
    // 2. Not yet clustered
    // 3. Have length <= seed chain length (so they fit in the batch)
    // 4. Be different from the seed chain
    std::vector<std::pair<t_pack_molecule*, size_t>> candidate_chains;  // (head_mol, chain_length)

    std::unordered_set<t_chain_info*> seen_chains;
    seen_chains.insert(seed_mol->chain_info.get());

    for (t_pack_molecule* mol = prepacker.get_molecules_vector().empty() ? nullptr : prepacker.get_molecules_vector()[0];
         mol != nullptr;
         mol = mol->next) {
        if (!mol->is_chain() || !mol->chain_info || !mol->chain_info->is_long_chain) {
            continue;
        }
        if (cluster_legalizer.is_mol_clustered(mol)) {
            continue;
        }
        if (seen_chains.count(mol->chain_info.get()) > 0) {
            continue;
        }

        size_t mol_chain_length = prepacker.calc_chain_length(mol);
        if (mol_chain_length > chain_length) {
            continue;  // Too long to fit
        }

        // Check if this chain's molecules are compatible with the cluster type
        bool compatible = true;
        auto chain_mols = prepacker.get_chain_molecules(mol);
        for (auto* chain_mol : chain_mols) {
            if (!cluster_legalizer.is_molecule_compatible(chain_mol, batch->clb_ids[0])) {
                compatible = false;
                break;
            }
        }

        if (compatible) {
            seen_chains.insert(mol->chain_info.get());
            candidate_chains.push_back({mol, mol_chain_length});
        }
    }

    // Sort candidates by length (longest first)
    std::sort(candidate_chains.begin(), candidate_chains.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1 && !candidate_chains.empty()) {
        VTR_LOG("Multi-chain packing: Found %zu candidate chains\n", candidate_chains.size());
    }

    // Try to pack each candidate chain into the batch
    for (const auto& [candidate_head, candidate_length] : candidate_chains) {
        auto candidate_mols = prepacker.get_chain_molecules(candidate_head);

        // Get candidate chain name for logging
        AtomBlockId candidate_root_atom = candidate_head->atom_block_ids[candidate_head->root];
        std::string candidate_name = atom_netlist_.block_name(candidate_root_atom);

        // Check if any molecule is already clustered
        bool any_clustered = false;
        for (auto* mol : candidate_mols) {
            if (cluster_legalizer.is_mol_clustered(mol)) {
                any_clustered = true;
                break;
            }
        }
        if (any_clustered) {
            if (MULTI_CHAIN_DEBUG && log_verbosity_ > 2) {
                VTR_LOG("Multi-chain packing: Skipping candidate '%s' (already clustered)\n",
                        candidate_name.c_str());
            }
            continue;
        }

        if (MULTI_CHAIN_DEBUG && log_verbosity_ > 2) {
            VTR_LOG("Multi-chain packing: Trying candidate chain '%s' (length %zu)\n",
                    candidate_name.c_str(), candidate_length);
        }

        // Try to pack all molecules in this candidate chain
        bool pack_success = true;
        std::vector<std::pair<size_t, t_pack_molecule*>> packed_mols;
        e_block_pack_status failure_status = e_block_pack_status::BLK_STATUS_UNDEFINED;
        size_t failure_mol_idx = 0;

        for (size_t i = 0; i < candidate_mols.size() && pack_success; i++) {
            t_pack_molecule* mol = candidate_mols[i];
            LegalizationClusterId clb_id = batch->clb_ids[i];

            e_block_pack_status status = cluster_legalizer.add_mol_to_cluster(mol, clb_id);
            if (status != e_block_pack_status::BLK_PASSED) {
                pack_success = false;
                failure_status = status;
                failure_mol_idx = i;
            } else {
                packed_mols.push_back({i, mol});
            }
        }

        if (pack_success) {
            // Successfully packed entire candidate chain
            batch->chain_molecules[candidate_head->chain_info.get()] = std::move(packed_mols);
            if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
                VTR_LOG("Multi-chain packing: Successfully packed candidate chain of length %zu\n",
                        candidate_length);
            }
        } else {
            // Log the failure reason
            if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
                const char* failure_reason = "unknown";
                switch (failure_status) {
                    case e_block_pack_status::BLK_FAILED_FEASIBLE:
                        failure_reason = "no feasible placement";
                        break;
                    case e_block_pack_status::BLK_FAILED_ROUTE:
                        failure_reason = "routing failed";
                        break;
                    case e_block_pack_status::BLK_FAILED_FLOORPLANNING:
                        failure_reason = "floorplanning constraint";
                        break;
                    case e_block_pack_status::BLK_FAILED_NOC_GROUP:
                        failure_reason = "NOC group mismatch";
                        break;
                    default:
                        break;
                }
                VTR_LOG("Multi-chain packing: Candidate '%s' failed at CLB[%zu]: %s\n",
                        candidate_name.c_str(), failure_mol_idx, failure_reason);

                // Get detailed failure info including placement attempts
                auto [num_tried, detail, attempts] = cluster_legalizer.get_last_molecule_failure_info();
                VTR_LOG("Multi-chain packing:   -> %d placements tried, detail: %s\n",
                        num_tried, detail.c_str());

                // Log each placement attempt
                for (size_t attempt_idx = 0; attempt_idx < attempts.size(); attempt_idx++) {
                    const auto& [prim_path, fail_reason] = attempts[attempt_idx];
                    VTR_LOG("Multi-chain packing:   -> Attempt %zu: tried '%s'",
                            attempt_idx + 1, prim_path.c_str());
                    if (!fail_reason.empty()) {
                        VTR_LOG(" - FAILED: %s", fail_reason.c_str());
                    }
                    VTR_LOG("\n");
                }
            }

            // Rollback partially packed molecules
            for (const auto& [clb_idx, mol] : packed_mols) {
                // Add to batch tracking so rollback can find them
                batch->chain_molecules[candidate_head->chain_info.get()].push_back({clb_idx, mol});
            }
            cluster_legalizer.rollback_chain_from_batch(*batch, candidate_head->chain_info.get());
        }
    }

    // Fill remaining space in each CLB with regular (non-chain) molecules
    int high_fanout_threshold = high_fanout_thresholds_.get_threshold(chosen_type->name);

    for (size_t clb_idx = 0; clb_idx < batch->clb_ids.size(); clb_idx++) {
        LegalizationClusterId clb_id = batch->clb_ids[clb_idx];
        if (!clb_id.is_valid()) continue;

        // Update cluster stats for the seed molecule in this CLB (for gain calculation)
        t_pack_molecule* seed_in_clb = seed_chain_mols[clb_idx];
        update_cluster_stats(seed_in_clb,
                             cluster_legalizer,
                             is_clock_,
                             is_global_,
                             packer_opts_.global_clocks,
                             packer_opts_.alpha, packer_opts_.beta,
                             packer_opts_.timing_driven, packer_opts_.connection_driven,
                             high_fanout_threshold,
                             timing_info,
                             attraction_groups,
                             net_output_feeds_driving_block_input_);

        // Get candidate molecules for filling
        int num_unrelated_clustering_attempts = 0;
        t_pack_molecule* candidate_mol = get_molecule_for_cluster(
            cluster_legalizer.get_cluster_pb(clb_id),
            attraction_groups,
            allow_unrelated_clustering,
            packer_opts_.prioritize_transitive_connectivity,
            packer_opts_.transitive_fanout_threshold,
            packer_opts_.feasible_block_array_size,
            &num_unrelated_clustering_attempts,
            prepacker,
            cluster_legalizer,
            clb_inter_blk_nets,
            clb_id,
            log_verbosity_,
            clustering_data.unclustered_list_head,
            clustering_data.unclustered_list_head_size,
            primitive_candidate_block_types_);

        int max_num_repeated_molecules = 1;
        if (attraction_groups.num_attraction_groups() > 0)
            max_num_repeated_molecules = attraction_groups_max_repeated_molecules_;

        int num_repeated_molecules = 0;
        while (candidate_mol != nullptr && num_repeated_molecules < max_num_repeated_molecules) {
            // Skip chain molecules (we handle them separately)
            if (candidate_mol->is_chain() && candidate_mol->chain_info &&
                candidate_mol->chain_info->is_long_chain) {
                // Get next candidate without incrementing repeated counter
                candidate_mol = get_molecule_for_cluster(
                    cluster_legalizer.get_cluster_pb(clb_id),
                    attraction_groups,
                    allow_unrelated_clustering,
                    packer_opts_.prioritize_transitive_connectivity,
                    packer_opts_.transitive_fanout_threshold,
                    packer_opts_.feasible_block_array_size,
                    &num_unrelated_clustering_attempts,
                    prepacker,
                    cluster_legalizer,
                    clb_inter_blk_nets,
                    clb_id,
                    log_verbosity_,
                    clustering_data.unclustered_list_head,
                    clustering_data.unclustered_list_head_size,
                    primitive_candidate_block_types_);
                continue;
            }

            bool success = try_add_candidate_mol_to_cluster(candidate_mol, clb_id, cluster_legalizer);

            if (success) {
                update_cluster_stats(candidate_mol,
                                     cluster_legalizer,
                                     is_clock_,
                                     is_global_,
                                     packer_opts_.global_clocks,
                                     packer_opts_.alpha,
                                     packer_opts_.beta,
                                     packer_opts_.timing_driven,
                                     packer_opts_.connection_driven,
                                     high_fanout_threshold,
                                     timing_info,
                                     attraction_groups,
                                     net_output_feeds_driving_block_input_);
                num_unrelated_clustering_attempts = 0;
            }

            t_pack_molecule* prev_candidate_mol = candidate_mol;
            candidate_mol = get_molecule_for_cluster(
                cluster_legalizer.get_cluster_pb(clb_id),
                attraction_groups,
                allow_unrelated_clustering,
                packer_opts_.prioritize_transitive_connectivity,
                packer_opts_.transitive_fanout_threshold,
                packer_opts_.feasible_block_array_size,
                &num_unrelated_clustering_attempts,
                prepacker,
                cluster_legalizer,
                clb_inter_blk_nets,
                clb_id,
                log_verbosity_,
                clustering_data.unclustered_list_head,
                clustering_data.unclustered_list_head_size,
                primitive_candidate_block_types_);

            if (candidate_mol == prev_candidate_mol)
                num_repeated_molecules++;
        }
    }

    // Check legality of all CLBs if using SKIP_INTRA_LB_ROUTE strategy
    bool all_legal = true;
    if (strategy != ClusterLegalizationStrategy::FULL) {
        for (LegalizationClusterId clb_id : batch->clb_ids) {
            if (!clb_id.is_valid()) continue;
            if (!cluster_legalizer.check_cluster_legality(clb_id)) {
                all_legal = false;
                if (MULTI_CHAIN_DEBUG && log_verbosity_ > 1) {
                    VTR_LOG("Multi-chain packing: CLB %zu failed final legality check\n",
                            (size_t)clb_id);
                }
            }
        }
    }

    if (!all_legal) {
        // Batch failed - destroy and fall back
        num_used_type_instances[chosen_type] -= chain_length;
        cluster_legalizer.destroy_batch(*batch);
        cluster_legalizer.compress();

        // Try again with FULL strategy
        if (strategy != ClusterLegalizationStrategy::FULL) {
            return try_grow_multi_chain_cluster(seed_mol,
                                                 ClusterLegalizationStrategy::FULL,
                                                 cluster_legalizer,
                                                 prepacker,
                                                 allow_unrelated_clustering,
                                                 balance_block_type_utilization,
                                                 timing_info,
                                                 clb_inter_blk_nets,
                                                 clustering_data,
                                                 attraction_groups,
                                                 num_used_type_instances,
                                                 mutable_device_ctx);
        }
        return created_clusters;  // Empty
    }

    // Store inter-block nets for each CLB (for gain calculation)
    for (LegalizationClusterId clb_id : batch->clb_ids) {
        if (!clb_id.is_valid()) continue;
        t_pb* cur_pb = cluster_legalizer.get_cluster_pb(clb_id);
        t_pb_stats* pb_stats = cur_pb->pb_stats;
        for (const AtomNetId mnet_id : pb_stats->marked_nets) {
            int external_terminals = atom_netlist_.net_pins(mnet_id).size() - pb_stats->num_pins_of_net_in_pb[mnet_id];
            if (external_terminals < packer_opts_.transitive_fanout_threshold && external_terminals > 0) {
                clb_inter_blk_nets[clb_id].push_back(mnet_id);
            }
        }
    }

    // Log CLB creation success for each CLB in the batch
    if (g_clustering_history_logger && g_clustering_history_logger->is_enabled()) {
        for (LegalizationClusterId clb_id : batch->clb_ids) {
            if (!clb_id.is_valid()) continue;
            g_clustering_history_logger->log_clb_success(
                clb_id,
                cluster_legalizer.get_cluster_pb(clb_id),
                cluster_legalizer.get_cluster_molecules(clb_id),
                strategy);
        }
    }

    // Finalize the batch (this calls clean_cluster on each CLB)
    cluster_legalizer.finalize_batch(*batch);

    // Record finalized CLBs for profile summary
    if (g_clustering_history_logger && g_clustering_history_logger->is_profile_enabled()) {
        for (LegalizationClusterId clb_id : batch->clb_ids) {
            if (!clb_id.is_valid()) continue;
            g_clustering_history_logger->record_finalized_clb(
                clb_id,
                cluster_legalizer.get_cluster_pb(clb_id),
                cluster_legalizer.get_cluster_molecules(clb_id));
        }
    }

    // Return all created cluster IDs
    for (LegalizationClusterId clb_id : batch->clb_ids) {
        if (clb_id.is_valid()) {
            created_clusters.push_back(clb_id);
        }
    }

    if (MULTI_CHAIN_DEBUG && log_verbosity_ > 0) {
        VTR_LOG("Multi-chain packing: Created %zu CLBs with %zu chains\n",
                created_clusters.size(), batch->chain_molecules.size());
    }

    return created_clusters;
}

LegalizationClusterId GreedyClusterer::start_new_cluster(
    t_pack_molecule* seed_mol,
    ClusterLegalizer& cluster_legalizer,
    bool balance_block_type_utilization,
    std::map<t_logical_block_type_ptr, size_t>& num_used_type_instances,
    DeviceContext& mutable_device_ctx) {
    /* Allocate a dummy initial cluster and load a atom block as a seed and check if it is legal */
    AtomBlockId root_atom = seed_mol->atom_block_ids[seed_mol->root];
    const std::string& root_atom_name = atom_netlist_.block_name(root_atom);
    const t_model* root_model = atom_netlist_.block_model(root_atom);

    auto itr = primitive_candidate_block_types_.find(root_model);
    VTR_ASSERT(itr != primitive_candidate_block_types_.end());
    std::vector<t_logical_block_type_ptr> candidate_types = itr->second;

    if (balance_block_type_utilization) {
        //We sort the candidate types in ascending order by their current utilization.
        //This means that the packer will prefer to use types with lower utilization.
        //This is a naive approach to try balancing utilization when multiple types can
        //support the same primitive(s).
        std::stable_sort(candidate_types.begin(), candidate_types.end(),
                         [&](t_logical_block_type_ptr lhs, t_logical_block_type_ptr rhs) {
                             int lhs_num_instances = 0;
                             int rhs_num_instances = 0;
                             // Count number of instances for each type
                             for (auto type : lhs->equivalent_tiles)
                                 lhs_num_instances += mutable_device_ctx.grid.num_instances(type, -1);
                             for (auto type : rhs->equivalent_tiles)
                                 rhs_num_instances += mutable_device_ctx.grid.num_instances(type, -1);

                             float lhs_util = vtr::safe_ratio<float>(num_used_type_instances[lhs], lhs_num_instances);
                             float rhs_util = vtr::safe_ratio<float>(num_used_type_instances[rhs], rhs_num_instances);
                             //Lower util first
                             return lhs_util < rhs_util;
                         });
    }

    //Try packing into each candidate type
    bool success = false;
    t_logical_block_type_ptr block_type;
    LegalizationClusterId new_cluster_id;
    for (auto type : candidate_types) {
        //Try packing into each mode
        e_block_pack_status pack_result = e_block_pack_status::BLK_STATUS_UNDEFINED;
        for (int j = 0; j < type->pb_graph_head->pb_type->num_modes && !success; j++) {
            std::tie(pack_result, new_cluster_id) = cluster_legalizer.start_new_cluster(seed_mol, type, j);
            success = (pack_result == e_block_pack_status::BLK_PASSED);
        }

        if (success) {
            // VTR_LOG("PASSED_SEED: Block Type %s\n", type->name.c_str());
            // If clustering succeeds return the new_cluster_id and type.
            block_type = type;
            break;
        } else {
            // VTR_LOG("FAILED_SEED: Block Type %s\n", type->name.c_str());
        }
    }

    if (!success) {
        //Explored all candidates
        if (seed_mol->type == MOLECULE_FORCED_PACK) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "Can not find any logic block that can implement molecule.\n"
                            "\tPattern %s %s\n",
                            seed_mol->pack_pattern->name,
                            root_atom_name.c_str());
        } else {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "Can not find any logic block that can implement molecule.\n"
                            "\tAtom %s (%s)\n",
                            root_atom_name.c_str(), root_model->name);
        }
    }

    VTR_ASSERT(success);
    VTR_ASSERT(new_cluster_id.is_valid());
    //Progress dot for seed-block
    fflush(stdout);

    // TODO: Below may make more sense in its own method.

    // Successfully created cluster
    num_used_type_instances[block_type]++;

    /* Expand FPGA size if needed */
    // Check used type instances against the possible equivalent physical locations
    unsigned int num_instances = 0;
    for (auto equivalent_tile : block_type->equivalent_tiles) {
        num_instances += mutable_device_ctx.grid.num_instances(equivalent_tile, -1);
    }

    if (num_used_type_instances[block_type] > num_instances) {
        mutable_device_ctx.grid = create_device_grid(packer_opts_.device_layout,
                                                     arch_.grid_layouts,
                                                     num_used_type_instances,
                                                     packer_opts_.target_device_utilization);
    }

    return new_cluster_id;
}

bool GreedyClusterer::try_add_candidate_mol_to_cluster(t_pack_molecule* candidate_mol,
                                                       LegalizationClusterId legalization_cluster_id,
                                                       ClusterLegalizer& cluster_legalizer) {
    // VTR_LOG("try_add_candidate_mol_to_cluster: entered\n");
    VTR_ASSERT(candidate_mol != nullptr);
    // VTR_LOG("try_add_candidate_mol_to_cluster: checking is_mol_clustered\n");
    VTR_ASSERT(!cluster_legalizer.is_mol_clustered(candidate_mol));
    // VTR_LOG("try_add_candidate_mol_to_cluster: checking legalization_cluster_id\n");
    VTR_ASSERT(legalization_cluster_id.is_valid());

    // VTR_LOG("try_add_candidate_mol_to_cluster: accessing atom_block_ids\n");
    AtomBlockId blk_id = candidate_mol->atom_block_ids[candidate_mol->root];
    // VTR_LOG("try_add_candidate_mol_to_cluster: checking blk_id valid\n");
    VTR_ASSERT(blk_id.is_valid());
    // VTR_LOG("try_add_candidate_mol_to_cluster: getting block name\n");
    std::string blk_name = atom_netlist_.block_name(blk_id);
    // VTR_LOG("try_add_candidate_mol_to_cluster: copied block name: %s\n", blk_name.c_str());

    // VTR_LOG("try_add_candidate_mol_to_cluster: checking pack_pattern %p\n", candidate_mol->pack_pattern);
    // if (candidate_mol->pack_pattern) {
    //     VTR_LOG("try_add_candidate_mol_to_cluster: pack_pattern name: %s\n", candidate_mol->pack_pattern->name);
    // } else {
    //     VTR_LOG("try_add_candidate_mol_to_cluster: pack_pattern is NULL\n");
    // }

    // VTR_LOG("Attempting to add candidate molecule %s (root: %s) to cluster\n",
    //         (candidate_mol->pack_pattern ? candidate_mol->pack_pattern->name : "NULL"),
    //         blk_name.c_str());

    e_block_pack_status pack_status = cluster_legalizer.add_mol_to_cluster(candidate_mol,
                                                                           legalization_cluster_id);

    // Record candidate failure for logging
    if (pack_status != e_block_pack_status::BLK_PASSED) {
        if (g_clustering_history_logger) {
            AtomBlockId blk_id = candidate_mol->atom_block_ids[candidate_mol->root];
            const t_model* blk_model = atom_netlist_.block_model(blk_id);
            g_clustering_history_logger->record_candidate_failure(pack_status, blk_model->name);
        }
    }

    // Print helpful debugging log messages.
    if (log_verbosity_ > 2) {
        switch (pack_status) {
            case e_block_pack_status::BLK_PASSED:
                VTR_LOG("\tPassed: ");
                break;
            case e_block_pack_status::BLK_FAILED_ROUTE:
                VTR_LOG("\tNO_ROUTE: ");
                break;
            case e_block_pack_status::BLK_FAILED_FLOORPLANNING:
                VTR_LOG("\tFAILED_FLOORPLANNING_CONSTRAINTS_CHECK: ");
                break;
            case e_block_pack_status::BLK_FAILED_FEASIBLE:
                VTR_LOG("\tFAILED_FEASIBILITY_CHECK: ");
                break;
            case e_block_pack_status::BLK_FAILED_NOC_GROUP:
                VTR_LOG("\tFAILED_NOC_GROUP_CHECK: ");
                break;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_PACK, "Unknown pack status thrown.");
                break;
        }
        // Get the block name and model name
        AtomBlockId blk_id = candidate_mol->atom_block_ids[candidate_mol->root];
        VTR_ASSERT(blk_id.is_valid());
        std::string blk_name = atom_netlist_.block_name(blk_id);
        const t_model* blk_model = atom_netlist_.block_model(blk_id);
        VTR_LOG("'%s' (%s)", blk_name.c_str(), blk_model->name);
        VTR_LOGV(candidate_mol->pack_pattern, " molecule %s molecule_size %zu",
                 candidate_mol->pack_pattern->name,
                 candidate_mol->atom_block_ids.size());
        VTR_LOG("\n");
        fflush(stdout);
    }

    return pack_status == e_block_pack_status::BLK_PASSED;
}

void GreedyClusterer::report_le_physical_block_usage(const ClusterLegalizer& cluster_legalizer) {
    // find the cluster type that has lut primitives
    auto logic_block_type = identify_logic_block_type(primitive_candidate_block_types_);
    // find a LE pb_type within the found logic_block_type
    auto le_pb_type = identify_le_block_type(logic_block_type);

    // If this architecture does not have an LE physical block, cannot report
    // its usage.
    if (le_pb_type == nullptr)
        return;

    // Track the number of Logic Elements (LEs) used. This is populated only for
    // architectures which has LEs. The architecture is assumed to have LEs iff
    // it has a logic block that contains LUT primitives and is the first
    // pb_block to have more than one instance from the top of the hierarchy
    // (All parent pb_block have one instance only and one mode only).

    // The number of LEs that are used for logic (LUTs/adders) only.
    int num_logic_le = 0;
    // The number of LEs that are used for registers only.
    int num_reg_le = 0;
    // The number of LEs that are used for both logic (LUTs/adders) and registers.
    int num_logic_and_reg_le = 0;

    for (LegalizationClusterId cluster_id : cluster_legalizer.clusters()) {
        // Update the data structure holding the LE counts
        update_le_count(cluster_legalizer.get_cluster_pb(cluster_id),
                        logic_block_type,
                        le_pb_type,
                        num_logic_le,
                        num_reg_le,
                        num_logic_and_reg_le);
    }

    // if this architecture has LE physical block, report its usage
    if (le_pb_type) {
        print_le_count(num_logic_le, num_reg_le, num_logic_and_reg_le, le_pb_type);
    }
}
