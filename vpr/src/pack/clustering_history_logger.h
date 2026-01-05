/**
 * @file
 * @brief   Logging of CLB creation history during clustering for debugging.
 *
 * This file provides the ClusteringHistoryLogger class which logs detailed
 * information about CLB creation attempts, including:
 * - Route type attempted (SKIP_INTRA_LB_ROUTE vs FULL)
 * - Success/failure status
 * - Molecules packed and their placements
 * - Congestion details on routing failures
 * - Timing information for each CLB creation
 */

#ifndef CLUSTERING_HISTORY_LOGGER_H
#define CLUSTERING_HISTORY_LOGGER_H

#include <chrono>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "atom_netlist_fwd.h"
#include "cluster_legalizer.h"

/**
 * @brief Active mode categories for FLE utilization tracking.
 *
 * Each FLE can have multiple active modes simultaneously (e.g., LUT5 + chain).
 */
enum class FleActiveMode {
    LUT5,          ///< LUT5 primitive is being used
    SIMPLE_CHAIN,  ///< Atoms belong to a simple_chain molecule
    CHAIN          ///< Atoms belong to a chain molecule (not simple_chain)
};

// Forward declarations
class t_pack_molecule;
struct t_pb;
struct t_pb_graph_node;
struct t_lb_router_data;
struct t_lb_type_rr_node;
struct t_lb_rr_node_stats;

// PlacementAttemptInfo is defined in cluster_legalizer.h (included above)

/**
 * @brief Logger for recording CLB creation history during clustering.
 *
 * This class provides methods to log detailed information about the clustering
 * process, including molecule packing attempts, routing paths used, and
 * congestion information on failures.
 */
class ClusteringHistoryLogger {
public:
    /// Recursive structure for pb hierarchy within a BLE5
    struct PbNodeInfo {
        std::string pb_type_name;  ///< Name of this pb_type (e.g., "ble5", "arithmetic_1chain", "adder")
        int pb_index = -1;         ///< Index within parent
        std::string mode;          ///< Mode selected (empty string for primitives)
        bool is_primitive = false; ///< True if this is a primitive (leaf node)
        std::string atom_name;     ///< For primitives, the atom placed here
        std::map<std::string, std::vector<std::string>> input_pins;  ///< Pin name -> list of targets it drives
        std::vector<PbNodeInfo> children;  ///< Child pb nodes
    };

    /// Information about a BLE5's pin utilization (now includes full hierarchy)
    struct Ble5Utilization {
        int ble5_index = -1;
        std::string ble5_mode;  ///< Mode selected (e.g., "blut5", "arithmetic")
        std::vector<std::string> atoms;  ///< Atoms placed in this BLE5 (flat list for quick reference)
        std::map<std::string, std::vector<std::string>> input_pins;  ///< BLE5-level pin usage
        std::vector<PbNodeInfo> hierarchy;  ///< Full pb hierarchy within this BLE5
    };

    /// Information about a FLE's BLE5 utilization
    struct FleUtilization {
        int fle_index = -1;
        std::string fle_mode;  ///< Mode selected for this FLE (e.g., "n2_lut5", "arithmetic")
        std::map<int, Ble5Utilization> ble5_usage;  ///< BLE5 index -> utilization
    };

    /**
     * @brief Construct a new ClusteringHistoryLogger.
     *
     * Opens the echo file for writing if echo is enabled.
     */
    ClusteringHistoryLogger();

    /**
     * @brief Destructor - closes the log file if open.
     */
    ~ClusteringHistoryLogger();

    // Prevent copying
    ClusteringHistoryLogger(const ClusteringHistoryLogger&) = delete;
    ClusteringHistoryLogger& operator=(const ClusteringHistoryLogger&) = delete;

    /**
     * @brief Check if history logging is enabled.
     * @return true if the history echo file is enabled and open.
     */
    bool is_enabled() const { return file_.is_open(); }

    /**
     * @brief Check if profile logging is enabled.
     * @return true if the profile echo file is enabled and open.
     */
    bool is_profile_enabled() const { return profile_file_.is_open(); }

    /**
     * @brief Log the start of a new packing iteration.
     *
     * Call this when VPR starts a new packing attempt (e.g., after a failed iteration).
     *
     * @param iteration     The iteration number (1-based).
     * @param strategy      The strategy being used for this iteration.
     */
    void log_iteration_start(int iteration, ClusterLegalizationStrategy strategy);

    /**
     * @brief Log the start of a new CLB creation attempt.
     *
     * @param cluster_id        The ID of the cluster being created.
     * @param cluster_type_name The name of the cluster type (e.g., "clb").
     * @param seed_molecule     The seed molecule used to start the cluster.
     */
    void log_clb_start(LegalizationClusterId cluster_id,
                       const std::string& cluster_type_name,
                       const t_pack_molecule* seed_molecule);

    /**
     * @brief Log a molecule packing attempt.
     *
     * @param molecule          The molecule being packed.
     * @param strategy          The legalization strategy (SKIP or FULL).
     * @param status            The result of the packing attempt.
     * @param primitives_list   The primitive placements for atoms in the molecule.
     * @param molecule_size     The number of atoms in the molecule.
     * @param elapsed_us        Time elapsed for this attempt in microseconds.
     */
    void log_molecule_attempt(const t_pack_molecule* molecule,
                              ClusterLegalizationStrategy strategy,
                              e_block_pack_status status,
                              t_pb_graph_node** primitives_list,
                              int molecule_size,
                              double elapsed_us);

    /**
     * @brief Log routing failure with congestion details.
     *
     * @param molecule          The molecule that failed routing.
     * @param router_data       The router data containing congestion info.
     * @param primitives_list   The primitive placements attempted.
     * @param molecule_size     The number of atoms in the molecule.
     */
    void log_routing_failure(const t_pack_molecule* molecule,
                             const t_lb_router_data* router_data,
                             t_pb_graph_node** primitives_list,
                             int molecule_size);

    /**
     * @brief Log detailed routing paths for all nets in a failed routing attempt.
     *
     * For each net, this traces the routing tree from source to sinks,
     * showing which RR nodes (pins) are used. This helps debug routing
     * congestion by showing exactly how signals are routed through the CLB.
     *
     * Format for each net:
     *   NET: <net_name>
     *     source: <atom_pin> @ <pb_pin>
     *       -> <intermediate_pin>
     *       -> <intermediate_pin>
     *       -> sink: <atom_pin> @ <pb_pin>
     *
     * @param router_data  The router data containing nets and routing trees.
     */
    void log_routing_paths(const t_lb_router_data* router_data);

    /**
     * @brief Log successful CLB finalization.
     *
     * @param cluster_id    The ID of the finalized cluster.
     * @param cluster_pb    The pb structure of the cluster.
     * @param molecules     The molecules packed into the cluster.
     * @param strategy      The legalization strategy used.
     */
    void log_clb_success(LegalizationClusterId cluster_id,
                         const t_pb* cluster_pb,
                         const std::vector<t_pack_molecule*>& molecules,
                         ClusterLegalizationStrategy strategy);

    /**
     * @brief Log CLB creation failure.
     *
     * @param cluster_id    The ID of the failed cluster attempt.
     * @param reason        Human-readable reason for failure.
     * @param strategy      The legalization strategy used.
     * @param cluster_pb    Optional: The pb structure of the cluster (to show attempted placements).
     * @param molecules     Optional: The molecules that were packed before failure.
     */
    void log_clb_failure(LegalizationClusterId cluster_id,
                         const std::string& reason,
                         ClusterLegalizationStrategy strategy,
                         const t_pb* cluster_pb = nullptr,
                         const std::vector<t_pack_molecule*>* molecules = nullptr);

    /**
     * @brief Log timing for the current CLB creation.
     *
     * Call this when done with a CLB to record elapsed time.
     */
    void log_clb_timing();

    /**
     * @brief Log details about feasibility failure for a molecule.
     *
     * @param molecule              The molecule that failed to pack.
     * @param num_placements_tried  How many primitive placements were attempted.
     * @param last_failure_reason   Description of the last failure reason.
     */
    void log_feasibility_failure(const t_pack_molecule* molecule,
                                  int num_placements_tried,
                                  const std::string& last_failure_reason);

    /**
     * @brief Log details about feasibility failure for a molecule with placement attempts.
     *
     * @param molecule              The molecule that failed to pack.
     * @param num_placements_tried  How many primitive placements were attempted.
     * @param last_failure_reason   Description of the last failure reason.
     * @param placement_attempts    Details of each primitive placement attempted.
     */
    void log_feasibility_failure(const t_pack_molecule* molecule,
                                  int num_placements_tried,
                                  const std::string& last_failure_reason,
                                  const std::vector<PlacementAttemptInfo>& placement_attempts);

    /**
     * @brief Get primitive placement description (hierarchical path).
     *
     * Returns a string like "clb[0]/fle[0]/ble5[0]/lut5[0]" describing
     * the location of a primitive in the pb_graph hierarchy.
     *
     * @param primitive  The pb_graph_node to describe.
     * @return           Hierarchical path string.
     */
    std::string get_placement_description(t_pb_graph_node* primitive) const;

    /**
     * @brief Record a candidate molecule packing failure.
     *
     * @param status          The failure status (reason).
     * @param primitive_type  The type of primitive (e.g., "lut", "ff", "adder").
     */
    void record_candidate_failure(e_block_pack_status status, const std::string& primitive_type);

    /**
     * @brief Reset the candidate failure stats for a new CLB.
     */
    void reset_candidate_stats();

    /**
     * @brief Log the accumulated candidate failure statistics.
     */
    void log_candidate_failure_stats();

    /// Source of routing attempt
    enum class RoutingSource {
        MOLECULE,   ///< From try_pack_molecule (per-molecule routing)
        LEGALITY    ///< From check_cluster_legality (final check)
    };

    /**
     * @brief Record statistics from a single intra-LB routing attempt.
     *
     * Uses the current routing source set via set_routing_source().
     *
     * @param lb_rr_graph_size   Number of nodes in the LB RR graph.
     * @param num_nets           Number of nets routed.
     * @param num_iterations     Number of pathfinder iterations used.
     * @param elapsed_us         Time elapsed in microseconds.
     * @param success            Whether routing succeeded.
     * @param is_impossible      Whether routing hit impossible state.
     */
    void record_routing_attempt(size_t lb_rr_graph_size,
                                 size_t num_nets,
                                 int num_iterations,
                                 double elapsed_us,
                                 bool success,
                                 bool is_impossible);

    /**
     * @brief Record a mode retry event.
     */
    void record_mode_retry();

    /**
     * @brief Set the current routing source for subsequent routing attempts.
     */
    void set_routing_source(RoutingSource source) { current_routing_source_ = source; }

    /**
     * @brief Get the current routing source.
     */
    RoutingSource get_routing_source() const { return current_routing_source_; }

    /**
     * @brief Set the current molecule type being routed.
     *
     * Call this before routing to track which primitive types are failing.
     * @param molecule_type  The type name (e.g., "adder", "latch", "lut6")
     */
    void set_current_molecule_type(const std::string& molecule_type) { current_molecule_type_ = molecule_type; }

    /**
     * @brief Log the accumulated routing statistics.
     */
    void log_routing_stats();

    /**
     * @brief Record a finalized CLB for the end-of-clustering summary.
     *
     * This stores CLB information including molecules, atom placements, and
     * pin utilization for output in the final summary.
     *
     * @param cluster_id    The ID of the finalized cluster.
     * @param cluster_pb    The pb structure of the cluster.
     * @param molecules     The molecules packed into the cluster.
     */
    void record_finalized_clb(LegalizationClusterId cluster_id,
                              const t_pb* cluster_pb,
                              const std::vector<t_pack_molecule*>& molecules);

    /**
     * @brief Write the final summary of all finalized CLBs.
     *
     * Call this at the end of clustering to output a summary showing all
     * CLBs with their molecules, atom placements, and pin utilization.
     */
    void write_summary();

private:
    /// Information about a molecule and its atom placements
    struct MoleculeInfo {
        std::string root_atom_name;   ///< Name of the root atom
        std::string pattern_name;     ///< Pack pattern name (empty if single atom)
        int num_blocks = 0;           ///< Number of blocks in molecule
        std::vector<std::pair<std::string, std::string>> atom_placements;  ///< (atom_name, placement_path)
    };

    /// Information about a finalized CLB for the summary output
    struct FinalizedClbInfo {
        LegalizationClusterId cluster_id;
        std::string cluster_name;
        std::string cluster_type;
        std::vector<MoleculeInfo> molecules;  ///< Molecules with their atoms grouped
        std::map<int, FleUtilization> fle_utilization;  ///< FLE index -> utilization info
        int total_fles = 0;
        int used_fles = 0;
    };

    /// Storage for all finalized CLBs
    std::vector<FinalizedClbInfo> finalized_clbs_;

    /// The output file stream for history (real-time logging)
    std::ofstream file_;

    /// The output file stream for profile (finalized CLB summary)
    std::ofstream profile_file_;

    /// Timer for tracking CLB creation time
    std::chrono::high_resolution_clock::time_point clb_start_time_;

    /// Counter for CLB attempts
    size_t clb_count_ = 0;

    /// Current routing source (set by caller before routing)
    RoutingSource current_routing_source_ = RoutingSource::MOLECULE;

    /// Current molecule type being routed (set by caller before routing)
    std::string current_molecule_type_;

    /// Candidate failure stats: reason -> (primitive_type -> count)
    std::map<e_block_pack_status, std::map<std::string, int>> candidate_failure_stats_;

    /// Routing statistics for the current CLB
    struct RoutingStats {
        size_t lb_rr_graph_size = 0;      ///< Size of the LB RR graph (set once)
        size_t total_attempts = 0;         ///< Total routing attempts
        size_t successful_attempts = 0;    ///< Successful routing attempts
        size_t failed_attempts = 0;        ///< Failed routing attempts
        size_t impossible_attempts = 0;    ///< Attempts that hit impossible state
        size_t total_nets_routed = 0;      ///< Total nets across all attempts
        int total_iterations = 0;          ///< Total pathfinder iterations
        int max_iterations = 0;            ///< Max iterations in single attempt
        double total_time_us = 0.0;        ///< Total routing time in microseconds
        double max_time_us = 0.0;          ///< Max time for single attempt
        // Source tracking
        size_t attempts_from_molecule = 0;     ///< From try_pack_molecule
        size_t attempts_from_legality = 0;     ///< From check_cluster_legality
        size_t mode_retry_count = 0;           ///< Mode conflict retries
        // Failure breakdown by primitive type
        std::map<std::string, size_t> failed_by_type;      ///< Failed attempts by molecule type
        std::map<std::string, size_t> impossible_by_type;  ///< Impossible attempts by molecule type
    };
    RoutingStats routing_stats_;

    /**
     * @brief Get a string representation of the pack status.
     */
    static std::string status_to_string(e_block_pack_status status);

    /**
     * @brief Get a string representation of the legalization strategy.
     */
    static std::string strategy_to_string(ClusterLegalizationStrategy strategy);

    /**
     * @brief Get atom name from molecule at given index.
     */
    std::string get_atom_name(const t_pack_molecule* molecule, int index) const;

    /**
     * @brief Get primitive placement description with mode information.
     *
     * This version walks up the t_pb hierarchy to show which mode was
     * selected at each level (e.g., "ble5[0][arithmetic_2chains]").
     */
    std::string get_placement_description_with_mode(const t_pb* atom_pb) const;

    /**
     * @brief Get congestion description for routing nodes.
     */
    std::string describe_congestion(const t_lb_router_data* router_data) const;
};

/**
 * @brief Global instance of the clustering history logger.
 *
 * This is initialized when the ClusterLegalizer is constructed and can be
 * accessed from the packing code to log events.
 */
extern ClusteringHistoryLogger* g_clustering_history_logger;

/**
 * @brief Build a detailed description of why pin feasibility check failed.
 *
 * This function examines the lookahead pin usage in a cluster and returns
 * a detailed string describing which pin classes exceeded their capacity
 * and which nets are competing for those pins.
 *
 * Similar to describe_congestion() for routing failures, this provides
 * visibility into why the pin feasibility filter rejected a packing attempt.
 *
 * @param cur_pb               The pb (cluster) that failed the pin feasibility check.
 * @param max_external_pin_util The maximum external pin utilization factor used.
 * @return                     A formatted string describing the failure details.
 */
std::string describe_pin_feasibility_failure(const t_pb* cur_pb, t_ext_pin_util max_external_pin_util);

#endif // CLUSTERING_HISTORY_LOGGER_H
