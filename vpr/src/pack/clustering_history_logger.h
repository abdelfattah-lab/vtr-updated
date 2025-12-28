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

// Forward declarations
class t_pack_molecule;
struct t_pb;
struct t_pb_graph_node;
struct t_lb_router_data;
struct t_lb_type_rr_node;
struct t_lb_rr_node_stats;

/**
 * @brief Logger for recording CLB creation history during clustering.
 *
 * This class provides methods to log detailed information about the clustering
 * process, including molecule packing attempts, routing paths used, and
 * congestion information on failures.
 */
class ClusteringHistoryLogger {
public:
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
     * @brief Check if logging is enabled.
     * @return true if the echo file is enabled and open.
     */
    bool is_enabled() const { return file_.is_open(); }

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
     * @param all_placement_attempts  All placement attempts with their results.
     */
    void log_feasibility_failure(const t_pack_molecule* molecule,
                                  int num_placements_tried,
                                  const std::string& last_failure_reason,
                                  const std::vector<std::string>& all_placement_attempts);

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

private:
    /// The output file stream
    std::ofstream file_;

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
     * @brief Get primitive placement description (without mode info).
     */
    std::string get_placement_description(t_pb_graph_node* primitive) const;

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

#endif // CLUSTERING_HISTORY_LOGGER_H
