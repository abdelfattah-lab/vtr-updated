/**
 * @file cluster_profiler.h
 * @brief Profiling utilities for clustering performance analysis
 *
 * This file provides a ClusterProfiler class that tracks timing and statistics
 * for various phases of the clustering algorithm to help identify bottlenecks.
 */

#pragma once

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>
#include <map>

/**
 * @brief Statistics for a single cluster's creation
 */
struct ClusterStats {
    size_t cluster_id = 0;
    std::string cluster_name;
    std::string cluster_type;

    // Timing (in microseconds)
    double total_time_us = 0;
    double start_cluster_time_us = 0;
    double molecule_selection_time_us = 0;
    double molecule_packing_time_us = 0;
    double stats_update_time_us = 0;
    double legality_check_time_us = 0;

    // Counts
    int molecules_tried = 0;
    int molecules_packed = 0;
    int molecules_failed = 0;
    int chain_molecules_tried = 0;
    int chain_molecules_packed = 0;

    // Packing attempt details
    int total_primitive_candidates = 0;
    int early_chain_rejects = 0;
    int placement_failures = 0;
    int routing_failures = 0;
    int expansion_failures = 0;  // Failed in expand_forced_pack_molecule_placement

    // Per-molecule stats within this cluster
    double avg_primitives_per_molecule = 0;
    double avg_pack_attempts_per_molecule = 0;

    // Primitive type breakdown
    int chain_primitive_candidates = 0;
    int non_chain_primitive_candidates = 0;

    // Primitive type name breakdown (name -> count)
    std::map<std::string, int> primitives_by_type;
};

/**
 * @brief Global profiler for clustering performance analysis
 */
class ClusterProfiler {
public:
    // Singleton access
    static ClusterProfiler& instance() {
        static ClusterProfiler profiler;
        return profiler;
    }

    // Enable/disable profiling
    void enable() { enabled_ = true; }
    void disable() { enabled_ = false; }
    bool is_enabled() const { return enabled_; }

    // Start/end cluster tracking
    void start_cluster(size_t cluster_id, const std::string& name, const std::string& type) {
        if (!enabled_) return;
        current_stats_ = ClusterStats();
        current_stats_.cluster_id = cluster_id;
        current_stats_.cluster_name = name;
        current_stats_.cluster_type = type;
        cluster_start_time_ = now();
    }

    void end_cluster() {
        if (!enabled_) return;
        current_stats_.total_time_us = elapsed_us(cluster_start_time_);
        if (current_stats_.molecules_tried > 0) {
            current_stats_.avg_primitives_per_molecule =
                static_cast<double>(current_stats_.total_primitive_candidates) / current_stats_.molecules_tried;
        }
        all_cluster_stats_.push_back(current_stats_);
    }

    // Phase timing
    void start_phase(const std::string& phase) {
        if (!enabled_) return;
        phase_start_times_[phase] = now();
    }

    void end_phase(const std::string& phase) {
        if (!enabled_) return;
        auto it = phase_start_times_.find(phase);
        if (it != phase_start_times_.end()) {
            double elapsed = elapsed_us(it->second);
            if (phase == "start_cluster") {
                current_stats_.start_cluster_time_us += elapsed;
            } else if (phase == "molecule_selection") {
                current_stats_.molecule_selection_time_us += elapsed;
            } else if (phase == "molecule_packing") {
                current_stats_.molecule_packing_time_us += elapsed;
            } else if (phase == "stats_update") {
                current_stats_.stats_update_time_us += elapsed;
            } else if (phase == "legality_check") {
                current_stats_.legality_check_time_us += elapsed;
            }
        }
    }

    // Molecule tracking
    void record_molecule_attempt(bool is_chain) {
        if (!enabled_) return;
        current_stats_.molecules_tried++;
        if (is_chain) current_stats_.chain_molecules_tried++;
    }

    void record_molecule_success(bool is_chain) {
        if (!enabled_) return;
        current_stats_.molecules_packed++;
        if (is_chain) current_stats_.chain_molecules_packed++;
    }

    void record_molecule_failure() {
        if (!enabled_) return;
        current_stats_.molecules_failed++;
    }

    // Detailed packing stats
    void record_primitive_candidate(bool is_chain) {
        if (!enabled_) return;
        current_stats_.total_primitive_candidates++;
        if (is_chain) {
            current_stats_.chain_primitive_candidates++;
        } else {
            current_stats_.non_chain_primitive_candidates++;
        }
    }

    void record_primitive_candidate_by_type(const std::string& type_name) {
        if (!enabled_) return;
        current_stats_.primitives_by_type[type_name]++;
    }

    void record_early_chain_reject() {
        if (!enabled_) return;
        current_stats_.early_chain_rejects++;
    }

    void record_placement_failure() {
        if (!enabled_) return;
        current_stats_.placement_failures++;
    }

    void record_routing_failure() {
        if (!enabled_) return;
        current_stats_.routing_failures++;
    }

    void record_expansion_failure() {
        if (!enabled_) return;
        current_stats_.expansion_failures++;
    }

    // Output results
    void write_report(const char* filename) {
        FILE* fp = fopen(filename, "w");
        if (!fp) return;

        fprintf(fp, "=================================================================\n");
        fprintf(fp, "Clustering Performance Profile\n");
        fprintf(fp, "=================================================================\n\n");

        // Summary statistics
        double total_time = 0;
        double total_start_cluster = 0;
        double total_mol_selection = 0;
        double total_mol_packing = 0;
        double total_stats_update = 0;
        double total_legality = 0;
        int total_molecules_tried = 0;
        int total_molecules_packed = 0;
        int total_primitives = 0;
        int total_early_rejects = 0;
        int total_placement_fails = 0;
        int total_routing_fails = 0;

        for (const auto& stats : all_cluster_stats_) {
            total_time += stats.total_time_us;
            total_start_cluster += stats.start_cluster_time_us;
            total_mol_selection += stats.molecule_selection_time_us;
            total_mol_packing += stats.molecule_packing_time_us;
            total_stats_update += stats.stats_update_time_us;
            total_legality += stats.legality_check_time_us;
            total_molecules_tried += stats.molecules_tried;
            total_molecules_packed += stats.molecules_packed;
            total_primitives += stats.total_primitive_candidates;
            total_early_rejects += stats.early_chain_rejects;
            total_placement_fails += stats.placement_failures;
            total_routing_fails += stats.routing_failures;
        }

        // Also compute chain vs non-chain breakdown
        int total_chain_prims = 0;
        int total_non_chain_prims = 0;
        int total_expansion_fails = 0;
        std::map<std::string, long long> total_by_type;
        for (const auto& stats : all_cluster_stats_) {
            total_chain_prims += stats.chain_primitive_candidates;
            total_non_chain_prims += stats.non_chain_primitive_candidates;
            total_expansion_fails += stats.expansion_failures;
            for (const auto& kv : stats.primitives_by_type) {
                total_by_type[kv.first] += kv.second;
            }
        }

        fprintf(fp, "SUMMARY\n");
        fprintf(fp, "-------\n");
        fprintf(fp, "Total clusters created: %zu\n", all_cluster_stats_.size());
        fprintf(fp, "Total time: %.2f ms\n", total_time / 1000.0);
        fprintf(fp, "\n");

        fprintf(fp, "Time Breakdown:\n");
        fprintf(fp, "  Start cluster:      %8.2f ms (%5.1f%%)\n",
                total_start_cluster / 1000.0, 100.0 * total_start_cluster / total_time);
        fprintf(fp, "  Molecule selection: %8.2f ms (%5.1f%%)\n",
                total_mol_selection / 1000.0, 100.0 * total_mol_selection / total_time);
        fprintf(fp, "  Molecule packing:   %8.2f ms (%5.1f%%)\n",
                total_mol_packing / 1000.0, 100.0 * total_mol_packing / total_time);
        fprintf(fp, "  Stats update:       %8.2f ms (%5.1f%%)\n",
                total_stats_update / 1000.0, 100.0 * total_stats_update / total_time);
        fprintf(fp, "  Legality check:     %8.2f ms (%5.1f%%)\n",
                total_legality / 1000.0, 100.0 * total_legality / total_time);
        fprintf(fp, "\n");

        fprintf(fp, "Molecule Statistics:\n");
        fprintf(fp, "  Total molecules tried:  %d\n", total_molecules_tried);
        fprintf(fp, "  Total molecules packed: %d\n", total_molecules_packed);
        fprintf(fp, "  Success rate: %.1f%%\n",
                total_molecules_tried > 0 ? 100.0 * total_molecules_packed / total_molecules_tried : 0);
        fprintf(fp, "\n");

        fprintf(fp, "Placement Statistics:\n");
        fprintf(fp, "  Total primitive candidates tried: %d\n", total_primitives);
        fprintf(fp, "    - For chain molecules:     %d (%.1f%%)\n",
                total_chain_prims, total_primitives > 0 ? 100.0 * total_chain_prims / total_primitives : 0);
        fprintf(fp, "    - For non-chain molecules: %d (%.1f%%)\n",
                total_non_chain_prims, total_primitives > 0 ? 100.0 * total_non_chain_prims / total_primitives : 0);
        fprintf(fp, "  Early chain rejects: %d (%.1f%% of chain candidates)\n",
                total_early_rejects,
                total_chain_prims > 0 ? 100.0 * total_early_rejects / total_chain_prims : 0);
        fprintf(fp, "  Expansion failures: %d\n", total_expansion_fails);
        fprintf(fp, "  Placement failures: %d\n", total_placement_fails);
        fprintf(fp, "  Routing failures: %d\n", total_routing_fails);
        if (total_molecules_tried > 0) {
            fprintf(fp, "  Avg primitives per molecule: %.1f\n",
                    static_cast<double>(total_primitives) / total_molecules_tried);
        }
        fprintf(fp, "\n");

        // Primitives by type breakdown
        if (!total_by_type.empty()) {
            fprintf(fp, "Primitive Candidates by Type:\n");
            // Sort by count descending
            std::vector<std::pair<std::string, long long>> sorted_types(
                total_by_type.begin(), total_by_type.end());
            std::sort(sorted_types.begin(), sorted_types.end(),
                      [](const auto& a, const auto& b) {
                          return a.second > b.second;
                      });
            for (const auto& kv : sorted_types) {
                fprintf(fp, "  %-25s %12lld (%5.1f%%)\n",
                        kv.first.c_str(), kv.second,
                        total_primitives > 0 ? 100.0 * kv.second / total_primitives : 0);
            }
            fprintf(fp, "\n");
        }

        // Per-cluster details for slow clusters
        fprintf(fp, "=================================================================\n");
        fprintf(fp, "PER-CLUSTER DETAILS (sorted by time, top 20)\n");
        fprintf(fp, "=================================================================\n\n");

        // Sort by time (descending)
        std::vector<const ClusterStats*> sorted_stats;
        for (const auto& stats : all_cluster_stats_) {
            sorted_stats.push_back(&stats);
        }
        std::sort(sorted_stats.begin(), sorted_stats.end(),
                  [](const ClusterStats* a, const ClusterStats* b) {
                      return a->total_time_us > b->total_time_us;
                  });

        int count = 0;
        for (const auto* stats : sorted_stats) {
            if (count++ >= 20) break;

            fprintf(fp, "Cluster %zu: %s (type: %s)\n",
                    stats->cluster_id, stats->cluster_name.c_str(), stats->cluster_type.c_str());
            fprintf(fp, "  Total time: %.2f ms\n", stats->total_time_us / 1000.0);
            fprintf(fp, "  Breakdown: start=%.2fms sel=%.2fms pack=%.2fms stats=%.2fms legal=%.2fms\n",
                    stats->start_cluster_time_us / 1000.0,
                    stats->molecule_selection_time_us / 1000.0,
                    stats->molecule_packing_time_us / 1000.0,
                    stats->stats_update_time_us / 1000.0,
                    stats->legality_check_time_us / 1000.0);
            fprintf(fp, "  Molecules: tried=%d packed=%d failed=%d\n",
                    stats->molecules_tried, stats->molecules_packed, stats->molecules_failed);
            fprintf(fp, "  Chains: tried=%d packed=%d\n",
                    stats->chain_molecules_tried, stats->chain_molecules_packed);
            fprintf(fp, "  Primitives: total=%d (chain=%d non-chain=%d)\n",
                    stats->total_primitive_candidates,
                    stats->chain_primitive_candidates, stats->non_chain_primitive_candidates);
            fprintf(fp, "  Failures: early_rej=%d expand=%d place=%d route=%d\n",
                    stats->early_chain_rejects, stats->expansion_failures,
                    stats->placement_failures, stats->routing_failures);
            fprintf(fp, "\n");
        }

        // All clusters summary table
        fprintf(fp, "=================================================================\n");
        fprintf(fp, "ALL CLUSTERS (time in ms)\n");
        fprintf(fp, "=================================================================\n");
        fprintf(fp, "%-6s %-30s %-15s %8s %6s %6s %6s %8s\n",
                "ID", "Name", "Type", "Time", "Tried", "Packed", "Prims", "EarlyRej");
        fprintf(fp, "%-6s %-30s %-15s %8s %6s %6s %6s %8s\n",
                "------", "------------------------------", "---------------",
                "--------", "------", "------", "------", "--------");

        for (const auto& stats : all_cluster_stats_) {
            std::string short_name = stats.cluster_name;
            if (short_name.length() > 30) {
                short_name = short_name.substr(0, 27) + "...";
            }
            std::string short_type = stats.cluster_type;
            if (short_type.length() > 15) {
                short_type = short_type.substr(0, 12) + "...";
            }
            fprintf(fp, "%-6zu %-30s %-15s %8.2f %6d %6d %6d %8d\n",
                    stats.cluster_id, short_name.c_str(), short_type.c_str(),
                    stats.total_time_us / 1000.0,
                    stats.molecules_tried, stats.molecules_packed,
                    stats.total_primitive_candidates, stats.early_chain_rejects);
        }

        fclose(fp);
    }

    void reset() {
        all_cluster_stats_.clear();
        current_stats_ = ClusterStats();
        phase_start_times_.clear();
    }

private:
    ClusterProfiler() = default;

    using TimePoint = std::chrono::high_resolution_clock::time_point;

    TimePoint now() {
        return std::chrono::high_resolution_clock::now();
    }

    double elapsed_us(TimePoint start) {
        auto end = now();
        return std::chrono::duration<double, std::micro>(end - start).count();
    }

    bool enabled_ = false;
    ClusterStats current_stats_;
    TimePoint cluster_start_time_;
    std::map<std::string, TimePoint> phase_start_times_;
    std::vector<ClusterStats> all_cluster_stats_;
};

// Convenience macros for profiling
#define CLUSTER_PROFILE_ENABLED() ClusterProfiler::instance().is_enabled()

#define CLUSTER_PROFILE_START_CLUSTER(id, name, type) \
    ClusterProfiler::instance().start_cluster(id, name, type)

#define CLUSTER_PROFILE_END_CLUSTER() \
    ClusterProfiler::instance().end_cluster()

#define CLUSTER_PROFILE_START_PHASE(phase) \
    ClusterProfiler::instance().start_phase(phase)

#define CLUSTER_PROFILE_END_PHASE(phase) \
    ClusterProfiler::instance().end_phase(phase)

#define CLUSTER_PROFILE_MOLECULE_ATTEMPT(is_chain) \
    ClusterProfiler::instance().record_molecule_attempt(is_chain)

#define CLUSTER_PROFILE_MOLECULE_SUCCESS(is_chain) \
    ClusterProfiler::instance().record_molecule_success(is_chain)

#define CLUSTER_PROFILE_MOLECULE_FAILURE() \
    ClusterProfiler::instance().record_molecule_failure()

#define CLUSTER_PROFILE_PRIMITIVE_CANDIDATE(is_chain) \
    ClusterProfiler::instance().record_primitive_candidate(is_chain)

#define CLUSTER_PROFILE_PRIMITIVE_CANDIDATE_BY_TYPE(type_name) \
    ClusterProfiler::instance().record_primitive_candidate_by_type(type_name)

#define CLUSTER_PROFILE_EARLY_CHAIN_REJECT() \
    ClusterProfiler::instance().record_early_chain_reject()

#define CLUSTER_PROFILE_PLACEMENT_FAILURE() \
    ClusterProfiler::instance().record_placement_failure()

#define CLUSTER_PROFILE_ROUTING_FAILURE() \
    ClusterProfiler::instance().record_routing_failure()

#define CLUSTER_PROFILE_EXPANSION_FAILURE() \
    ClusterProfiler::instance().record_expansion_failure()
