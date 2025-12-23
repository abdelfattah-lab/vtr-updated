/*
* Added by Junius Pun (juniuspun00@gmail.com)
* 
* Implementation of compressor trees to shrink multi-level additions into Boolean logic that can be packed into LUTs.
*/

#include "compressor.h"
#include "netlist_utils.h"
#include "node_utils.h"
#include "odin_util.h"

#include "adder.h"
#include "vtr_list.h"
#include "log.h"

#include <vector>
#include <tuple>

using vtr::insert_in_vptr_list;

// helper functions.
static signal_list_t *ranks_to_adder_chain(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> &ranks);
static npin_t *make_output_pin(nnode_t *node, int idx);

// implement basic gates.
static npin_t *implement_AND(nnode_t *node, short mark, npin_t *a, npin_t *b);

// make full adder and half adder equivalent boolean gates.
static std::pair<npin_t *, npin_t *> implement_FA(nnode_t *node, short mark, npin_t *a, npin_t *b, npin_t *c);
static std::pair<npin_t *, npin_t *> implement_HA(nnode_t *node, short mark, npin_t *a, npin_t *b);

// compressor tree implementations.
static signal_list_t *implement_compressor_tree_wallace(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);
static signal_list_t *implement_compressor_tree_wallace_ternary(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);
static signal_list_t *implement_compressor_tree_dadda(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);
static signal_list_t *implement_compressor_tree_cascade(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);
static signal_list_t *implement_compressor_tree_ternary(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);

// helper function for ternary final stage - handles up to 3 rows with chained adders.
static signal_list_t *ranks_to_ternary_adder_chain(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> &ranks);

/*-----------------------------------------------------------
* (function: implement_compressor_tree)
* 
* @brief compresses a given multi-level addition, arranged by rank, into a single row of output pins.
*
* @note this uses a compressor tree approach as provided.
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
signal_list_t *implement_compressor_tree(compressor_tree_type_e tree_type, nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks)
{
    switch (tree_type) {
        case compressor_tree_type_e::WALLACE:
            // use Wallace tree.
            return implement_compressor_tree_wallace(node, mark, netlist, ranks);
        case compressor_tree_type_e::WALLACE_TERNARY:
            // use Wallace tree with ternary final stage.
            return implement_compressor_tree_wallace_ternary(node, mark, netlist, ranks);
        case compressor_tree_type_e::WALLACE_TERNARY_EXP:
            // experimental version now uses the same implementation as WALLACE_TERNARY.
            return implement_compressor_tree_wallace_ternary(node, mark, netlist, ranks);
        case compressor_tree_type_e::DADDA:
            // use Dadda tree.
            return implement_compressor_tree_dadda(node, mark, netlist, ranks);
        case compressor_tree_type_e::CASCADE:
            // use cascade-friendly sequential accumulation.
            return implement_compressor_tree_cascade(node, mark, netlist, ranks);
        case compressor_tree_type_e::TERNARY_TREE:
            // use ternary adder tree for DCC3 chain topology.
            return implement_compressor_tree_ternary(node, mark, netlist, ranks);
        default:
            // invalid type; throw an error.
            Yosys::log_error("Unrecognized compressor tree type. Please use one of the following in compressor_tree_type_e enum, defined under 'parmys-plugin/core/compressor.h'.");
    }
}

/*-----------------------------------------------------------
* (function: implement_compressor_tree_wallace)
* 
* @brief compresses a given multi-level addition, arranged by rank, into a single row of output pins.
*
* @note this uses Asif & Kong's Proposed Wallace tree approach (https://doi.org/10.1155/2014/343960).
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
static signal_list_t *implement_compressor_tree_wallace(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks) 
{
    std::vector<std::vector<npin_t *>> temp;
    int i, rank_size, max_rank_size = 0, cur_ranks_size = ranks.size();

    // get the maximum rank size.
    for (i = 0; i < cur_ranks_size; i++) {
        rank_size = ranks[i].size();
        if (rank_size > max_rank_size) max_rank_size = rank_size;
    }

    // reduce ranks until only 2 rows remain.
    while (max_rank_size > 2) {
        // initialize tracker variables for this scope.
        int target_rank_size = (max_rank_size / 3) * 2 + (max_rank_size % 3);
        bool is_first_reducible_rank = true;
        int last_adder_count = 0;
        int new_ranks_size = 0;

        // reduce current ranks.
        for (i = 0; i < cur_ranks_size; i++) {
            int cur_adder_count = 0;

            // get rank and size.
            rank_size = ranks[i].size();

            // make new row for this rank.
            if (new_ranks_size < i + 1) {
                std::vector<npin_t *> r0;
                temp.push_back(r0);
                new_ranks_size++;
            }
            if (rank_size < 2) {
                // skip if there is no need to reduce.
                last_adder_count = 0;
                continue;
            }
            // make new row for generated carries.
            if (new_ranks_size < i + 2) {
                std::vector<npin_t *> r1;
                temp.push_back(r1);
                new_ranks_size++;
            }

            // make as many FAs as possible.
            while (rank_size >= 3) {
                // make FA with last 3 pins.
                npin_t *sum, *carry, *a, *b, *c;
                a = ranks[i].back();
                ranks[i].pop_back();
                b = ranks[i].back();
                ranks[i].pop_back();
                c = ranks[i].back();
                ranks[i].pop_back();
                std::tie(sum, carry) = implement_FA(node, mark, a, b, c);

                // add FA output pins to new ranks.
                temp[i].push_back(sum);
                temp[i+1].push_back(carry);

                // reduce size.
                rank_size -= 3;

                // add to current adder count, and mark first rank flag as invalid.
                cur_adder_count++;
                is_first_reducible_rank = false;
            }

            // insert HA only if (a) target rows need to be met, or (b) first rank with size >= 2.
            if (rank_size == 2 // check for HA eligibility (if rank_size > 2, then an FA would have been inserted.)
                && (
                    is_first_reducible_rank // is first reducible rank.
                    || cur_adder_count + last_adder_count + rank_size > target_rank_size // requires reduction to target size.
                )
            ) {
                // make HA with last 2 pins.
                npin_t *sum, *carry, *a, *b;
                a = ranks[i].back();
                ranks[i].pop_back();
                b = ranks[i].back();
                ranks[i].pop_back();
                std::tie(sum, carry) = implement_HA(node, mark, a, b);

                // add HA output pins to new ranks.
                temp[i].push_back(sum);
                temp[i+1].push_back(carry);

                // reduce size.
                rank_size -= 2;

                // add to current adder count, and mark first rank flag as invalid.
                cur_adder_count++;
                is_first_reducible_rank = false;
            }
            // assign last_adder_count for next rank.
            last_adder_count = cur_adder_count;
        }


        max_rank_size = 0;
        // re-append all from temp.
        for (i = 0; i < new_ranks_size; i++) {
            std::vector<npin_t *> temp_rank = temp[i];
            if (temp_rank.size()) {
                // copy all elements over to corresponding rank, else add new rank.
                if (i < cur_ranks_size) {
                    ranks[i].insert(ranks[i].end(), temp_rank.begin(), temp_rank.end());
                }
                else {
                    ranks.push_back(temp_rank);
                }
            }
            else if (i >= cur_ranks_size) {
                // rank will be invalid (carry vector created but no pins added.)
                continue;
            }
            
            // check for max size.
            rank_size = ranks[i].size();
            if (rank_size > max_rank_size) {
                max_rank_size = rank_size;
            }
        }

        // clear temp vector.
        temp.clear();
        
        // re-assign rank size.
        cur_ranks_size = ranks.size();
    }

    // return final rows combined with adder chain.
    return ranks_to_adder_chain(node, mark, netlist, ranks);
}

/*-----------------------------------------------------------
* (function: implement_compressor_tree_wallace_ternary)
*
* @brief compresses a given multi-level addition, arranged by rank, into a single row of output pins.
*
* @note This is a modified Wallace tree that reduces to max height 3 instead of 2,
*       then uses a ternary adder chain (A+B)+C for the final stage.
*       This is optimized for architectures with ternary adder support (e.g., DCC3).
*
*       Reduction strategy with carry tracking:
*       - Track carries from previous columns to determine effective height
*       - Use FA when effective_height >= 4 (more aggressive than before)
*       - Use HA when effective_height == 4 and only 2 pins remain after FAs
*       - Target: reduce each column to height <= 3
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
static signal_list_t *implement_compressor_tree_wallace_ternary(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks)
{
    std::vector<std::vector<npin_t *>> temp;
    int i, rank_size, max_rank_size = 0, cur_ranks_size = ranks.size();

    // get the maximum rank size.
    for (i = 0; i < cur_ranks_size; i++) {
        rank_size = ranks[i].size();
        if (rank_size > max_rank_size) max_rank_size = rank_size;
    }

    // reduce ranks until only 3 rows remain (instead of 2 for standard Wallace).
    // Key improvement: track carries between columns for better reduction decisions.
    while (max_rank_size > 3) {
        int new_ranks_size = 0;
        int last_carry_count = 0;  // carries generated by previous column

        // reduce current ranks.
        for (i = 0; i < cur_ranks_size; i++) {
            int cur_carry_count = 0;  // carries we'll generate for next column

            // get rank and size.
            rank_size = ranks[i].size();

            // Effective height includes incoming carries from previous column
            int effective_height = rank_size + last_carry_count;

            // make new row for this rank.
            if (new_ranks_size < i + 1) {
                std::vector<npin_t *> r0;
                temp.push_back(r0);
                new_ranks_size++;
            }
            if (effective_height < 2) {
                // skip if there is no need to reduce.
                last_carry_count = 0;
                continue;
            }
            // make new row for generated carries.
            if (new_ranks_size < i + 2) {
                std::vector<npin_t *> r1;
                temp.push_back(r1);
                new_ranks_size++;
            }

            // Improved reduction strategy for ternary target (height 3):
            // - FA: consumes 3 pins, produces 1 sum (stays) + 1 carry (goes to next column)
            //   Net effect on this column: -2 height
            // - HA: consumes 2 pins, produces 1 sum (stays) + 1 carry (goes to next column)
            //   Net effect on this column: -1 height
            //
            // To reach height 3 from height h:
            //   h=4: need to reduce by 1 → 1 HA, or 1 FA if we have 3+ pins (leaves 1 pin + incoming carries)
            //   h=5: need to reduce by 2 → 1 FA
            //   h=6: need to reduce by 3 → 1 FA + 1 HA, or 2 FAs if we have 6+ pins
            //   h=7: need to reduce by 4 → 2 FAs
            //   etc.
            //
            // Strategy: Use FAs as much as possible (they're more efficient), then HA if needed.

            // Calculate target: we want final height <= 3
            // After reduction, height = rank_size - 2*num_FA - num_HA + incoming_carries_from_temp
            // But incoming_carries_from_temp is what we're generating now, so it's complex.
            // Simpler: just reduce aggressively while height > 3.

            // Use FAs while we have 3+ pins AND effective height > 3
            while (rank_size >= 3 && (rank_size + cur_carry_count) > 3) {
                // make FA with last 3 pins.
                npin_t *sum, *carry, *a, *b, *c;
                a = ranks[i].back();
                ranks[i].pop_back();
                b = ranks[i].back();
                ranks[i].pop_back();
                c = ranks[i].back();
                ranks[i].pop_back();
                std::tie(sum, carry) = implement_FA(node, mark, a, b, c);

                // add FA output pins to new ranks.
                temp[i].push_back(sum);
                temp[i+1].push_back(carry);

                // reduce size and track carry.
                rank_size -= 3;
                cur_carry_count++;
            }

            // Use HA if we still have height > 3 and have 2+ pins
            // This happens when rank_size == 2 and cur_carry_count >= 2
            while (rank_size >= 2 && (rank_size + cur_carry_count) > 3) {
                // make HA with last 2 pins.
                npin_t *sum, *carry, *a, *b;
                a = ranks[i].back();
                ranks[i].pop_back();
                b = ranks[i].back();
                ranks[i].pop_back();
                std::tie(sum, carry) = implement_HA(node, mark, a, b);

                // add HA output pins to new ranks.
                temp[i].push_back(sum);
                temp[i+1].push_back(carry);

                // reduce size and track carry.
                rank_size -= 2;
                cur_carry_count++;
            }

            // Update carry count for next column
            last_carry_count = cur_carry_count;
        }


        max_rank_size = 0;
        // re-append all from temp.
        for (i = 0; i < new_ranks_size; i++) {
            std::vector<npin_t *> temp_rank = temp[i];
            if (temp_rank.size()) {
                // copy all elements over to corresponding rank, else add new rank.
                if (i < cur_ranks_size) {
                    ranks[i].insert(ranks[i].end(), temp_rank.begin(), temp_rank.end());
                }
                else {
                    ranks.push_back(temp_rank);
                }
            }
            else if (i >= cur_ranks_size) {
                // rank will be invalid (carry vector created but no pins added.)
                continue;
            }

            // check for max size.
            rank_size = ranks[i].size();
            if (rank_size > max_rank_size) {
                max_rank_size = rank_size;
            }
        }

        // clear temp vector.
        temp.clear();

        // re-assign rank size.
        cur_ranks_size = ranks.size();
    }

    // return final rows combined with ternary adder chain (handles up to 3 rows).
    return ranks_to_ternary_adder_chain(node, mark, netlist, ranks);
}

/*-----------------------------------------------------------
* (function: ranks_to_ternary_adder_chain)
*
* @brief Converts ranks with height <= 3 to a ternary adder chain.
*
* @note For ternary adder architectures (e.g., DCC3), this creates (A+B)+C pattern:
*       - Height 1: Direct connection (no adder needed)
*       - Height 2: Binary adder (A+B)
*       - Height 3: Ternary chain - first adder computes (A+B),
*                   its sumout feeds second adder which computes (A+B)+C
*
*       Chain order optimization for DCC3:
*       The key insight is that sumout of adder1 directly feeds input of adder2.
*       We want to minimize:
*       1. Total adder bit-width (W1 + W2)
*       2. "Wasted" bits - zero-padded positions that don't carry useful data
*       3. Misalignment between sumout and C's input span
*
*       Strategy: Try all 3 orderings, compute a cost metric, pick the best.
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
static signal_list_t *ranks_to_ternary_adder_chain(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> &ranks)
{
    int cur_ranks_size = ranks.size();

    // Find maximum rank height.
    int max_rank_size = 0;
    for (int i = 0; i < cur_ranks_size; i++) {
        int rank_size = ranks[i].size();
        if (rank_size > max_rank_size) max_rank_size = rank_size;
    }

    // If max height <= 2, fall back to standard binary adder chain.
    if (max_rank_size <= 2) {
        return ranks_to_adder_chain(node, mark, netlist, ranks);
    }

    // max_rank_size == 3: Use ternary adder chain.
    // Transpose ranks into rows: row_i[bit] = ranks[bit][i]
    std::vector<std::vector<npin_t *>> rows;
    for (int row_idx = 0; row_idx < max_rank_size; row_idx++) {
        std::vector<npin_t *> row;
        for (int bit_pos = 0; bit_pos < cur_ranks_size; bit_pos++) {
            if (row_idx < (int)ranks[bit_pos].size()) {
                row.push_back(ranks[bit_pos][row_idx]);
            } else {
                row.push_back(get_zero_pin(netlist));
            }
        }
        rows.push_back(row);
    }

    // Assertion: We should have at most 3 rows at this point
    oassert(max_rank_size <= 3 && "WALLACE_TERNARY should reduce to at most 3 rows");

    // Helper: Compute the actual bit span of a row (first and last non-zero positions)
    // Returns {first_nonzero, last_nonzero} or {-1, -1} if all zeros
    nnet_t *gnd_net = netlist->zero_net;
    auto get_span = [&gnd_net](const std::vector<npin_t *> &row) -> std::pair<int, int> {
        int first = -1, last = -1;
        for (int i = 0; i < (int)row.size(); i++) {
            if (row[i] && row[i]->net && strcmp(row[i]->net->name, gnd_net->name) != 0) {
                if (first == -1) first = i;
                last = i;
            }
        }
        return {first, last};
    };

    // Get spans for all 3 rows
    std::pair<int, int> spans[3];
    for (int r = 0; r < 3; r++) {
        spans[r] = get_span(rows[r]);
        // If a row is all zeros, treat it as spanning [0, 0] for calculation purposes
        if (spans[r].first == -1) {
            spans[r] = {0, 0};
        }
    }

    // Helper: Calculate cost for a given chain order (a, b, c) meaning (rows[a] + rows[b]) + rows[c]
    // Cost considers:
    // 1. Total adder width (W1 + W2) - primary metric
    // 2. Sumout-to-C alignment bonus - prefer C that aligns well with sumout
    // 3. Chain utilization - prefer orderings where sumout bits are "useful"
    auto calc_chain_cost = [&](int a, int b, int c) -> int {
        // First adder: rows[a] + rows[b]
        // Span of first adder output starts at min of the two row starts
        int start_ab = std::min(spans[a].first, spans[b].first);
        int end_ab = std::max(spans[a].second, spans[b].second);
        // Width includes +1 for potential carry
        int width_ab = (end_ab - start_ab + 1) + 1;

        // Second adder: sumout_ab + rows[c]
        // sumout_ab has span [start_ab, start_ab + width_ab - 1]
        int sumout_start = start_ab;
        int sumout_end = start_ab + width_ab - 1;

        int start_abc = std::min(sumout_start, spans[c].first);
        int end_abc = std::max(sumout_end, spans[c].second);
        int width_abc = (end_abc - start_abc + 1) + 1;

        // Primary cost: total adder width
        int total_width = width_ab + width_abc;

        // Secondary: penalize misalignment between sumout and C
        // If C's span doesn't overlap well with sumout, we're wasting bits
        int overlap_start = std::max(sumout_start, spans[c].first);
        int overlap_end = std::min(sumout_end, spans[c].second);
        int overlap = std::max(0, overlap_end - overlap_start + 1);

        // Calculate "wasted" bits in second adder:
        // - sumout bits that don't overlap with C
        // - C bits that don't overlap with sumout
        int sumout_len = sumout_end - sumout_start + 1;
        int c_len = spans[c].second - spans[c].first + 1;
        int wasted = (sumout_len - overlap) + (c_len - overlap);

        // Final cost: weighted sum (total_width is primary, wasted is secondary)
        // Lower is better
        return total_width * 10 + wasted;
    };

    // Try all 3 orderings and find the best
    // Order 0: (0+1)+2, Order 1: (0+2)+1, Order 2: (1+2)+0
    int orderings[3][3] = {
        {0, 1, 2},  // (row0 + row1) + row2
        {0, 2, 1},  // (row0 + row2) + row1
        {1, 2, 0}   // (row1 + row2) + row0
    };

    int best_order = 0;
    int best_cost = calc_chain_cost(orderings[0][0], orderings[0][1], orderings[0][2]);

    for (int o = 1; o < 3; o++) {
        int cost = calc_chain_cost(orderings[o][0], orderings[o][1], orderings[o][2]);
        if (cost < best_cost) {
            best_cost = cost;
            best_order = o;
        }
    }

    int first_a = orderings[best_order][0];
    int first_b = orderings[best_order][1];
    int second_c = orderings[best_order][2];

    std::vector<npin_t *> &row_a = rows[first_a];
    std::vector<npin_t *> &row_b = rows[first_b];
    std::vector<npin_t *> &row_c = rows[second_c];

    int width_a = row_a.size();
    int width_b = row_b.size();
    int width_c = row_c.size();

    // First adder: A + B
    int width_ab = std::max(width_a, width_b) + 1;
    nnode_t *add1 = make_2port_gate(ADD, width_ab, width_ab, width_ab, node, mark);
    add_list = insert_in_vptr_list(add_list, add1);

    // Connect row_a to port A of first adder.
    for (int bit = 0; bit < width_ab; bit++) {
        npin_t *pin_a = (bit < width_a) ? row_a[bit] : get_zero_pin(netlist);
        add_input_pin_to_node(add1, pin_a, bit);
    }
    // Connect row_b to port B of first adder.
    for (int bit = 0; bit < width_ab; bit++) {
        npin_t *pin_b = (bit < width_b) ? row_b[bit] : get_zero_pin(netlist);
        add_input_pin_to_node(add1, pin_b, width_ab + bit);
    }

    // Get sumout from first adder - this feeds input of second adder.
    std::vector<npin_t *> sumout_ab;
    for (int bit = 0; bit < width_ab; bit++) {
        sumout_ab.push_back(make_output_pin(add1, bit));
    }

    // Second adder: (A+B) + C
    // KEY: sumout of adder1 feeds input of adder2 -> ternary chain pattern!
    int width_abc = std::max(width_ab, width_c) + 1;
    nnode_t *add2 = make_2port_gate(ADD, width_abc, width_abc, width_abc, node, mark);
    add_list = insert_in_vptr_list(add_list, add2);

    // Connect sumout_ab to port A of second adder.
    for (int bit = 0; bit < width_abc; bit++) {
        npin_t *pin_a = (bit < width_ab) ? sumout_ab[bit] : get_zero_pin(netlist);
        add_input_pin_to_node(add2, pin_a, bit);
    }
    // Connect row_c to port B of second adder.
    for (int bit = 0; bit < width_abc; bit++) {
        npin_t *pin_b = (bit < width_c) ? row_c[bit] : get_zero_pin(netlist);
        add_input_pin_to_node(add2, pin_b, width_abc + bit);
    }

    // Build output signal list from second adder's output.
    signal_list_t *ret = init_signal_list();
    for (int bit = 0; bit < width_abc; bit++) {
        add_pin_to_signal_list(ret, make_output_pin(add2, bit));
    }

    return ret;
}

/*-----------------------------------------------------------
* (function: implement_compressor_tree_dadda)
*
* @brief compresses a given multi-level addition, arranged by rank, into a single row of output pins.
*
* @note this uses a Dadda tree approach.
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
static signal_list_t *implement_compressor_tree_dadda(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks) {
    std::vector<std::vector<npin_t *>> temp;
    int i, rank_size, max_rank_size = 0, cur_ranks_size = ranks.size();

    // get the maximum rank size.
    for (i = 0; i < cur_ranks_size; i++) {
        rank_size = ranks[i].size();
        if (rank_size > max_rank_size) max_rank_size = rank_size;
    }

    // get all d-factors.
    std::vector<int> d_factors;
    int d = 2;
    while (d < max_rank_size) {
        d_factors.push_back(d);
        d = d * 3 / 2;
    }

    if (!d_factors.empty()) {
        // roll d back by one.
        d = d_factors.back();
        d_factors.pop_back();
    }

    /* reduce ranks according to Dadda's algorithm:
        0. Define rank_size' as rank_size + last_carry_count + cur_adder_count.
        1. if rank_size' <= d, then move to next rank.
        2. if rank_size' = d+1, then combine 2 elements with HA, cur_adder_count += 1, then move to next rank.
        3. else combine 3 elements with FA, cur_adder_count += 1, then repeat from 1.
        4. assign last_carry_count = cur_adder_count.
        5. repeat steps 1-4 until all ranks have <= 2 elements.
    */
    while (max_rank_size > 2) {
        int new_ranks_size = 0;
        
        int last_carry_count = 0;
        // reduce current ranks.
        for (i = 0; i < cur_ranks_size; i++) {
            int cur_adder_count = 0;
            // get rank and size.
            rank_size = ranks[i].size();

            // make new row for this rank.
            if (new_ranks_size < i + 1) {
                std::vector<npin_t *> r0;
                temp.push_back(r0);
                new_ranks_size++;
            }
            if (rank_size + last_carry_count <= d) {
                // skip if there is no need to reduce.
                continue;
            }
            // make new row for generated carries.
            if (new_ranks_size < i + 2) {
                std::vector<npin_t *> r1;
                temp.push_back(r1);
                new_ranks_size++;
            }

            // make as many FAs as possible.
            while (rank_size + last_carry_count + cur_adder_count > d+1 && rank_size >= 3) {
                // make FA with last 3 pins.
                npin_t *sum, *carry, *a, *b, *c;
                a = ranks[i].back();
                ranks[i].pop_back();
                b = ranks[i].back();
                ranks[i].pop_back();
                c = ranks[i].back();
                ranks[i].pop_back();
                std::tie(sum, carry) = implement_FA(node, mark, a, b, c);

                // add FA output pins to new ranks.
                temp[i].push_back(sum);
                temp[i+1].push_back(carry);

                // reduce size.
                rank_size -= 3;

                // add to carry count.
                cur_adder_count++;
            }

            // insert HA if rank_size = d+1.
            if (rank_size + last_carry_count + cur_adder_count == d+1 && rank_size >= 2) {
                // make FA with last 3 pins.
                npin_t *sum, *carry, *a, *b;
                a = ranks[i].back();
                ranks[i].pop_back();
                b = ranks[i].back();
                ranks[i].pop_back();
                std::tie(sum, carry) = implement_HA(node, mark, a, b);

                // add FA output pins to new ranks.
                temp[i].push_back(sum);
                temp[i+1].push_back(carry);

                // reduce size.
                rank_size -= 2;
                
                // add to carry count.
                cur_adder_count++;
            }

            // set carry count for next rank.
            last_carry_count = cur_adder_count;
        }

        max_rank_size = 0;
        // re-append all from temp.
        for (i = 0; i < new_ranks_size; i++) {
            std::vector<npin_t *> temp_rank = temp[i];
            if (temp_rank.size()) {
                // copy all elements over to corresponding rank, else add new rank.
                if (i < cur_ranks_size) {
                    ranks[i].insert(ranks[i].end(), temp_rank.begin(), temp_rank.end());
                }
                else {
                    ranks.push_back(temp_rank);
                }
            }
            
            // check for max size.
            rank_size = ranks[i].size();
            if (rank_size > max_rank_size) {
                max_rank_size = rank_size;
            }
        }

        // clear temp vector.
        temp.clear();
        
        // re-assign rank size.
        cur_ranks_size = ranks.size();

        // reduce d if possible.
        while (!d_factors.empty()) {
            d = d_factors.back();
            d_factors.pop_back();
            if (d < max_rank_size) break;
        }
    }

    // return final rows combined with adder chain.
    return ranks_to_adder_chain(node, mark, netlist, ranks);
}

/*-----------------------------------------------------------
* (function: implement_compressor_tree_cascade)
*
* @brief compresses a given multi-level addition using cascade-friendly sequential accumulation.
*
* @note This approach creates chains of adders where the sumout of one adder feeds the input
*       of the next adder, enabling efficient mapping to double-carry-chain architectures.
*       Instead of parallel reduction (Wallace/Dadda), partial products at each bit position
*       are accumulated sequentially: acc = PP[0] + PP[1], then acc = acc + PP[2], etc.
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
static signal_list_t *implement_compressor_tree_cascade(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks)
{
    int cur_ranks_size = ranks.size();
    if (cur_ranks_size == 0) {
        return init_signal_list();
    }

    // Find the maximum number of partial products at any bit position
    int max_rank_size = 0;
    for (int i = 0; i < cur_ranks_size; i++) {
        int rank_size = ranks[i].size();
        if (rank_size > max_rank_size) max_rank_size = rank_size;
    }

    // If max rank size is <= 2, we can directly use the final adder chain
    if (max_rank_size <= 2) {
        return ranks_to_adder_chain(node, mark, netlist, ranks);
    }

    // Cascade approach: accumulate partial products sequentially
    // For each pair of partial product rows, create an adder chain where:
    // - First adder: row[0] + row[1]
    // - Second adder: sumout_of_first + row[2]
    // - Third adder: sumout_of_second + row[3]
    // - etc.

    // Transpose: collect partial products by row index (instead of by bit position)
    // row_i contains the i-th partial product from each bit position
    std::vector<std::vector<npin_t *>> rows;
    for (int row_idx = 0; row_idx < max_rank_size; row_idx++) {
        std::vector<npin_t *> row;
        for (int bit_pos = 0; bit_pos < cur_ranks_size; bit_pos++) {
            if (row_idx < (int)ranks[bit_pos].size()) {
                row.push_back(ranks[bit_pos][row_idx]);
            } else {
                row.push_back(get_zero_pin(netlist));
            }
        }
        rows.push_back(row);
    }

    // Now chain the rows: acc = rows[0] + rows[1], then acc = acc + rows[2], etc.
    // The key is that sumout of each adder feeds the input of the next adder

    std::vector<npin_t *> accumulator = rows[0];
    int acc_width = accumulator.size();

    for (int row_idx = 1; row_idx < max_rank_size; row_idx++) {
        std::vector<npin_t *> &current_row = rows[row_idx];
        int row_width = current_row.size();

        // Determine the output width (may grow by 1 bit for carry)
        int output_width = std::max(acc_width, row_width) + 1;

        // Create an adder: accumulator + current_row
        nnode_t *add_node = make_2port_gate(ADD, output_width, output_width, output_width, node, mark);
        add_list = insert_in_vptr_list(add_list, add_node);

        // Connect accumulator to port A (first input)
        for (int bit = 0; bit < output_width; bit++) {
            npin_t *pin_a;
            if (bit < acc_width) {
                pin_a = accumulator[bit];
            } else {
                pin_a = get_zero_pin(netlist);
            }
            add_input_pin_to_node(add_node, pin_a, bit);
        }

        // Connect current_row to port B (second input)
        for (int bit = 0; bit < output_width; bit++) {
            npin_t *pin_b;
            if (bit < row_width) {
                pin_b = current_row[bit];
            } else {
                pin_b = get_zero_pin(netlist);
            }
            add_input_pin_to_node(add_node, pin_b, output_width + bit);
        }

        // Create output pins for the adder (these become the new accumulator)
        // This is the key: the sumout pins will be used as inputs to the next adder
        accumulator.clear();
        for (int bit = 0; bit < output_width; bit++) {
            accumulator.push_back(make_output_pin(add_node, bit));
        }
        acc_width = output_width;
    }

    // Build the output signal list from the final accumulator
    signal_list_t *ret = init_signal_list();
    for (int bit = 0; bit < acc_width; bit++) {
        add_pin_to_signal_list(ret, accumulator[bit]);
    }

    return ret;
}

/*-----------------------------------------------------------
* (function: implement_compressor_tree_ternary)
*
* @brief compresses a given multi-level addition using a ternary adder tree.
*
* @note This approach creates a tree of ternary reduction units, where each unit
*       combines 3 inputs using 2 chained adders: (A + B) + C.
*       The sumout of the first adder feeds the input of the second adder,
*       creating the chain pattern that maps efficiently to DCC3 architecture.
*
*       Tree structure for 9 rows:
*       Level 0: [A,B,C] [D,E,F] [G,H,I]  -> 3 ternary units
*       Level 1: [R0, R1, R2]             -> 1 ternary unit
*       Level 2: [Final]
*
*       Each ternary unit: adder1(A,B) -> sumout -> adder2(sumout,C)
*
* @param ranks the i-th vector (0-indexed) contains the pins with weight 2^i.
* @returns output signal list.
* ---------------------------------------------------------*/
static signal_list_t *implement_compressor_tree_ternary(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks)
{
    int cur_ranks_size = ranks.size();
    if (cur_ranks_size == 0) {
        return init_signal_list();
    }

    // Find the maximum number of partial products at any bit position
    int max_rank_size = 0;
    for (int i = 0; i < cur_ranks_size; i++) {
        int rank_size = ranks[i].size();
        if (rank_size > max_rank_size) max_rank_size = rank_size;
    }

    // If max rank size is <= 2, we can directly use the final adder chain
    if (max_rank_size <= 2) {
        return ranks_to_adder_chain(node, mark, netlist, ranks);
    }

    // Transpose: collect partial products by row index (instead of by bit position)
    // row_i contains the i-th partial product from each bit position
    std::vector<std::vector<npin_t *>> rows;
    for (int row_idx = 0; row_idx < max_rank_size; row_idx++) {
        std::vector<npin_t *> row;
        for (int bit_pos = 0; bit_pos < cur_ranks_size; bit_pos++) {
            if (row_idx < (int)ranks[bit_pos].size()) {
                row.push_back(ranks[bit_pos][row_idx]);
            } else {
                row.push_back(get_zero_pin(netlist));
            }
        }
        rows.push_back(row);
    }

    // Reduce in ternary groups until only 1 row remains
    while (rows.size() > 1) {
        std::vector<std::vector<npin_t *>> new_rows;
        int num_rows = rows.size();

        for (int i = 0; i < num_rows; i += 3) {
            if (i + 2 < num_rows) {
                // Full ternary group: (A + B) + C
                // This creates the chain pattern: adder1.sumout -> adder2.input
                std::vector<npin_t *> &row_a = rows[i];
                std::vector<npin_t *> &row_b = rows[i + 1];
                std::vector<npin_t *> &row_c = rows[i + 2];

                int width_a = row_a.size();
                int width_b = row_b.size();
                int width_c = row_c.size();

                // First adder: A + B
                int width_ab = std::max(width_a, width_b) + 1;
                nnode_t *add1 = make_2port_gate(ADD, width_ab, width_ab, width_ab, node, mark);
                add_list = insert_in_vptr_list(add_list, add1);

                // Connect row_a to port A
                for (int bit = 0; bit < width_ab; bit++) {
                    npin_t *pin_a = (bit < width_a) ? row_a[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add1, pin_a, bit);
                }
                // Connect row_b to port B
                for (int bit = 0; bit < width_ab; bit++) {
                    npin_t *pin_b = (bit < width_b) ? row_b[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add1, pin_b, width_ab + bit);
                }

                // Get sumout from first adder
                std::vector<npin_t *> sumout_ab;
                for (int bit = 0; bit < width_ab; bit++) {
                    sumout_ab.push_back(make_output_pin(add1, bit));
                }

                // Second adder: sumout_ab + C
                // This is the KEY: sumout of adder1 feeds input of adder2 -> chain pattern!
                int width_abc = std::max(width_ab, width_c) + 1;
                nnode_t *add2 = make_2port_gate(ADD, width_abc, width_abc, width_abc, node, mark);
                add_list = insert_in_vptr_list(add_list, add2);

                // Connect sumout_ab to port A
                for (int bit = 0; bit < width_abc; bit++) {
                    npin_t *pin_a = (bit < width_ab) ? sumout_ab[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add2, pin_a, bit);
                }
                // Connect row_c to port B
                for (int bit = 0; bit < width_abc; bit++) {
                    npin_t *pin_b = (bit < width_c) ? row_c[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add2, pin_b, width_abc + bit);
                }

                // Get final sumout
                std::vector<npin_t *> result;
                for (int bit = 0; bit < width_abc; bit++) {
                    result.push_back(make_output_pin(add2, bit));
                }
                new_rows.push_back(result);

            } else if (i + 1 < num_rows) {
                // Binary pair: A + B (only 2 rows left)
                std::vector<npin_t *> &row_a = rows[i];
                std::vector<npin_t *> &row_b = rows[i + 1];

                int width_a = row_a.size();
                int width_b = row_b.size();
                int width_ab = std::max(width_a, width_b) + 1;

                nnode_t *add_node = make_2port_gate(ADD, width_ab, width_ab, width_ab, node, mark);
                add_list = insert_in_vptr_list(add_list, add_node);

                // Connect row_a to port A
                for (int bit = 0; bit < width_ab; bit++) {
                    npin_t *pin_a = (bit < width_a) ? row_a[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add_node, pin_a, bit);
                }
                // Connect row_b to port B
                for (int bit = 0; bit < width_ab; bit++) {
                    npin_t *pin_b = (bit < width_b) ? row_b[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add_node, pin_b, width_ab + bit);
                }

                // Get sumout
                std::vector<npin_t *> result;
                for (int bit = 0; bit < width_ab; bit++) {
                    result.push_back(make_output_pin(add_node, bit));
                }
                new_rows.push_back(result);

            } else {
                // Single row: pass through unchanged
                new_rows.push_back(rows[i]);
            }
        }

        rows = new_rows;
    }

    // Build the output signal list from the final row
    signal_list_t *ret = init_signal_list();
    if (!rows.empty()) {
        for (int bit = 0; bit < (int)rows[0].size(); bit++) {
            add_pin_to_signal_list(ret, rows[0][bit]);
        }
    }

    return ret;
}

// converts ranks with height <= 2 to a final adder chain (if required).
static signal_list_t *ranks_to_adder_chain(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> &ranks) {
    // make output list and combine with adders if required.
    signal_list_t *ret = init_signal_list();
    nnode_t *add_node;
    bool makingAdderChain = false;
    int adder_start_i, adder_input_size;

    int cur_ranks_size = ranks.size();
    for (int i = 0; i < cur_ranks_size; i++) {
        int rank_size = ranks[i].size();
        
        // Make the adder chain if not yet instantiated.
        if (rank_size > 1 && !makingAdderChain) {
            makingAdderChain = true;

            // calculate required input and output size.
            adder_input_size = cur_ranks_size - i;

            // make adder.
            add_node = make_2port_gate(ADD, adder_input_size, adder_input_size, adder_input_size + 1, node, mark);
            add_list = insert_in_vptr_list(add_list, add_node);

            // set adder start index.
            adder_start_i = i;
        }

        npin_t *pin;
        if (rank_size) {
            // at least one pin to insert.

            if (makingAdderChain) {
                // insert pins into adder.
                int adder_idx = i - adder_start_i;
                // insert first pin.
                add_input_pin_to_node(add_node, ranks[i].back(), adder_idx);
                ranks[i].pop_back();
                // insert second pin (if any), else zero pin.
                npin_t *second_input;
                if (rank_size > 1) {
                    second_input = ranks[i].back();
                    ranks[i].pop_back();
                }
                else {
                    second_input = get_zero_pin(netlist);
                }
                add_input_pin_to_node(add_node, second_input, adder_input_size + adder_idx);

                // assign pin to add to signal list as adder output.
                pin = make_output_pin(add_node, adder_idx);
            }
            else {
                // this branch is only reached if the rank has only one pin, before the adder chain; add it directly.
                pin = ranks[i].back();
                ranks[i].pop_back();
            }
        }
        else {
            // no pins to insert at this rank; attach a '0'.
            pin = get_zero_pin(netlist);
        }

        // add pin to signal list.
        add_pin_to_signal_list(ret, pin);
    }

    // add the last adder output of the carry chain (if it exists).
    if (makingAdderChain) {
        add_pin_to_signal_list(ret, make_output_pin(add_node, adder_input_size));
    }

    // return signal list.
    return ret;
}

// helper function to make new output pins and net, and connect them to the node.
static npin_t *make_output_pin(nnode_t *node, int idx)
{
    // make required pins and net.
    npin_t *node_out = allocate_npin(), *ret_out = allocate_npin();
    nnet_t *node_net = allocate_nnet();

    // assign names.
    node_net->name = make_full_ref_name(NULL, NULL, NULL, node->name, idx);
    ret_out->name = node_net->name;

    // attach output pin to AND node.
    add_output_pin_to_node(node, node_out, idx);
    // hook pins to net.
    add_driver_pin_to_net(node_net, node_out);
    add_fanout_pin_to_net(node_net, ret_out);

    return ret_out;
}

/*-----------------------------------------------------------
* (function: implement_AND)
* 
* @brief implements out = ab.
*
* @param a, b inputs to the AND gate.
* @returns out.
* ---------------------------------------------------------*/
static npin_t *implement_AND(nnode_t *node, short mark, npin_t *a, npin_t *b)
{
    nnode_t *and_node = make_2port_gate(LOGICAL_AND, 1, 1, 1, node, mark);

    // tie inputs.
    add_input_pin_to_node(and_node, copy_input_npin(a), 0);
    add_input_pin_to_node(and_node, copy_input_npin(b), 1);

    // tie output.
    return make_output_pin(and_node, 0);
}

/*-----------------------------------------------------------
* (function: implement_FA)
* 
* @brief converts provided pins as input to a Full Adder into sum and carry, using boolean functions sum = a^b^c and carry = ab+bc+ac.
*
* @param a, b, c inputs to the FA.
* @returns pair of output pins as (sum, carry).
* ---------------------------------------------------------*/
static std::pair<npin_t *, npin_t *> implement_FA(nnode_t *node, short mark, npin_t *a, npin_t *b, npin_t *c)
{
    // sum node.
    nnode_t *sum_node = make_3port_gate(LOGICAL_XOR, 1, 1, 1, 1, node, mark);
    add_input_pin_to_node(sum_node, copy_input_npin(a), 0);
    add_input_pin_to_node(sum_node, copy_input_npin(b), 1);
    add_input_pin_to_node(sum_node, copy_input_npin(c), 2);
    npin_t *sum = make_output_pin(sum_node, 0);

    // carry node.    
    nnode_t *carry_node = make_3port_gate(LOGICAL_OR, 1, 1, 1, 1, node, mark);
    add_input_pin_to_node(carry_node, implement_AND(node, mark, a, b), 0);
    add_input_pin_to_node(carry_node, implement_AND(node, mark, a, c), 1);
    add_input_pin_to_node(carry_node, implement_AND(node, mark, b, c), 2);
    npin_t *carry = make_output_pin(carry_node, 0);

    // return pair.
    return { sum, carry };
}

/*-----------------------------------------------------------
* (function: implement_HA)
* 
* @brief converts provided pins as input to a Half Adder into sum and carry, using boolean functions sum = a^b and carry = a+b.
*
* @param a, b inputs to the HA.
* @returns pair of output pins as (sum, carry).
* ---------------------------------------------------------*/
static std::pair<npin_t *, npin_t *> implement_HA(nnode_t *node, short mark, npin_t *a, npin_t *b)
{
    // sum node.
    nnode_t *sum_node = make_2port_gate(LOGICAL_XOR, 1, 1, 1, node, mark);
    add_input_pin_to_node(sum_node, copy_input_npin(a), 0);
    add_input_pin_to_node(sum_node, copy_input_npin(b), 1);
    npin_t *sum = make_output_pin(sum_node, 0);

    // carry node.
    npin_t *carry = implement_AND(node, mark, a, b);

    // return pair.
    return { sum, carry };
}