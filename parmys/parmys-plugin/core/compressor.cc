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
*       This is optimized for architectures with ternary adder support.
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
    // For ternary: FA reduces height by 2 (h → h-2), HA reduces by 1 (h → h-1).
    // Strategy: use FAs only when height >= 5, use HA when height == 4.
    while (max_rank_size > 3) {
        int new_ranks_size = 0;

        // reduce current ranks.
        for (i = 0; i < cur_ranks_size; i++) {
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
                continue;
            }
            // make new row for generated carries.
            if (new_ranks_size < i + 2) {
                std::vector<npin_t *> r1;
                temp.push_back(r1);
                new_ranks_size++;
            }

            // For ternary: use FAs only when height >= 5.
            // FA: consumes 3, produces sum + carry → net height reduction of 2.
            // h=5: 1 FA → remaining=2, sums=1, final=3 ✓
            // h=6: 1 FA → remaining=3, sums=1, final=4 (next pass handles)
            while (rank_size >= 5) {
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
            }

            // For ternary: use HA when height == 4 to reduce to 3.
            // HA: consumes 2, produces sum + carry → net height reduction of 1.
            // h=4: 1 HA → remaining=2, sums=1, final=3 ✓
            // Don't use HA for height <= 3 (would reduce below target).
            if (rank_size == 4) {
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
            }
            // rank_size is now <= 3 (target height for ternary)
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
* @note For ternary adder architectures, this creates (A+B)+C pattern:
*       - Height 1: Direct connection (no adder needed)
*       - Height 2: Binary adder (A+B)
*       - Height 3: Ternary chain - first adder computes (A+B),
*                   its sumout feeds second adder which computes (A+B)+C
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

    // FIX: Chain order optimization - try all 3 orderings and pick the best
    // Count actual non-zero pins in each row to estimate "density"
    // The row with the most zeros should be added last (as C in (A+B)+C)
    // to minimize the work done in the first adder
    auto count_nonzero = [&netlist](const std::vector<npin_t *> &row) {
        int count = 0;
        nnet_t *gnd_net = netlist->zero_net;
        for (auto *pin : row) {
            if (pin && pin->net && strcmp(pin->net->name, gnd_net->name) != 0) {
                count++;
            }
        }
        return count;
    };

    int density[3] = {
        count_nonzero(rows[0]),
        count_nonzero(rows[1]),
        count_nonzero(rows[2])
    };

    // Choose chain order: put the sparsest row as C (added last)
    // This minimizes the width of the first adder
    int first_a, first_b, second_c;
    if (density[2] <= density[0] && density[2] <= density[1]) {
        // row[2] is sparsest -> (row[0] + row[1]) + row[2]
        first_a = 0; first_b = 1; second_c = 2;
    } else if (density[1] <= density[0] && density[1] <= density[2]) {
        // row[1] is sparsest -> (row[0] + row[2]) + row[1]
        first_a = 0; first_b = 2; second_c = 1;
    } else {
        // row[0] is sparsest -> (row[1] + row[2]) + row[0]
        first_a = 1; first_b = 2; second_c = 0;
    }

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

    // Connect row_c to port A of second adder.
    for (int bit = 0; bit < width_abc; bit++) {
        npin_t *pin_a = (bit < width_c) ? row_c[bit] : get_zero_pin(netlist);
        add_input_pin_to_node(add2, pin_a, bit);
    }
    // Connect sumout_ab to port B of second adder.
    // This matches the pack pattern: adder[0].sumout -> adder[1].b
    for (int bit = 0; bit < width_abc; bit++) {
        npin_t *pin_b = (bit < width_ab) ? sumout_ab[bit] : get_zero_pin(netlist);
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

                // Connect row_c to port A
                for (int bit = 0; bit < width_abc; bit++) {
                    npin_t *pin_a = (bit < width_c) ? row_c[bit] : get_zero_pin(netlist);
                    add_input_pin_to_node(add2, pin_a, bit);
                }
                // Connect sumout_ab to port B
                // This matches the pack pattern: adder[0].sumout -> adder[1].b
                for (int bit = 0; bit < width_abc; bit++) {
                    npin_t *pin_b = (bit < width_ab) ? sumout_ab[bit] : get_zero_pin(netlist);
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