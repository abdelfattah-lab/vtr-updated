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
static signal_list_t *implement_compressor_tree_dadda(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks) ;

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
        case compressor_tree_type_e::DADDA:
            // use Dadda tree.
            return implement_compressor_tree_dadda(node, mark, netlist, ranks);
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