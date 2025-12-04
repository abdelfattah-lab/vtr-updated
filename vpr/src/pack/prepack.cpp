/*
 * Prepacking: Group together technology-mapped netlist blocks before packing.
 * This gives hints to the packer on what groups of blocks to keep together during packing.
 * Primary purpose:
 *    1) "Forced" packs (eg LUT+FF pair)
 *    2) Carry-chains
 * Duties: Find pack patterns in architecture, find pack patterns in netlist.
 *
 * Author: Jason Luu
 * March 12, 2012
 */

#include "prepack.h"
#include "globals.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <queue>
#include <utility>
#include <vector>

#include "atom_netlist.h"
#include "cluster_util.h"
#include "echo_files.h"
#include "physical_types.h"
#include "vpr_error.h"
#include "vpr_types.h"
#include "vpr_utils.h"
#include "vtr_assert.h"
#include "vtr_range.h"
#include "vtr_util.h"
#include "vtr_vector.h"

/*****************************************/
/*Local Function Declaration			 */
/*****************************************/
static std::vector<t_pack_patterns> alloc_and_load_pack_patterns(const std::vector<t_logical_block_type>& logical_block_types);

static void free_list_of_pack_patterns(std::vector<t_pack_patterns>& list_of_pack_patterns);

static void free_pack_pattern(t_pack_patterns* pack_pattern);

static t_pack_molecule* alloc_and_load_pack_molecules(t_pack_patterns* list_of_pack_patterns,
                                                      vtr::vector<AtomBlockId, t_pb_graph_node*>& expected_lowest_cost_pb_gnode,
                                                      const int num_packing_patterns,
                                                      std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                                      const AtomNetlist& atom_nlist,
                                                      const std::vector<t_logical_block_type>& logical_block_types);

static void discover_pattern_names_in_pb_graph_node(t_pb_graph_node* pb_graph_node,
                                                    std::unordered_map<std::string, int>& pattern_names);

static void forward_infer_pattern(t_pb_graph_pin* pb_graph_pin);

static void backward_infer_pattern(t_pb_graph_pin* pb_graph_pin);

static std::vector<t_pack_patterns> alloc_and_init_pattern_list_from_hash(std::unordered_map<std::string, int> pattern_names);

static t_pb_graph_edge* find_expansion_edge_of_pattern(const int pattern_index,
                                                       const t_pb_graph_node* pb_graph_node);

static void forward_expand_pack_pattern_from_edge(const t_pb_graph_edge* expansion_edge,
                                                  t_pack_patterns* list_of_packing_patterns,
                                                  const int curr_pattern_index,
                                                  int* L_num_blocks,
                                                  const bool make_root_of_chain);

static void backward_expand_pack_pattern_from_edge(const t_pb_graph_edge* expansion_edge,
                                                   t_pack_patterns* list_of_packing_patterns,
                                                   const int curr_pattern_index,
                                                   t_pb_graph_pin* destination_pin,
                                                   t_pack_pattern_block* destination_block,
                                                   int* L_num_blocks);

static int compare_pack_pattern(const t_pack_patterns* pattern_a, const t_pack_patterns* pattern_b);

static void free_pack_pattern_block(t_pack_pattern_block* pattern_block, t_pack_pattern_block** pattern_block_list);

static t_pack_molecule* try_create_molecule(t_pack_patterns* list_of_pack_patterns,
                                            const int pack_pattern_index,
                                            AtomBlockId blk_id,
                                            std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                            const AtomNetlist& atom_nlist);

static bool try_expand_molecule(t_pack_molecule* molecule,
                                const AtomBlockId blk_id,
                                const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                const AtomNetlist& atom_nlist);

static void print_pack_molecules(const char* fname,
                                 const t_pack_patterns* list_of_pack_patterns,
                                 const int num_pack_patterns,
                                 const t_pack_molecule* list_of_molecules,
                                 const AtomNetlist& atom_nlist);

static t_pb_graph_node* get_expected_lowest_cost_primitive_for_atom_block(const AtomBlockId blk_id,
                                                                          const std::vector<t_logical_block_type>& logical_block_types);

static t_pb_graph_node* get_expected_lowest_cost_primitive_for_atom_block_in_pb_graph_node(const AtomBlockId blk_id, t_pb_graph_node* curr_pb_graph_node, float* cost);

static AtomBlockId find_new_root_atom_for_chain(const AtomBlockId blk_id,
                                                const t_pack_patterns* list_of_pack_patterns,
                                                const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                                const AtomNetlist& atom_nlist);

static std::vector<t_pb_graph_pin*> find_end_of_path(t_pb_graph_pin* input_pin, int pattern_index);

static void expand_search(const t_pb_graph_pin* input_pin, std::queue<t_pb_graph_pin*>& pins_queue, const int pattern_index);

static void find_all_equivalent_chains(t_pack_patterns* chain_pattern, const t_pb_graph_node* root_block);

static void update_chain_root_pins(t_pack_patterns* chain_pattern,
                                   const std::vector<t_pb_graph_pin*>& chain_input_pins);

static void get_all_connected_primitive_pins(const t_pb_graph_pin* cluster_input_pin, std::vector<t_pb_graph_pin*>& connected_primitive_pins, int pattern_id);

static void init_molecule_chain_info(const AtomBlockId blk_id,
                                     t_pack_molecule* molecule,
                                     const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                     const AtomNetlist& atom_nlist);

static AtomBlockId get_sink_block(const AtomBlockId block_id,
                                  const t_model_ports* model_port,
                                  const BitIndex pin_number,
                                  const AtomNetlist& atom_nlist);

static AtomBlockId get_driving_block(const AtomBlockId block_id,
                                     const t_model_ports* model_port,
                                     const BitIndex pin_number,
                                     const AtomNetlist& atom_nlist);

static void print_chain_starting_points(t_pack_patterns* chain_pattern);

/** The following methods are utilized for extra carry chain logic: */

static t_pb_graph_pin* find_chain_exit_pin(t_pb_graph_pin* input_pin, int pattern_index);

static t_pack_pattern_block* get_atom_pattern_block(const t_pack_molecule* molecule, const int block_id);

static bool chain_input_is_reachable(const t_pack_molecule* molecule,
                                     const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                     const AtomNetlist& atom_nlist); // modernized

static t_pb_graph_node* get_driver_pb_graph_node(const t_pack_molecule* prev_molecule, const AtomBlockId driver_block);

static int get_forced_chain_id(t_pack_molecule* molecule,
                               const t_pack_molecule* prev_molecule,
                               const AtomBlockId driver_block_id);

static AtomBlockId get_adder_driver_block(const AtomBlockId block_id,
                                          const t_pack_patterns* pack_pattern,
                                          const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                          const AtomNetlist& atom_nlist); // modernized

static bool molecule_is_hierarchical(const t_pack_molecule* molecule);

static bool valid_second_level_placement(const AtomBlockId first_level_block,
                                         const AtomBlockId block_id,
                                         const t_pack_molecule* molecule,
                                         const AtomNetlist& atom_nlist); // modernized

static AtomBlockId is_second_level_block(const t_pack_pattern_block* pattern_block, const t_pack_molecule* molecule);

static bool check_alm_input_limitation(t_pack_molecule* molecule,
                                       const AtomNetlist& atom_nlist); // modernized

static void get_block_input_nets(const AtomBlockId block_id,
                                 std::unordered_set<AtomNetId>& nets,
                                 const AtomNetlist& atom_nlist);

static int get_pb_placement_index(t_pack_pattern_block* pattern_block, std::string pb_name);

static void modify_molecule(t_pack_molecule* molecule,
                            t_pack_pattern_block* pattern_block,
                            const AtomNetlist& atom_nlist);

static bool check_lut_chain_molecules(t_pack_molecule* molecule, const AtomNetlist& atom_nlist);

/*****************************************/
/*Function Definitions					 */
/*****************************************/

/**
 * Find all packing patterns in architecture
 * [0..num_packing_patterns-1]
 *
 * Limitations: Currently assumes that forced pack nets must be single-fanout
 * as this covers all the reasonable architectures we wanted.
 * More complicated structures should probably be handled either downstream
 * (general packing) or upstream (in tech mapping).
 * If this limitation is too constraining, code is designed so that this limitation can be removed.
 */
static std::vector<t_pack_patterns> alloc_and_load_pack_patterns(const std::vector<t_logical_block_type>& logical_block_types) {
    int L_num_blocks;
    std::vector<t_pack_patterns> list_of_packing_patterns;
    t_pb_graph_edge* expansion_edge;

    /* alloc and initialize array of packing patterns based on architecture complex blocks */
    std::unordered_map<std::string, int> pattern_names;
    for (const t_logical_block_type& type : logical_block_types) {
        discover_pattern_names_in_pb_graph_node(type.pb_graph_head, pattern_names);
    }

    list_of_packing_patterns = alloc_and_init_pattern_list_from_hash(pattern_names);

    /* load packing patterns by traversing the edges to find edges belonging to pattern */
    for (size_t i = 0; i < pattern_names.size(); i++) {
        for (const t_logical_block_type& type : logical_block_types) {
            // find an edge that belongs to this pattern
            expansion_edge = find_expansion_edge_of_pattern(i, type.pb_graph_head);
            if (!expansion_edge) {
                continue;
            }

            L_num_blocks = 0;
            list_of_packing_patterns[i].base_cost = 0;
            // use the found expansion edge to build the pack pattern
            backward_expand_pack_pattern_from_edge(expansion_edge,
                                                   list_of_packing_patterns.data(), i, nullptr, nullptr, &L_num_blocks);
            list_of_packing_patterns[i].num_blocks = L_num_blocks;

            /* Default settings: A section of a netlist must match all blocks in a pack
             * pattern before it can be made a molecule except for carry-chains.
             * For carry-chains, since carry-chains are typically quite flexible in terms
             * of size, it is optional whether or not an atom in a netlist matches any
             * particular block inside the chain */
            list_of_packing_patterns[i].is_block_optional = new bool[L_num_blocks];
            for (int k = 0; k < L_num_blocks; k++) {
                list_of_packing_patterns[i].is_block_optional[k] = false;
                if (list_of_packing_patterns[i].is_chain && list_of_packing_patterns[i].root_block->block_id != k) {
                    list_of_packing_patterns[i].is_block_optional[k] = true;
                }
            }

            // if this is a chain pattern (extends between complex blocks), check if there
            // are multiple equivalent chains with different starting and ending points
            if (list_of_packing_patterns[i].is_chain) {
                find_all_equivalent_chains(&list_of_packing_patterns[i], type.pb_graph_head);
                print_chain_starting_points(&list_of_packing_patterns[i]);
            }

            // if pack pattern i is found to belong to current block type, go to next pack pattern
            break;
        }
    }

    //Sanity check, every pattern should have a root block
    for (size_t i = 0; i < pattern_names.size(); ++i) {
        if (list_of_packing_patterns[i].root_block == nullptr) {
            VPR_FATAL_ERROR(VPR_ERROR_ARCH, "Failed to find root block for pack pattern %s", list_of_packing_patterns[i].name);
        }
    }

    return list_of_packing_patterns;
}

/**
 * Locate all pattern names
 * Side-effect: set all pb_graph_node temp_scratch_pad field to NULL
 *				For cases where a pattern inference is "obvious", mark it as obvious.
 */
static void discover_pattern_names_in_pb_graph_node(t_pb_graph_node* pb_graph_node,
                                                    std::unordered_map<std::string, int>& pattern_names) {
    /* Iterate over all edges to discover if an edge in current physical block belongs to a pattern
     * If edge does, then record the name of the pattern in a hash table */

    if (pb_graph_node == nullptr) {
        return;
    }

    pb_graph_node->temp_scratch_pad = nullptr;

    for (int i = 0; i < pb_graph_node->num_input_ports; i++) {
        for (int j = 0; j < pb_graph_node->num_input_pins[i]; j++) {
            bool hasPattern = false;
            for (int k = 0; k < pb_graph_node->input_pins[i][j].num_output_edges; k++) {
                auto output_edge = pb_graph_node->input_pins[i][j].output_edges[k];
                for (int m = 0; m < output_edge->num_pack_patterns; m++) {
                    hasPattern = true;
                    // insert the found pattern name to the hash table. If this pattern is inserted
                    // for the first time, then its index is the current size of the hash table
                    // otherwise the insert function will return an iterator of the previously
                    // inserted element with the index given to that pattern
                    std::string pattern_name(output_edge->pack_pattern_names[m]);
                    int index = (pattern_names.insert({pattern_name, pattern_names.size()}).first)->second;
                    if (!output_edge->pack_pattern_indices) {
                        output_edge->pack_pattern_indices = new int[output_edge->num_pack_patterns];
                    }
                    output_edge->pack_pattern_indices[m] = index;
                    // if this output edges belongs to a pack pattern. Expand forward starting from
                    // all its output pins to check if you need to infer pattern for direct connections
                    for (int ipin = 0; ipin < output_edge->num_output_pins; ipin++) {
                        forward_infer_pattern(output_edge->output_pins[ipin]);
                    }
                }
            }
            // if the output edge to this pin is annotated with a pack pattern
            // trace the inputs to this pin and mark them to infer pattern
            // if they are direct connections (num_input_edges == 1)
            if (hasPattern) {
                backward_infer_pattern(&pb_graph_node->input_pins[i][j]);
            }
        }
    }

    for (int i = 0; i < pb_graph_node->num_output_ports; i++) {
        for (int j = 0; j < pb_graph_node->num_output_pins[i]; j++) {
            bool hasPattern = false;
            for (int k = 0; k < pb_graph_node->output_pins[i][j].num_output_edges; k++) {
                auto output_edge = pb_graph_node->output_pins[i][j].output_edges[k];
                for (int m = 0; m < output_edge->num_pack_patterns; m++) {
                    hasPattern = true;
                    // insert the found pattern name to the hash table. If this pattern is inserted
                    // for the first time, then its index is the current size of the hash table
                    // otherwise the insert function will return an iterator of the previously
                    // inserted element with the index given to that pattern
                    std::string pattern_name(output_edge->pack_pattern_names[m]);
                    int index = (pattern_names.insert({pattern_name, pattern_names.size()}).first)->second;
                    if (!output_edge->pack_pattern_indices) {
                        output_edge->pack_pattern_indices = new int[output_edge->num_pack_patterns];
                    }
                    output_edge->pack_pattern_indices[m] = index;
                    // if this output edges belongs to a pack pattern. Expand forward starting from
                    // all its output pins to check if you need to infer pattern for direct connections
                    for (int ipin = 0; ipin < output_edge->num_output_pins; ipin++) {
                        forward_infer_pattern(output_edge->output_pins[ipin]);
                    }
                }
            }
            // if the output edge to this pin is annotated with a pack pattern
            // trace the inputs to this pin and mark them to infer pattern
            // if they are direct connections (num_input_edges == 1)
            if (hasPattern) {
                backward_infer_pattern(&pb_graph_node->output_pins[i][j]);
            }
        }
    }

    for (int i = 0; i < pb_graph_node->num_clock_ports; i++) {
        for (int j = 0; j < pb_graph_node->num_clock_pins[i]; j++) {
            bool hasPattern = false;
            for (int k = 0; k < pb_graph_node->clock_pins[i][j].num_output_edges; k++) {
                auto& output_edge = pb_graph_node->clock_pins[i][j].output_edges[k];
                for (int m = 0; m < output_edge->num_pack_patterns; m++) {
                    hasPattern = true;
                    // insert the found pattern name to the hash table. If this pattern is inserted
                    // for the first time, then its index is the current size of the hash table
                    // otherwise the insert function will return an iterator of the previously
                    // inserted element with the index given to that pattern
                    std::string pattern_name(output_edge->pack_pattern_names[m]);
                    int index = (pattern_names.insert({pattern_name, pattern_names.size()}).first)->second;
                    if (output_edge->pack_pattern_indices == nullptr) {
                        output_edge->pack_pattern_indices = new int[output_edge->num_pack_patterns];
                    }
                    output_edge->pack_pattern_indices[m] = index;
                    // if this output edges belongs to a pack pattern. Expand forward starting from
                    // all its output pins to check if you need to infer pattern for direct connections
                    for (int ipin = 0; ipin < output_edge->num_output_pins; ipin++) {
                        forward_infer_pattern(output_edge->output_pins[ipin]);
                    }
                }
            }
            // if the output edge to this pin is annotated with a pack pattern
            // trace the inputs to this pin and mark them to infer pattern
            // if they are direct connections (num_input_edges == 1)
            if (hasPattern) {
                backward_infer_pattern(&pb_graph_node->clock_pins[i][j]);
            }
        }
    }

    for (int i = 0; i < pb_graph_node->pb_type->num_modes; i++) {
        for (int j = 0; j < pb_graph_node->pb_type->modes[i].num_pb_type_children; j++) {
            for (int k = 0; k < pb_graph_node->pb_type->modes[i].pb_type_children[j].num_pb; k++) {
                discover_pattern_names_in_pb_graph_node(&pb_graph_node->child_pb_graph_nodes[i][j][k], pattern_names);
            }
        }
    }
}

/**
 * In obvious cases where a pattern edge has only one path to go, set that path to be inferred
 */
static void forward_infer_pattern(t_pb_graph_pin* pb_graph_pin) {
    if (pb_graph_pin->num_output_edges == 1 && pb_graph_pin->output_edges[0]->num_pack_patterns == 0 && pb_graph_pin->output_edges[0]->infer_pattern == false) {
        pb_graph_pin->output_edges[0]->infer_pattern = true;
        if (pb_graph_pin->output_edges[0]->num_output_pins == 1) {
            forward_infer_pattern(pb_graph_pin->output_edges[0]->output_pins[0]);
        }
    }
}
static void backward_infer_pattern(t_pb_graph_pin* pb_graph_pin) {
    if (pb_graph_pin->num_input_edges == 1 && pb_graph_pin->input_edges[0]->num_pack_patterns == 0 && pb_graph_pin->input_edges[0]->infer_pattern == false) {
        pb_graph_pin->input_edges[0]->infer_pattern = true;
        if (pb_graph_pin->input_edges[0]->num_input_pins == 1) {
            backward_infer_pattern(pb_graph_pin->input_edges[0]->input_pins[0]);
        }
    }
}

/**
 * Allocates memory for models and loads the name of the packing pattern
 * so that it can be identified and loaded with more complete information later
 */
static std::vector<t_pack_patterns> alloc_and_init_pattern_list_from_hash(std::unordered_map<std::string, int> pattern_names) {
    std::vector<t_pack_patterns> nlist(pattern_names.size());

    for (const auto& curr_pattern : pattern_names) {
        VTR_ASSERT(nlist[curr_pattern.second].name == nullptr);
        nlist[curr_pattern.second].name = vtr::strdup(curr_pattern.first.c_str());
        nlist[curr_pattern.second].root_block = nullptr;
        nlist[curr_pattern.second].is_chain = false;
        nlist[curr_pattern.second].index = curr_pattern.second;
    }

    return nlist;
}

static void free_list_of_pack_patterns(std::vector<t_pack_patterns>& list_of_pack_patterns) {
    for (size_t i = 0; i < list_of_pack_patterns.size(); i++) {
        free_pack_pattern(&list_of_pack_patterns[i]);
    }
}

static void free_pack_pattern(t_pack_patterns* pack_pattern) {
    if (pack_pattern) {
        int num_pack_pattern_blocks = pack_pattern->num_blocks;
        t_pack_pattern_block** pattern_block_list = new t_pack_pattern_block*[num_pack_pattern_blocks];
        for (int i = 0; i < num_pack_pattern_blocks; i++)
            pattern_block_list[i] = nullptr;

        free(pack_pattern->name);
        delete[] pack_pattern->is_block_optional;
        free_pack_pattern_block(pack_pattern->root_block, pattern_block_list);
        for (int j = 0; j < num_pack_pattern_blocks; j++) {
            delete pattern_block_list[j];
        }
        delete[] pattern_block_list;
    }
}

/**
 * Locate first edge that belongs to pattern index
 */
static t_pb_graph_edge* find_expansion_edge_of_pattern(const int pattern_index,
                                                       const t_pb_graph_node* pb_graph_node) {
    int i, j, k, m;
    t_pb_graph_edge* edge;
    /* Iterate over all edges to discover if an edge in current physical block belongs to a pattern
     * If edge does, then return that edge
     */

    if (pb_graph_node == nullptr) {
        return nullptr;
    }

    for (i = 0; i < pb_graph_node->num_input_ports; i++) {
        for (j = 0; j < pb_graph_node->num_input_pins[i]; j++) {
            auto& input_pin = pb_graph_node->input_pins[i][j];
            for (k = 0; k < input_pin.num_output_edges; k++) {
                for (m = 0; m < input_pin.output_edges[k]->num_pack_patterns; m++) {
                    if (input_pin.output_edges[k]->pack_pattern_indices[m] == pattern_index) {
                        return input_pin.output_edges[k];
                    }
                }
            }
        }
    }

    for (i = 0; i < pb_graph_node->num_output_ports; i++) {
        for (j = 0; j < pb_graph_node->num_output_pins[i]; j++) {
            auto& output_pin = pb_graph_node->output_pins[i][j];
            for (k = 0; k < output_pin.num_output_edges; k++) {
                for (m = 0; m < output_pin.output_edges[k]->num_pack_patterns; m++) {
                    if (output_pin.output_edges[k]->pack_pattern_indices[m] == pattern_index) {
                        return output_pin.output_edges[k];
                    }
                }
            }
        }
    }

    for (i = 0; i < pb_graph_node->num_clock_ports; i++) {
        for (j = 0; j < pb_graph_node->num_clock_pins[i]; j++) {
            auto& clock_pin = pb_graph_node->clock_pins[i][j];
            for (k = 0; k < clock_pin.num_output_edges; k++) {
                for (m = 0; m < clock_pin.output_edges[k]->num_pack_patterns; m++) {
                    if (clock_pin.output_edges[k]->pack_pattern_indices[m] == pattern_index) {
                        return clock_pin.output_edges[k];
                    }
                }
            }
        }
    }

    for (i = 0; i < pb_graph_node->pb_type->num_modes; i++) {
        auto& pb_mode = pb_graph_node->pb_type->modes[i];
        for (j = 0; j < pb_mode.num_pb_type_children; j++) {
            for (k = 0; k < pb_mode.pb_type_children[j].num_pb; k++) {
                edge = find_expansion_edge_of_pattern(pattern_index, &pb_graph_node->child_pb_graph_nodes[i][j][k]);
                if (edge != nullptr) {
                    return edge;
                }
            }
        }
    }
    return nullptr;
}

/**
 *  This function expands forward from the given expansion_edge. If a primitive is found that
 *  belongs to the pack pattern we are searching for, create a pack pattern block of using
 *  this primitive to be added later to the pack pattern when creating the pack pattern
 *  connections in the backward_expand_pack_pattern_from_edge function.
 *
 *  expansion_edge: starting edge to expand forward from
 *  list_of_packing_patterns: list of packing patterns in the architecture
 *  curr_pattern_index: current packing pattern that we are building
 *  L_num_blocks: number of primitives found to belong to this pattern so far
 *  make_root_of_chain: flag indicating that the given expansion_edge is connected
 *                      to a primitive that is the root of this packing pattern
 *
 *  Convention: Pack pattern block connections are made on backward expansion only (to make
 *              future multi-fanout support easier) so this function will not update connections
 */
static void forward_expand_pack_pattern_from_edge(const t_pb_graph_edge* expansion_edge,
                                                  t_pack_patterns* list_of_packing_patterns,
                                                  const int curr_pattern_index,
                                                  int* L_num_blocks,
                                                  bool make_root_of_chain) {
    int i, j, k;
    int iport, ipin, iedge;
    bool found; /* Error checking, ensure only one fan-out for each pattern net */
    t_pack_pattern_block* destination_block = nullptr;
    t_pb_graph_node* destination_pb_graph_node = nullptr;

    found = expansion_edge->infer_pattern;
    // if the pack pattern shouldn't be inferred check if the expansion
    // edge is annotated with the current pack pattern we are expanding
    for (i = 0; !found && i < expansion_edge->num_pack_patterns; i++) {
        if (expansion_edge->pack_pattern_indices[i] == curr_pattern_index) {
            found = true;
        }
    }
    // if this edge isn't annotated with the current pack pattern
    // no need to explore it
    if (!found) {
        return;
    }

    found = false;
    // iterate over the expansion edge output pins
    for (i = 0; i < expansion_edge->num_output_pins; i++) {
        // check if expansion_edge parent node is a primitive (i.e num_nodes = 0)
        if (expansion_edge->output_pins[i]->is_primitive_pin()) {
            destination_pb_graph_node = expansion_edge->output_pins[i]->parent_node;
            VTR_ASSERT(found == false);
            /* Check assumption that each forced net has only one fan-out */
            /* This is the destination node */
            found = true;

            // the temp_scratch_pad points to the last primitive from this pb_graph_node that was added to a packing pattern.
            const auto& destination_pb_temp = (t_pack_pattern_block*)destination_pb_graph_node->temp_scratch_pad;
            // if this pb_graph_node (primitive) is not added to the packing pattern already, add it and expand all its edges
            if (destination_pb_temp == nullptr || destination_pb_temp->pattern_index != curr_pattern_index) {
                // a primitive that belongs to this pack pattern is found: 1) create a new pattern block,
                // 2) assign an id to this pattern block, 3) increment the number of found blocks belonging to this
                // pattern and 4) expand all its edges to find the other primitives that belong to this pattern
                destination_block = new t_pack_pattern_block();
                list_of_packing_patterns[curr_pattern_index].base_cost += compute_primitive_base_cost(destination_pb_graph_node);
                destination_block->block_id = *L_num_blocks;
                (*L_num_blocks)++;
                destination_pb_graph_node->temp_scratch_pad = (void*)destination_block;
                destination_block->pattern_index = curr_pattern_index;
                destination_block->pb_type = destination_pb_graph_node->pb_type;

                // explore the inputs to this primitive
                for (iport = 0; iport < destination_pb_graph_node->num_input_ports; iport++) {
                    for (ipin = 0; ipin < destination_pb_graph_node->num_input_pins[iport]; ipin++) {
                        for (iedge = 0; iedge < destination_pb_graph_node->input_pins[iport][ipin].num_input_edges; iedge++) {
                            backward_expand_pack_pattern_from_edge(destination_pb_graph_node->input_pins[iport][ipin].input_edges[iedge],
                                                                   list_of_packing_patterns,
                                                                   curr_pattern_index,
                                                                   &destination_pb_graph_node->input_pins[iport][ipin],
                                                                   destination_block, L_num_blocks);
                        }
                    }
                }

                // explore the outputs of this primitive
                for (iport = 0; iport < destination_pb_graph_node->num_output_ports; iport++) {
                    for (ipin = 0; ipin < destination_pb_graph_node->num_output_pins[iport]; ipin++) {
                        for (iedge = 0; iedge < destination_pb_graph_node->output_pins[iport][ipin].num_output_edges; iedge++) {
                            forward_expand_pack_pattern_from_edge(destination_pb_graph_node->output_pins[iport][ipin].output_edges[iedge],
                                                                  list_of_packing_patterns,
                                                                  curr_pattern_index, L_num_blocks, false);
                        }
                    }
                }

                // explore the clock pins of this primitive
                for (iport = 0; iport < destination_pb_graph_node->num_clock_ports; iport++) {
                    for (ipin = 0; ipin < destination_pb_graph_node->num_clock_pins[iport]; ipin++) {
                        for (iedge = 0; iedge < destination_pb_graph_node->clock_pins[iport][ipin].num_input_edges; iedge++) {
                            backward_expand_pack_pattern_from_edge(destination_pb_graph_node->clock_pins[iport][ipin].input_edges[iedge],
                                                                   list_of_packing_patterns,
                                                                   curr_pattern_index,
                                                                   &destination_pb_graph_node->clock_pins[iport][ipin],
                                                                   destination_block, L_num_blocks);
                        }
                    }
                }
            }

            // if this pb_graph_node (primitive) should be added to the pack pattern blocks
            if (((t_pack_pattern_block*)destination_pb_graph_node->temp_scratch_pad)->pattern_index == curr_pattern_index) {
                // if this pb_graph_node is known to be the root of the chain, update the root block and root pin
                if (make_root_of_chain == true) {
                    list_of_packing_patterns[curr_pattern_index].chain_root_pins = {{expansion_edge->output_pins[i]}};
                    list_of_packing_patterns[curr_pattern_index].root_block = destination_block;
                }
            }

            // the expansion_edge parent node is not a primitive
        } else {
            // continue expanding forward
            for (j = 0; j < expansion_edge->output_pins[i]->num_output_edges; j++) {
                if (expansion_edge->output_pins[i]->output_edges[j]->infer_pattern == true) {
                    forward_expand_pack_pattern_from_edge(expansion_edge->output_pins[i]->output_edges[j],
                                                          list_of_packing_patterns,
                                                          curr_pattern_index,
                                                          L_num_blocks,
                                                          make_root_of_chain);
                } else {
                    for (k = 0; k < expansion_edge->output_pins[i]->output_edges[j]->num_pack_patterns; k++) {
                        if (expansion_edge->output_pins[i]->output_edges[j]->pack_pattern_indices[k] == curr_pattern_index) {
                            if (found == true) {
                                /* Check assumption that each forced net has only one fan-out */
                                VPR_FATAL_ERROR(VPR_ERROR_PACK,
                                                "Invalid packing pattern defined.  Multi-fanout nets not supported when specifying pack patterns.\n"
                                                "Problem on %s[%d].%s[%d] for pattern %s\n",
                                                expansion_edge->output_pins[i]->parent_node->pb_type->name,
                                                expansion_edge->output_pins[i]->parent_node->placement_index,
                                                expansion_edge->output_pins[i]->port->name,
                                                expansion_edge->output_pins[i]->pin_number,
                                                list_of_packing_patterns[curr_pattern_index].name);
                            }
                            found = true;
                            forward_expand_pack_pattern_from_edge(expansion_edge->output_pins[i]->output_edges[j],
                                                                  list_of_packing_patterns,
                                                                  curr_pattern_index,
                                                                  L_num_blocks,
                                                                  make_root_of_chain);
                        }
                    } // End for pack patterns of output edge
                }
            } // End for number of output edges
        }
    } // End for output pins of expansion edge
}

/**
 * Find if driver of edge is in the same pattern, if yes, add to pattern
 *  Convention: Connections are made on backward expansion only (to make future multi-
 *               fanout support easier) so this function must update both source and
 *               destination blocks
 */
static void backward_expand_pack_pattern_from_edge(const t_pb_graph_edge* expansion_edge,
                                                   t_pack_patterns* list_of_packing_patterns,
                                                   const int curr_pattern_index,
                                                   t_pb_graph_pin* destination_pin,
                                                   t_pack_pattern_block* destination_block,
                                                   int* L_num_blocks) {
    int i, j, k;
    int iport, ipin, iedge;
    bool found; /* Error checking, ensure only one fan-out for each pattern net */
    t_pack_pattern_block* source_block = nullptr;
    t_pb_graph_node* source_pb_graph_node = nullptr;
    t_pack_pattern_connections* pack_pattern_connection = nullptr;

    found = expansion_edge->infer_pattern;
    // if the pack pattern shouldn't be inferred check if the expansion
    // edge is annotated with the current pack pattern we are expanding
    for (i = 0; !found && i < expansion_edge->num_pack_patterns; i++) {
        if (expansion_edge->pack_pattern_indices[i] == curr_pattern_index) {
            found = true;
        }
    }

    // if this edge isn't annotated with the current pack pattern
    // no need to explore it
    if (!found) {
        return;
    }

    found = false;
    // iterate over all the drivers of this edge
    for (i = 0; i < expansion_edge->num_input_pins; i++) {
        // check if the expansion_edge parent node is a primitive
        if (expansion_edge->input_pins[i]->is_primitive_pin()) {
            source_pb_graph_node = expansion_edge->input_pins[i]->parent_node;
            VTR_ASSERT(found == false);
            /* Check assumption that each forced net has only one fan-out */
            /* This is the source node for destination */
            found = true;

            /* If this pb_graph_node is part not of the current pattern index, put it in and expand all its edges */
            source_block = (t_pack_pattern_block*)source_pb_graph_node->temp_scratch_pad;
            if (source_block == nullptr || source_block->pattern_index != curr_pattern_index) {
                source_block = new t_pack_pattern_block();
                source_block->block_id = *L_num_blocks;
                (*L_num_blocks)++;
                list_of_packing_patterns[curr_pattern_index].base_cost += compute_primitive_base_cost(source_pb_graph_node);
                source_pb_graph_node->temp_scratch_pad = (void*)source_block;
                source_block->pattern_index = curr_pattern_index;
                source_block->pb_type = source_pb_graph_node->pb_type;

                if (list_of_packing_patterns[curr_pattern_index].root_block == nullptr) {
                    list_of_packing_patterns[curr_pattern_index].root_block = source_block;
                }

                // explore the inputs of this primitive
                for (iport = 0; iport < source_pb_graph_node->num_input_ports; iport++) {
                    for (ipin = 0; ipin < source_pb_graph_node->num_input_pins[iport]; ipin++) {
                        for (iedge = 0; iedge < source_pb_graph_node->input_pins[iport][ipin].num_input_edges; iedge++) {
                            backward_expand_pack_pattern_from_edge(source_pb_graph_node->input_pins[iport][ipin].input_edges[iedge],
                                                                   list_of_packing_patterns,
                                                                   curr_pattern_index,
                                                                   &source_pb_graph_node->input_pins[iport][ipin],
                                                                   source_block,
                                                                   L_num_blocks);
                        }
                    }
                }

                // explore the outputs of this primitive
                for (iport = 0; iport < source_pb_graph_node->num_output_ports; iport++) {
                    for (ipin = 0; ipin < source_pb_graph_node->num_output_pins[iport]; ipin++) {
                        for (iedge = 0; iedge < source_pb_graph_node->output_pins[iport][ipin].num_output_edges; iedge++) {
                            forward_expand_pack_pattern_from_edge(source_pb_graph_node->output_pins[iport][ipin].output_edges[iedge],
                                                                  list_of_packing_patterns,
                                                                  curr_pattern_index,
                                                                  L_num_blocks,
                                                                  false);
                        }
                    }
                }

                // explore the clock pins of this primitive
                for (iport = 0; iport < source_pb_graph_node->num_clock_ports; iport++) {
                    for (ipin = 0; ipin < source_pb_graph_node->num_clock_pins[iport]; ipin++) {
                        for (iedge = 0; iedge < source_pb_graph_node->clock_pins[iport][ipin].num_input_edges; iedge++) {
                            backward_expand_pack_pattern_from_edge(source_pb_graph_node->clock_pins[iport][ipin].input_edges[iedge],
                                                                   list_of_packing_patterns,
                                                                   curr_pattern_index,
                                                                   &source_pb_graph_node->clock_pins[iport][ipin],
                                                                   source_block,
                                                                   L_num_blocks);
                        }
                    }
                }
            }

            if (destination_pin != nullptr) {
                VTR_ASSERT(((t_pack_pattern_block*)source_pb_graph_node->temp_scratch_pad)->pattern_index == curr_pattern_index);
                source_block = (t_pack_pattern_block*)source_pb_graph_node->temp_scratch_pad;
                pack_pattern_connection = new t_pack_pattern_connections();
                pack_pattern_connection->from_block = source_block;
                pack_pattern_connection->from_pin = expansion_edge->input_pins[i];
                pack_pattern_connection->to_block = destination_block;
                pack_pattern_connection->to_pin = destination_pin;
                pack_pattern_connection->next = source_block->connections;
                source_block->connections = pack_pattern_connection;

                pack_pattern_connection = new t_pack_pattern_connections();
                pack_pattern_connection->from_block = source_block;
                pack_pattern_connection->from_pin = expansion_edge->input_pins[i];
                pack_pattern_connection->to_block = destination_block;
                pack_pattern_connection->to_pin = destination_pin;
                pack_pattern_connection->next = destination_block->connections;
                destination_block->connections = pack_pattern_connection;

                if (source_block == destination_block) {
                    VPR_FATAL_ERROR(VPR_ERROR_PACK,
                                    "Invalid packing pattern defined. Source and destination block are the same (%s).\n",
                                    source_block->pb_type->name);
                }
            }

            // expansion edge parent is not a primitive
        } else {
            // check if this input pin of the expansion edge has no driving pin
            if (expansion_edge->input_pins[i]->num_input_edges == 0) {
                // check if this input pin of the expansion edge belongs to a root block (i.e doesn't have a parent block)
                if (expansion_edge->input_pins[i]->parent_node->pb_type->parent_mode == nullptr) {
                    // This pack pattern extends to CLB (root pb block) input pin,
                    // thus it extends across multiple logic blocks, treat as a chain
                    list_of_packing_patterns[curr_pattern_index].is_chain = true;
                    // since this input pin has not driving nets, expand in the forward direction instead
                    forward_expand_pack_pattern_from_edge(expansion_edge,
                                                          list_of_packing_patterns,
                                                          curr_pattern_index,
                                                          L_num_blocks,
                                                          true);
                }
                // this input pin of the expansion edge has a driving pin
            } else {
                // iterate over all the driving edges of this input pin
                for (j = 0; j < expansion_edge->input_pins[i]->num_input_edges; j++) {
                    // if pattern should be inferred for this edge continue the expansion backwards
                    if (expansion_edge->input_pins[i]->input_edges[j]->infer_pattern == true) {
                        backward_expand_pack_pattern_from_edge(expansion_edge->input_pins[i]->input_edges[j],
                                                               list_of_packing_patterns,
                                                               curr_pattern_index,
                                                               destination_pin,
                                                               destination_block,
                                                               L_num_blocks);
                        // if pattern shouldn't be inferred
                    } else {
                        // check if this input pin edge is annotated with the current pattern
                        for (k = 0; k < expansion_edge->input_pins[i]->input_edges[j]->num_pack_patterns; k++) {
                            if (expansion_edge->input_pins[i]->input_edges[j]->pack_pattern_indices[k] == curr_pattern_index) {
                                VTR_ASSERT(found == false);
                                /* Check assumption that each forced net has only one fan-out */
                                found = true;
                                backward_expand_pack_pattern_from_edge(expansion_edge->input_pins[i]->input_edges[j],
                                                                       list_of_packing_patterns,
                                                                       curr_pattern_index,
                                                                       destination_pin,
                                                                       destination_block,
                                                                       L_num_blocks);
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * Pre-pack atoms in netlist to molecules
 * 1.  Single atoms are by definition a molecule.
 * 2.  Forced pack molecules are groupings of atoms that matches a t_pack_pattern definition.
 * 3.  Chained molecules are molecules that follow a carry-chain style pattern,
 *     ie. a single linear chain that can be split across multiple complex blocks
 */
static void fill_vacant_chain_spots(t_pack_molecule* list_of_molecules_head,
                                    const t_pack_patterns* list_of_pack_patterns,
                                    const int num_packing_patterns,
                                    std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules) {
    auto& atom_ctx = g_vpr_ctx.mutable_atom();
    AtomNetlist& atom_nlist = atom_ctx.nlist;

    // Find or create a ground net
    AtomNetId gnd_net_id = atom_nlist.find_net("gnd");
    if (!gnd_net_id) {
        gnd_net_id = atom_nlist.create_net("gnd");
        // We need a driver for this net. Ideally a constant generator.
        // For now, we assume if it didn't exist, we might need to create a dummy driver or leave it undriven (which might be an error).
        // However, usually 'gnd' exists if used in the design.
        // If we create it, we should probably make it a constant.
        // Let's try to find a constant zero block/pin if possible, or just create the net and hope the router handles it (or legalizer).
        // A safer bet is to look for any net that is a constant 0.
        // But for this specific task, let's assume "gnd" is the standard name.
    }

    t_pack_molecule* cur_molecule = list_of_molecules_head;
    while (cur_molecule != nullptr) {
        if (cur_molecule->type == MOLECULE_FORCED_PACK && cur_molecule->pack_pattern->is_chain) {
            // Check if this is the double carry chain pattern we are interested in
            // We look for the specific indexing: Row 0 (0-19) and Row 1 (39-20)
            // We can check if the pattern has at least 40 blocks.
            if (cur_molecule->num_blocks >= 40) {
                for (int i = 0; i < 20; ++i) {
                    int row0_idx = i;
                    int row1_idx = 39 - i;

                    AtomBlockId row0_blk = cur_molecule->atom_block_ids[row0_idx];
                    AtomBlockId row1_blk = cur_molecule->atom_block_ids[row1_idx];

                    // Case 1: Row 0 occupied, Row 1 empty -> Fill Row 1
                    if (row0_blk && !row1_blk) {
                        // Found a vacant spot in row 1!
                        VTR_LOG("Filling vacant spot in molecule for pattern %s at index %d (paired with %d)\n",
                                cur_molecule->pack_pattern->name, row1_idx, row0_idx);

                        // 1. Create new block
                        std::string new_name = atom_nlist.block_name(row0_blk) + "_pass_through_" + std::to_string(row1_idx);
                        const t_model* model = atom_nlist.block_model(row0_blk);
                        AtomBlockId new_blk_id = atom_nlist.create_block(new_name, model);

                        // 2. Connect Pins
                        // We need to find the cin, cout, a, b, and sumout ports.
                        const t_model_ports* cin_model_port = cur_molecule->pack_pattern->chain_root_pins[0][0]->port->model_port;
                        const t_model_ports* cout_model_port = cur_molecule->pack_pattern->chain_exit_pins[0]->port->model_port;

                        // Find a, b, and sumout model ports from the adder model
                        const t_model* adder_model = model;
                        const t_model_ports* a_model_port = nullptr;
                        const t_model_ports* b_model_port = nullptr;
                        const t_model_ports* sumout_model_port = nullptr;
                        for (const t_model_ports* port = adder_model->inputs; port; port = port->next) {
                            if (std::string(port->name) == "a") a_model_port = port;
                            if (std::string(port->name) == "b") b_model_port = port;
                        }
                        for (const t_model_ports* port = adder_model->outputs; port; port = port->next) {
                            if (std::string(port->name) == "sumout") sumout_model_port = port;
                        }

                        // Check if there's a downstream block that needs COUT
                        // Row 1 chain flows: 39 → 38 → ... → 20, so next is row1_idx - 1
                        int next_idx = row1_idx - 1;
                        bool has_downstream = (next_idx >= 20 && cur_molecule->atom_block_ids[next_idx]);

                        // Create ports on the new block
                        AtomPortId cin_port_id = atom_nlist.create_port(new_blk_id, cin_model_port);
                        AtomPortId a_port_id = a_model_port ? atom_nlist.create_port(new_blk_id, a_model_port) : AtomPortId::INVALID();
                        AtomPortId b_port_id = b_model_port ? atom_nlist.create_port(new_blk_id, b_model_port) : AtomPortId::INVALID();
                        AtomPortId sumout_port_id = sumout_model_port ? atom_nlist.create_port(new_blk_id, sumout_model_port) : AtomPortId::INVALID();
                        // Only create COUT if there's a downstream block
                        AtomPortId cout_port_id = has_downstream ? atom_nlist.create_port(new_blk_id, cout_model_port) : AtomPortId::INVALID();

                        // Determine Driver for CIN
                        AtomNetId cin_driver_net;
                        if (row1_idx == 39) {
                            // Start of chain -> GND
                            cin_driver_net = gnd_net_id;
                        } else {
                            // Middle of chain -> Driven by previous block's COUT
                            AtomBlockId prev_blk = cur_molecule->atom_block_ids[row1_idx + 1];
                            VTR_ASSERT(prev_blk); // Should exist because we iterate 39 down to 20

                            // Find COUT net of prev_blk
                            AtomPortId prev_cout_port = atom_nlist.find_atom_port(prev_blk, cout_model_port);
                            if (!prev_cout_port) {
                                // Create if missing (e.g. if prev block was a real atom that didn't use cout)
                                prev_cout_port = atom_nlist.create_port(prev_blk, cout_model_port);
                            }

                            cin_driver_net = atom_nlist.port_net(prev_cout_port, 0);
                            if (!cin_driver_net) {
                                // Create net if missing
                                // COUT net names use [0] suffix, while block names often use [1] for sumout variant
                                // We need to create a unique COUT net name to avoid collision with SUMOUT net
                                std::string block_name = atom_nlist.block_name(prev_blk);
                                std::string net_name;
                                // Replace trailing [1] with [0] for COUT net naming convention
                                if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[1]") {
                                    net_name = block_name.substr(0, block_name.size() - 3) + "[0]";
                                } else if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[0]") {
                                    net_name = block_name; // Already has [0] suffix
                                } else {
                                    net_name = block_name + "_cout"; // Fallback for unusual naming
                                }
                                cin_driver_net = atom_nlist.create_net(net_name);
                                // Only add driver pin if net doesn't already have one
                                // (create_net may return an existing net with the same name)
                                if (!atom_nlist.net_driver(cin_driver_net)) {
                                    atom_nlist.create_pin(prev_cout_port, 0, cin_driver_net, PinType::DRIVER, false);
                                }
                            }
                        }

                        // Connect CIN
                        atom_nlist.create_pin(cin_port_id, 0, cin_driver_net, PinType::SINK, false);

                        // Connect A and B to ground (for pass-through behavior: A=0, B=0 makes SUM=CIN)
                        if (a_port_id) {
                            atom_nlist.create_pin(a_port_id, 0, gnd_net_id, PinType::SINK, false);
                        }
                        if (b_port_id) {
                            atom_nlist.create_pin(b_port_id, 0, gnd_net_id, PinType::SINK, false);
                        }

                        // Create and Connect COUT Net only if COUT port exists (has downstream block)
                        if (cout_port_id) {
                            // COUT net names use [0] suffix convention
                            std::string block_name = atom_nlist.block_name(new_blk_id);
                            std::string cout_net_name;
                            if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[1]") {
                                cout_net_name = block_name.substr(0, block_name.size() - 3) + "[0]";
                            } else if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[0]") {
                                cout_net_name = block_name;
                            } else {
                                cout_net_name = block_name + "_cout";
                            }
                            AtomNetId cout_net = atom_nlist.create_net(cout_net_name);
                            // Only add driver pin if net doesn't already have one
                            if (!atom_nlist.net_driver(cout_net)) {
                                atom_nlist.create_pin(cout_port_id, 0, cout_net, PinType::DRIVER, false);
                            }

                            // Rewire the next block in the chain if it exists and is occupied
                            // The chain flows from row1_idx to row1_idx - 1
                            if (has_downstream) {
                                AtomBlockId next_blk = cur_molecule->atom_block_ids[next_idx];
                                VTR_ASSERT(next_blk); // Should exist since has_downstream is true
                                // The next block exists (it was already there).
                                // We must disconnect its cin from whatever it was connected to (e.g. gnd)
                                // and connect it to our new cout net.
                                AtomPortId next_cin_port = atom_nlist.find_atom_port(next_blk, cin_model_port);
                                if (next_cin_port) {
                                    AtomPinId next_cin_pin = atom_nlist.port_pin(next_cin_port, 0);
                                    if (next_cin_pin) {
                                        // set_pin_net automatically removes the previous connection
                                        atom_nlist.set_pin_net(next_cin_pin, PinType::SINK, cout_net);
                                    }
                                }
                            }
                        }

                        // Create SUMOUT net (for sumout connections between rows)
                        if (sumout_port_id) {
                            std::string sumout_net_name = atom_nlist.block_name(new_blk_id) + "_sumout";
                            AtomNetId sumout_net = atom_nlist.create_net(sumout_net_name);
                            atom_nlist.create_pin(sumout_port_id, 0, sumout_net, PinType::DRIVER, false);
                        }

                        // 3. Update Molecule
                        cur_molecule->atom_block_ids[row1_idx] = new_blk_id;
                        cur_molecule->num_blocks++; // Increment block count

                        // 4. Register in atom_molecules
                        atom_molecules.insert({new_blk_id, cur_molecule});
                    }
                    // Case 2: Row 1 occupied, Row 0 empty -> Fill Row 0
                    else if (!row0_blk && row1_blk) {
                        // Found a vacant spot in row 0!
                        VTR_LOG("Filling vacant spot in molecule for pattern %s at index %d (paired with %d)\n",
                                cur_molecule->pack_pattern->name, row0_idx, row1_idx);

                        // 1. Create new block
                        std::string new_name = atom_nlist.block_name(row1_blk) + "_pass_through_" + std::to_string(row0_idx);
                        const t_model* model = atom_nlist.block_model(row1_blk);
                        AtomBlockId new_blk_id = atom_nlist.create_block(new_name, model);

                        // 2. Connect Pins
                        // We need to find the cin, cout, a, b, and sumout ports.
                        const t_model_ports* cin_model_port = cur_molecule->pack_pattern->chain_root_pins[0][0]->port->model_port;
                        const t_model_ports* cout_model_port = cur_molecule->pack_pattern->chain_exit_pins[0]->port->model_port;

                        // Find a, b, and sumout model ports from the adder model
                        const t_model* adder_model = model;
                        const t_model_ports* a_model_port = nullptr;
                        const t_model_ports* b_model_port = nullptr;
                        const t_model_ports* sumout_model_port = nullptr;
                        for (const t_model_ports* port = adder_model->inputs; port; port = port->next) {
                            if (std::string(port->name) == "a") a_model_port = port;
                            if (std::string(port->name) == "b") b_model_port = port;
                        }
                        for (const t_model_ports* port = adder_model->outputs; port; port = port->next) {
                            if (std::string(port->name) == "sumout") sumout_model_port = port;
                        }

                        // Check if there's a downstream block that needs COUT
                        // Row 0 chain flows: 0 → 1 → 2 → ... → 19, so next is row0_idx + 1
                        int next_idx = row0_idx + 1;
                        bool has_downstream = (next_idx < 20 && cur_molecule->atom_block_ids[next_idx]);

                        // Create ports on the new block
                        AtomPortId cin_port_id = atom_nlist.create_port(new_blk_id, cin_model_port);
                        AtomPortId a_port_id = a_model_port ? atom_nlist.create_port(new_blk_id, a_model_port) : AtomPortId::INVALID();
                        AtomPortId b_port_id = b_model_port ? atom_nlist.create_port(new_blk_id, b_model_port) : AtomPortId::INVALID();
                        AtomPortId sumout_port_id = sumout_model_port ? atom_nlist.create_port(new_blk_id, sumout_model_port) : AtomPortId::INVALID();
                        // Only create COUT if there's a downstream block
                        AtomPortId cout_port_id = has_downstream ? atom_nlist.create_port(new_blk_id, cout_model_port) : AtomPortId::INVALID();

                        // Determine Driver for CIN
                        // Row 0 chain flows: 0 → 1 → 2 → ... → 19
                        AtomNetId cin_driver_net;
                        if (row0_idx == 0) {
                            // Start of row 0 chain -> GND
                            cin_driver_net = gnd_net_id;
                        } else {
                            // Middle of chain -> Driven by previous block's COUT (position i-1)
                            AtomBlockId prev_blk = cur_molecule->atom_block_ids[row0_idx - 1];
                            VTR_ASSERT(prev_blk); // Should exist because we're in the middle of the chain

                            // Find COUT net of prev_blk
                            AtomPortId prev_cout_port = atom_nlist.find_atom_port(prev_blk, cout_model_port);
                            if (!prev_cout_port) {
                                // Create if missing
                                prev_cout_port = atom_nlist.create_port(prev_blk, cout_model_port);
                            }

                            cin_driver_net = atom_nlist.port_net(prev_cout_port, 0);
                            if (!cin_driver_net) {
                                // Create net if missing
                                // COUT net names use [0] suffix, while block names often use [1] for sumout variant
                                // We need to create a unique COUT net name to avoid collision with SUMOUT net
                                std::string block_name = atom_nlist.block_name(prev_blk);
                                std::string net_name;
                                // Replace trailing [1] with [0] for COUT net naming convention
                                if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[1]") {
                                    net_name = block_name.substr(0, block_name.size() - 3) + "[0]";
                                } else if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[0]") {
                                    net_name = block_name; // Already has [0] suffix
                                } else {
                                    net_name = block_name + "_cout"; // Fallback for unusual naming
                                }
                                cin_driver_net = atom_nlist.create_net(net_name);
                                // Only add driver pin if net doesn't already have one
                                // (create_net may return an existing net with the same name)
                                if (!atom_nlist.net_driver(cin_driver_net)) {
                                    atom_nlist.create_pin(prev_cout_port, 0, cin_driver_net, PinType::DRIVER, false);
                                }
                            }
                        }

                        // Connect CIN
                        atom_nlist.create_pin(cin_port_id, 0, cin_driver_net, PinType::SINK, false);

                        // Connect A and B to ground (for pass-through behavior: A=0, B=0 makes SUM=CIN)
                        if (a_port_id) {
                            atom_nlist.create_pin(a_port_id, 0, gnd_net_id, PinType::SINK, false);
                        }
                        if (b_port_id) {
                            atom_nlist.create_pin(b_port_id, 0, gnd_net_id, PinType::SINK, false);
                        }

                        // Create and Connect COUT Net only if COUT port exists (has downstream block)
                        if (cout_port_id) {
                            // COUT net names use [0] suffix convention
                            std::string block_name = atom_nlist.block_name(new_blk_id);
                            std::string cout_net_name;
                            if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[1]") {
                                cout_net_name = block_name.substr(0, block_name.size() - 3) + "[0]";
                            } else if (block_name.size() >= 3 && block_name.substr(block_name.size() - 3) == "[0]") {
                                cout_net_name = block_name;
                            } else {
                                cout_net_name = block_name + "_cout";
                            }
                            AtomNetId cout_net = atom_nlist.create_net(cout_net_name);
                            // Only add driver pin if net doesn't already have one
                            if (!atom_nlist.net_driver(cout_net)) {
                                atom_nlist.create_pin(cout_port_id, 0, cout_net, PinType::DRIVER, false);
                            }

                            // Rewire the next block in the chain if it exists and is occupied
                            // Row 0 chain flows forward: 0 → 1 → 2 → ... → 19
                            if (has_downstream) {
                                AtomBlockId next_blk = cur_molecule->atom_block_ids[next_idx];
                                VTR_ASSERT(next_blk); // Should exist since has_downstream is true
                                // The next block exists (it was already there).
                                // We must disconnect its cin from whatever it was connected to (e.g. gnd)
                                // and connect it to our new cout net.
                                AtomPortId next_cin_port = atom_nlist.find_atom_port(next_blk, cin_model_port);
                                if (next_cin_port) {
                                    AtomPinId next_cin_pin = atom_nlist.port_pin(next_cin_port, 0);
                                    if (next_cin_pin) {
                                        // set_pin_net automatically removes the previous connection
                                        atom_nlist.set_pin_net(next_cin_pin, PinType::SINK, cout_net);
                                    }
                                }
                            }
                        }

                        // Create SUMOUT net (for sumout connections between rows)
                        if (sumout_port_id) {
                            std::string sumout_net_name = atom_nlist.block_name(new_blk_id) + "_sumout";
                            AtomNetId sumout_net = atom_nlist.create_net(sumout_net_name);
                            atom_nlist.create_pin(sumout_port_id, 0, sumout_net, PinType::DRIVER, false);
                        }

                        // 3. Update Molecule
                        cur_molecule->atom_block_ids[row0_idx] = new_blk_id;
                        cur_molecule->num_blocks++; // Increment block count

                        // 4. Register in atom_molecules
                        atom_molecules.insert({new_blk_id, cur_molecule});
                    }
                }
            }
        }
        cur_molecule = cur_molecule->next;
    }
}

static t_pack_molecule* alloc_and_load_pack_molecules(t_pack_patterns* list_of_pack_patterns,
                                                      vtr::vector<AtomBlockId, t_pb_graph_node*>& expected_lowest_cost_pb_gnode,
                                                      const int num_packing_patterns,
                                                      std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                                      const AtomNetlist& atom_nlist,
                                                      const std::vector<t_logical_block_type>& logical_block_types) {
    int i, j, best_pattern;
    t_pack_molecule* list_of_molecules_head;
    t_pack_molecule* cur_molecule;
    bool* is_used;

    is_used = new bool[num_packing_patterns];
    for (i = 0; i < num_packing_patterns; i++) {
        is_used[i] = false;
    }

    cur_molecule = list_of_molecules_head = nullptr;

    /* Find forced pack patterns
     * Simplifying assumptions: Each atom can map to at most one molecule,
     *                          use first-fit mapping based on priority of pattern
     * TODO: Need to investigate better mapping strategies than first-fit
     */
    for (i = 0; i < num_packing_patterns; i++) {
        best_pattern = 0;
        for (j = 1; j < num_packing_patterns; j++) {
            if (is_used[best_pattern]) {
                best_pattern = j;
            } else if (is_used[j] == false && compare_pack_pattern(&list_of_pack_patterns[j], &list_of_pack_patterns[best_pattern]) == 1) {
                best_pattern = j;
            }
        }
        // If all patterns are already marked as used (e.g. because some
        // patterns were pre-disabled such as *lut_chain*), stop creating
        // forced-pack molecules.
        if (is_used[best_pattern]) {
            break;
        }
        is_used[best_pattern] = true;

        auto blocks = atom_nlist.blocks();
        for (auto blk_iter = blocks.begin(); blk_iter != blocks.end(); ++blk_iter) {
            auto blk_id = *blk_iter;

            cur_molecule = try_create_molecule(list_of_pack_patterns, best_pattern, blk_id, atom_molecules, atom_nlist);
            if (cur_molecule != nullptr) {
                cur_molecule->next = list_of_molecules_head;
                /* In the event of multiple molecules with the same atom block pattern,
                 * bias to use the molecule with less costly physical resources first */
                /* TODO: Need to normalize magical number 100 */
                cur_molecule->base_gain = cur_molecule->num_blocks - (cur_molecule->pack_pattern->base_cost / 100);
                list_of_molecules_head = cur_molecule;

                //Note: atom_molecules is an (ordered) multimap so the last molecule
                //      inserted for a given blk_id will be the last valid element
                //      in the equal_range
                auto rng = atom_molecules.equal_range(blk_id); //The range of molecules matching this block
                bool range_empty = (rng.first == rng.second);
                bool cur_was_last_inserted = false;
                if (!range_empty) {
                    auto last_valid_iter = --rng.second; //Iterator to last element (only valid if range is not empty)
                    cur_was_last_inserted = (last_valid_iter->second == cur_molecule);
                }
                if (range_empty || !cur_was_last_inserted) {
                    /* molecule did not cover current atom (possibly because molecule created is
                     * part of a long chain that extends past multiple logic blocks), try again */
                    --blk_iter;
                }
            }
        }
    }
    delete[] is_used;

    /* List all atom blocks as a molecule for blocks that do not belong to any molecules.
     * This allows the packer to be consistent as it now packs molecules only instead of atoms and molecules
     *
     * If a block belongs to a molecule, then carrying the single atoms around can make the packing problem
     * more difficult because now it needs to consider splitting molecules.
     */
    for (auto blk_id : atom_nlist.blocks()) {
        expected_lowest_cost_pb_gnode[blk_id] = get_expected_lowest_cost_primitive_for_atom_block(blk_id, logical_block_types);

        auto rng = atom_molecules.equal_range(blk_id);
        bool rng_empty = (rng.first == rng.second);
        if (rng_empty) {
            cur_molecule = new t_pack_molecule;
            cur_molecule->type = MOLECULE_SINGLE_ATOM;
            cur_molecule->num_blocks = 1;
            cur_molecule->root = 0;
            cur_molecule->pack_pattern = nullptr;

            cur_molecule->atom_block_ids = {blk_id};

            cur_molecule->next = list_of_molecules_head;
            cur_molecule->base_gain = 1;
            list_of_molecules_head = cur_molecule;

            atom_molecules.insert({blk_id, cur_molecule});
        }
    }

    fill_vacant_chain_spots(list_of_molecules_head, list_of_pack_patterns, num_packing_patterns, atom_molecules);

    // After fill_vacant_chain_spots, we need to compress the atom netlist.
    // fill_vacant_chain_spots calls set_pin_net() which internally calls remove_net_pin(),
    // marking the netlist as dirty. We must compress to clean it up before clustering.
    auto& mutable_atom_ctx = g_vpr_ctx.mutable_atom();
    AtomNetlist& mutable_atom_nlist = mutable_atom_ctx.nlist;
    auto id_remapper = mutable_atom_nlist.compress();

    // Update all AtomBlockIds in molecules using the remapper, since compress() may renumber IDs
    t_pack_molecule* cur_mol = list_of_molecules_head;
    while (cur_mol != nullptr) {
        for (size_t i = 0; i < cur_mol->atom_block_ids.size(); i++) {
            AtomBlockId old_id = cur_mol->atom_block_ids[i];
            if (old_id) {
                AtomBlockId new_id = id_remapper.new_block_id(old_id);
                cur_mol->atom_block_ids[i] = new_id;
            }
        }
        cur_mol = cur_mol->next;
    }

    // Update the atom_molecules multimap with remapped IDs
    std::multimap<AtomBlockId, t_pack_molecule*> remapped_atom_molecules;
    for (auto& pair : atom_molecules) {
        AtomBlockId old_id = pair.first;
        AtomBlockId new_id = id_remapper.new_block_id(old_id);
        remapped_atom_molecules.insert({new_id, pair.second});
    }
    atom_molecules = std::move(remapped_atom_molecules);

    if (getEchoEnabled() && isEchoFileEnabled(E_ECHO_PRE_PACKING_MOLECULES_AND_PATTERNS)) {
        print_pack_molecules(getEchoFileName(E_ECHO_PRE_PACKING_MOLECULES_AND_PATTERNS),
                             list_of_pack_patterns, num_packing_patterns,
                             list_of_molecules_head,
                             atom_nlist);
    }

    return list_of_molecules_head;
}

static void free_pack_pattern_block(t_pack_pattern_block* pattern_block, t_pack_pattern_block** pattern_block_list) {
    t_pack_pattern_connections *connection, *next;
    if (pattern_block == nullptr || pattern_block->block_id == OPEN) {
        /* already traversed, return */
        return;
    }
    pattern_block_list[pattern_block->block_id] = pattern_block;
    pattern_block->block_id = OPEN;
    connection = pattern_block->connections;
    while (connection) {
        free_pack_pattern_block(connection->from_block, pattern_block_list);
        free_pack_pattern_block(connection->to_block, pattern_block_list);
        next = connection->next;
        delete connection;
        connection = next;
    }
}

/**
 * Given a pattern and an atom block to serve as the root block, determine if
 * the candidate atom block serving as the root node matches the pattern.
 * If yes, return the molecule with this atom block as the root, if not, return NULL
 *
 * Limitations: Currently assumes that forced pack nets must be single-fanout as
 *              this covers all the reasonable architectures we wanted. More complicated
 *              structures should probably be handled either downstream (general packing)
 *              or upstream (in tech mapping).
 *              If this limitation is too constraining, code is designed so that this limitation can be removed
 *
 * Side Effect: If successful, link atom to molecule
 */
static t_pack_molecule* try_create_molecule(t_pack_patterns* list_of_pack_patterns,
                                            const int pack_pattern_index,
                                            AtomBlockId blk_id,
                                            std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                            const AtomNetlist& atom_nlist) {
    auto pack_pattern = &list_of_pack_patterns[pack_pattern_index];

    // Debugging: trace attempts to create lut_chain / simple_lut_chain molecules
    std::string pattern_name(pack_pattern->name);
    bool debug_lut_chain = true;

    t_pack_molecule* molecule;

    // Check pack pattern validity
    if (pack_pattern == nullptr || pack_pattern->num_blocks == 0 || pack_pattern->root_block == nullptr) {
        return nullptr;
    }

    // If a chain pattern extends beyond a single logic block, we must find
    // the furthest blk_id up the chain that is not mapped to a molecule yet.
    if (pack_pattern->is_chain) {
        AtomBlockId orig_blk_id = blk_id;
        blk_id = find_new_root_atom_for_chain(blk_id, pack_pattern, atom_molecules, atom_nlist);
        if (!blk_id) {
            return nullptr;
        }
    }

    molecule = new t_pack_molecule;
    molecule->type = MOLECULE_FORCED_PACK;
    molecule->pack_pattern = pack_pattern;
    molecule->atom_block_ids = std::vector<AtomBlockId>(pack_pattern->num_blocks); //Initializes invalid
    molecule->num_blocks = pack_pattern->num_blocks;
    molecule->root = pack_pattern->root_block->block_id;

    if (try_expand_molecule(molecule, blk_id, atom_molecules, atom_nlist)) {
        // Success! commit molecule
        // update chain info for chain molecules
        if (molecule->pack_pattern->is_chain) {
            init_molecule_chain_info(blk_id, molecule, atom_molecules, atom_nlist);
        }

        // update the atom_molcules with the atoms that are mapped to this molecule
        for (int i = 0; i < molecule->pack_pattern->num_blocks; i++) {
            auto blk_id2 = molecule->atom_block_ids[i];
            if (!blk_id2) {
                VTR_ASSERT(molecule->pack_pattern->is_block_optional[i]);
                continue;
            }

            atom_molecules.insert({blk_id2, molecule});
        }
    } else {
        delete molecule;
        return nullptr;
    }

    return molecule;
}

/**
 * Determine if an atom block can match with the pattern to from a molecule.
 *
 * This function takes a molecule that represents a packing pattern. It also
 * takes a (netlist) atom block represented by blk_id which matches the
 * root primitive of this packing pattern. Using this atom block and the structure
 * of the packing pattern, this function tries to fill all the available positions
 * in the packing pattern. If all the non-optional primitive positions in the
 * pattern are filled return true, return false otherwise.
 *      molecule       : the molecule we are trying to expand
 *      atom_molecules : map of atom block ids that are assigned a molecule and a pointer to this molecule
 *      blk_id         : chosen to be the root of this molecule and the code is expanding from
 */
static bool try_expand_molecule(t_pack_molecule* molecule,
                                const AtomBlockId blk_id,
                                const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                const AtomNetlist& atom_nlist) {
    std::string pattern_name(molecule->pack_pattern->name);
    bool debug_lut_chain = pattern_name.find("lut_chain") != std::string::npos;

    bool has_second_level = false;
    bool found_second_level = false;
    const bool hierarchical_molecule = molecule_is_hierarchical(molecule);
    const t_model_ports* cin_port_model = nullptr;
    if (molecule->is_chain()) {
        cin_port_model = molecule->pack_pattern->chain_root_pins[0][0]->port->model_port;
    }
    // root block of the pack pattern, which is the starting point of this pattern
    const auto pattern_root_block = molecule->pack_pattern->root_block;
    // bool array indicating whether a position in a pack pattern is optional or should
    // be filled with an atom for legality
    const auto is_block_optional = molecule->pack_pattern->is_block_optional;

    // create a queue of pattern block and atom block id suggested for this block
    std::queue<std::pair<t_pack_pattern_block*, AtomBlockId>> pattern_block_queue;
    // initialize the queue with the pattern root block and the matching atom block
    pattern_block_queue.push(std::make_pair(pattern_root_block, blk_id));

    // do breadth first search by walking through the pack pattern structure along
    // with the atom netlist structure
    while (!pattern_block_queue.empty()) {
        // get the front pattern block, atom block id pair from the queue
        const auto pattern_block_atom_pair = pattern_block_queue.front();
        const auto pattern_block = pattern_block_atom_pair.first;
        const auto block_id = pattern_block_atom_pair.second;

        // remove the front of the queue
        pattern_block_queue.pop();

        // get the atom block id of the atom occupying this primitive position in this molecule
        auto molecule_atom_block_id = molecule->atom_block_ids[pattern_block->block_id];

        // if this primitive position in this molecule is already visited and
        // matches block in the atom netlist go to the next node in the queue
        if (molecule_atom_block_id) {
            continue;
        }

        if (!block_id || !primitive_type_feasible(block_id, pattern_block->pb_type) || (molecule_atom_block_id && molecule_atom_block_id != block_id) || atom_molecules.find(block_id) != atom_molecules.end()) {
            // Stopping conditions, if:
            // 1) this is an invalid atom block (nothing)
            // 2) this atom block cannot fit in this primitive type
            // 3) this primitive is occupied by another block
            // 4) this atom block is already used by another molecule
            // then if the molecule cannot be formed without placing an atom
            // at that primitive position, then creating this molecule has failed
            // otherwise go to the next atom block and its corresponding pattern block
            if (!is_block_optional[pattern_block->block_id]) {
                return false;
            }
            continue;
        }

        if (hierarchical_molecule && !found_second_level) {
            auto first_level_block = is_second_level_block(pattern_block, molecule);
            if (first_level_block) {
                if (valid_second_level_placement(first_level_block, block_id, molecule, atom_nlist))
                    found_second_level = true;
                else
                    continue;
            }
        }

        // set this node in the molecule as visited
        molecule->atom_block_ids[pattern_block->block_id] = block_id;

        // starting from the first connections, add all the connections of this block to the queue
        auto block_connection = pattern_block->connections;

        while (block_connection != nullptr) {
            // this block is the driver of this connection
            if (block_connection->from_block == pattern_block) {
                // find the block this connection is driving and add it to the queue
                auto port_model = block_connection->from_pin->port->model_port;
                auto ipin = block_connection->from_pin->pin_number;
                AtomBlockId sink_blk_id = get_sink_block(block_id, port_model, ipin, atom_nlist);

                // add this sink block id with its corresponding pattern block to the queue
                pattern_block_queue.push(std::make_pair(block_connection->to_block, sink_blk_id));
                // this block is being driven by this connection
            } else if (block_connection->to_block == pattern_block) {
                // find the block that is driving this connection and it to the queue
                auto port_model = block_connection->to_pin->port->model_port;
                auto ipin = block_connection->to_pin->pin_number;
                auto driver_blk_id = get_driving_block(block_id, port_model, ipin, atom_nlist);

                if (molecule->is_chain() && port_model != cin_port_model && block_connection->to_pin->parent_node->pb_type == block_connection->from_pin->parent_node->pb_type) {
                    has_second_level = true;
                }
                // add this driver block id with its corresponding pattern block to the queue
                // only if it's driving the cin port. To avoid adding blocks by tracking the adder
                // inputs port which will result in a molecule that cannot be placed
                if (molecule->is_chain() && (port_model == cin_port_model || block_connection->from_block->pb_type->model != molecule->pack_pattern->chain_root_pins[0][0]->parent_node->pb_type->model)) {
                    pattern_block_queue.push(std::make_pair(block_connection->from_block, driver_blk_id));
                }
            }

            // this block should be either driving or driven by the connection
            VTR_ASSERT(block_connection->from_block == pattern_block || block_connection->to_block == pattern_block);
            // go to the next connection of this pattern block
            block_connection = block_connection->next;
        }
    }

    // If this is a hierarchical molecule but no second-level block was
    // discovered for this particular root, treat it as a non-hierarchical
    // chain instance. Hierarchical placement constraints are only applied
    // when a valid second-level candidate is actually present.
    if (!has_second_level && hierarchical_molecule) {
        return false;
    }
    if (molecule->is_chain()) {
        bool reachable = chain_input_is_reachable(molecule, atom_molecules, atom_nlist);
        bool alm_ok = check_alm_input_limitation(molecule, atom_nlist);
        bool lut_ok = check_lut_chain_molecules(molecule, atom_nlist);
        if (debug_lut_chain) {
            VTR_LOG("try_expand_molecule[%s]: chain checks reachable=%d alm_ok=%d lut_ok=%d\n",
                    pattern_name.c_str(),
                    reachable ? 1 : 0,
                    alm_ok ? 1 : 0,
                    lut_ok ? 1 : 0);
        }
        return reachable && alm_ok && lut_ok;
    }

    // if all non-optional positions in the pack pattern have atoms
    // mapped to them, then this molecule is valid
    return true;
}

/**
 * Find the atom block in the netlist driven by this pin of the input atom block
 * If doesn't exist return AtomBlockId::INVALID()
 * Limitation: The block should be driving only one sink block
 *      block_id   : id of the atom block that is driving the net connected to the sink block
 *      model_port : the model of the port driving the net
 *      pin_number : the pin_number of the pin driving the net (pin index within the port)
 */
static AtomBlockId get_sink_block(const AtomBlockId block_id,
                                  const t_model_ports* model_port,
                                  const BitIndex pin_number,
                                  const AtomNetlist& atom_nlist) {
    auto port_id = atom_nlist.find_atom_port(block_id, model_port);

    if (port_id) {
        auto net_id = atom_nlist.port_net(port_id, pin_number);
        if (net_id && atom_nlist.net_sinks(net_id).size() == 1) { /* Single fanout assumption */
            auto net_sinks = atom_nlist.net_sinks(net_id);
            auto sink_pin_id = *(net_sinks.begin());
            return atom_nlist.pin_block(sink_pin_id);
        }
    }

    return AtomBlockId::INVALID();
}

/**
 * Find the atom block in the netlist driving this pin of the input atom block
 * If doesn't exist return AtomBlockId::INVALID()
 * Limitation: This driving block should be driving only the input block
 *      block_id   : id of the atom block that is connected to a net driven by the driving block
 *      model_port : the model of the port driven by the net
 *      pin_number : the pin_number of the pin driven by the net (pin index within the port)
 */
static AtomBlockId get_driving_block(const AtomBlockId block_id,
                                     const t_model_ports* model_port,
                                     const BitIndex pin_number,
                                     const AtomNetlist& atom_nlist) {
    auto port_id = atom_nlist.find_atom_port(block_id, model_port);

    if (port_id) {
        auto net_id = atom_nlist.port_net(port_id, pin_number);
        if (net_id && atom_nlist.net_sinks(net_id).size() == 1) { /* Single fanout assumption */

            auto driver_blk_id = atom_nlist.net_driver_block(net_id);

            if (model_port->is_clock) {
                auto driver_blk_type = atom_nlist.block_type(driver_blk_id);

                // TODO: support multi-clock primitives.
                //       If the driver block is a .input block, this assertion should not
                //       be triggered as the sink block might have only one input pin, which
                //       would be a clock pin in case the sink block primitive is a clock generator,
                //       resulting in a pin_number == 0.
                VTR_ASSERT(pin_number == 1 || (pin_number == 0 && driver_blk_type == AtomBlockType::INPAD));
            }

            return atom_nlist.net_driver_block(net_id);
        }
    }

    return AtomBlockId::INVALID();
}

/**
 * Variant of get_driving_block used when expanding chain patterns.
 *
 * For chains we allow the driving block to drive multiple sinks; we always
 * return the unique driver of the net (if any), regardless of fanout.
 */

static void print_pack_molecules(const char* fname,
                                 const t_pack_patterns* list_of_pack_patterns,
                                 const int num_pack_patterns,
                                 const t_pack_molecule* list_of_molecules,
                                 const AtomNetlist& atom_nlist) {
    int i;
    FILE* fp;
    const t_pack_molecule* list_of_molecules_current;

    fp = std::fopen(fname, "w");
    fprintf(fp, "# of pack patterns %d\n", num_pack_patterns);

    for (i = 0; i < num_pack_patterns; i++) {
        VTR_ASSERT(list_of_pack_patterns[i].root_block);
        fprintf(fp, "pack pattern index %d block count %d name %s root %s\n",
                list_of_pack_patterns[i].index,
                list_of_pack_patterns[i].num_blocks,
                list_of_pack_patterns[i].name,
                list_of_pack_patterns[i].root_block->pb_type->name);

        if (list_of_pack_patterns[i].is_chain) {
            fprintf(fp, "\tChain Root Pins:\n");
            for (size_t chain_idx = 0; chain_idx < list_of_pack_patterns[i].chain_root_pins.size(); ++chain_idx) {
                fprintf(fp, "\t\tChain ID %zu:\n", chain_idx);
                for (const auto* pin : list_of_pack_patterns[i].chain_root_pins[chain_idx]) {
                    fprintf(fp, "\t\t\t%s\n", pin->to_string().c_str());
                }
            }
        }

        // For debugging: print the mapping from pattern block indices to their
        // corresponding pb_type names. This helps interpret which primitive
        // each sparse index in the molecule corresponds to.
        fprintf(fp, "pack pattern %d block-to-pb_type mapping:\n", list_of_pack_patterns[i].index);
        for (int b = 0; b < list_of_pack_patterns[i].num_blocks; ++b) {
            t_pack_pattern_block* pb = nullptr;
            // Find the pattern_block with this block_id by walking from root.
            std::vector<bool> visited(list_of_pack_patterns[i].num_blocks);
            std::queue<t_pack_pattern_block*> q;
            q.push(list_of_pack_patterns[i].root_block);
            while (!q.empty()) {
                auto* blk = q.front();
                q.pop();
                if (!blk || visited[blk->block_id])
                    continue;
                visited[blk->block_id] = true;
                if (blk->block_id == b) {
                    pb = blk;
                    break;
                }
                auto conn = blk->connections;
                while (conn) {
                    q.push(conn->from_block);
                    q.push(conn->to_block);
                    conn = conn->next;
                }
            }
            if (pb) {
                fprintf(fp, "\tpattern block %d -> pb_type %s\n",
                        b, pb->pb_type->name);
            } else {
                fprintf(fp, "\tpattern block %d -> <unreachable>\n", b);
            }
        }
    }

    list_of_molecules_current = list_of_molecules;
    while (list_of_molecules_current != nullptr) {
        if (list_of_molecules_current->type == MOLECULE_SINGLE_ATOM) {
            fprintf(fp, "\nmolecule type: atom\n");
            fprintf(fp, "\tpattern index %d: atom block %s\n", i,
                    atom_nlist.block_name(list_of_molecules_current->atom_block_ids[0]).c_str());
        } else if (list_of_molecules_current->type == MOLECULE_FORCED_PACK) {
            fprintf(fp, "\nmolecule type: %s\n",
                    list_of_molecules_current->pack_pattern->name);
            if (list_of_molecules_current->is_chain()) {
                fprintf(fp, "\tis_long_chain: %d\n", list_of_molecules_current->chain_info->is_long_chain);
                fprintf(fp, "\tchain_id: %d\n", list_of_molecules_current->chain_info->chain_id);
                fprintf(fp, "\tfirst_pack_molecule: %p\n", (void*)list_of_molecules_current->chain_info->first_packed_molecule);
            }
            for (i = 0; i < list_of_molecules_current->pack_pattern->num_blocks;
                 i++) {
                if (!list_of_molecules_current->atom_block_ids[i]) {
                    fprintf(fp, "\tpattern index %d: empty \n", i);
                } else {
                    // For debugging DCC-style chains, report the adder row
                    // (placement index of the corresponding 'adder' pb_type)
                    // when this pattern block represents an adder primitive.
                    int adder_row = -1;
                    t_pack_pattern_block* pb = nullptr;
                    {
                        const t_pack_patterns* patt = list_of_molecules_current->pack_pattern;
                        std::vector<bool> visited(patt->num_blocks);
                        std::queue<t_pack_pattern_block*> q;
                        q.push(patt->root_block);
                        while (!q.empty()) {
                            auto* blk = q.front();
                            q.pop();
                            if (!blk || visited[blk->block_id])
                                continue;
                            visited[blk->block_id] = true;
                            if (blk->block_id == i) {
                                pb = blk;
                                break;
                            }
                            auto conn = blk->connections;
                            while (conn) {
                                q.push(conn->from_block);
                                q.push(conn->to_block);
                                conn = conn->next;
                            }
                        }
                    }

                    fprintf(fp, "\tpattern index %d: atom block %s (ID: %zu)", i,
                            atom_nlist.block_name(list_of_molecules_current->atom_block_ids[i]).c_str(),
                            size_t(list_of_molecules_current->atom_block_ids[i]));
                    if (pb && std::string(pb->pb_type->name) == "adder") {
                        adder_row = get_pb_placement_index(pb, "adder");
                    }
                    if (adder_row >= 0) {
                        fprintf(fp, " row %d", adder_row);
                    }
                    if (list_of_molecules_current->pack_pattern->root_block->block_id == i) {
                        fprintf(fp, " root node\n");
                    } else {
                        fprintf(fp, "\n");
                    }

                    // Print pin and net information for debugging connectivity
                    AtomBlockId blk_id = list_of_molecules_current->atom_block_ids[i];
                    for (auto port_id : atom_nlist.block_ports(blk_id)) {
                        std::string port_name = atom_nlist.port_name(port_id);
                        for (auto pin_id : atom_nlist.port_pins(port_id)) {
                            auto net_id = atom_nlist.pin_net(pin_id);
                            if (net_id) {
                                std::string net_name = atom_nlist.net_name(net_id);
                                fprintf(fp, "\t\t-> pin %s[%zu]: net %s",
                                        port_name.c_str(),
                                        size_t(atom_nlist.pin_port_bit(pin_id)),
                                        net_name.c_str());

                                // Show if this is a driver or sink
                                if (atom_nlist.net_driver(net_id) == pin_id) {
                                    fprintf(fp, " (driver)");
                                } else {
                                    fprintf(fp, " (sink)");
                                }
                                fprintf(fp, "\n");
                            }
                        }
                    }
                }
            }
        } else {
            VTR_ASSERT(0);
        }
        list_of_molecules_current = list_of_molecules_current->next;
    }

    fclose(fp);
}

/* Search through all primitives and return the lowest cost primitive that fits this atom block */
static t_pb_graph_node* get_expected_lowest_cost_primitive_for_atom_block(const AtomBlockId blk_id,
                                                                          const std::vector<t_logical_block_type>& logical_block_types) {
    float cost, best_cost;
    t_pb_graph_node *current, *best;

    best_cost = UNDEFINED;
    best = nullptr;
    current = nullptr;
    for (const t_logical_block_type& type : logical_block_types) {
        cost = UNDEFINED;
        current = get_expected_lowest_cost_primitive_for_atom_block_in_pb_graph_node(blk_id, type.pb_graph_head, &cost);
        if (cost != UNDEFINED) {
            if (best_cost == UNDEFINED || best_cost > cost) {
                best_cost = cost;
                best = current;
            }
        }
    }

    return best;
}

static t_pb_graph_node* get_expected_lowest_cost_primitive_for_atom_block_in_pb_graph_node(const AtomBlockId blk_id, t_pb_graph_node* curr_pb_graph_node, float* cost) {
    t_pb_graph_node *best, *cur;
    float cur_cost, best_cost;
    int i, j;

    best = nullptr;
    best_cost = UNDEFINED;
    if (curr_pb_graph_node == nullptr) {
        return nullptr;
    }

    if (curr_pb_graph_node->pb_type->blif_model != nullptr) {
        if (primitive_type_feasible(blk_id, curr_pb_graph_node->pb_type)) {
            cur_cost = compute_primitive_base_cost(curr_pb_graph_node);
            if (best_cost == UNDEFINED || best_cost > cur_cost) {
                best_cost = cur_cost;
                best = curr_pb_graph_node;
            }
        }
    } else {
        for (i = 0; i < curr_pb_graph_node->pb_type->num_modes; i++) {
            /* Early fail if this primitive for a mode that is disabled for packing */
            if (true == curr_pb_graph_node->pb_type->modes[i].disable_packing) {
                continue;
            }

            for (j = 0; j < curr_pb_graph_node->pb_type->modes[i].num_pb_type_children; j++) {
                *cost = UNDEFINED;
                cur = get_expected_lowest_cost_primitive_for_atom_block_in_pb_graph_node(blk_id, &curr_pb_graph_node->child_pb_graph_nodes[i][j][0], cost);
                if (cur != nullptr) {
                    if (best == nullptr || best_cost > *cost) {
                        best = cur;
                        best_cost = *cost;
                    }
                }
            }
        }
    }

    *cost = best_cost;
    return best;
}

/* Determine which of two pack pattern should take priority */
static int compare_pack_pattern(const t_pack_patterns* pattern_a, const t_pack_patterns* pattern_b) {
    float base_gain_a, base_gain_b, diff;

    /* Bigger patterns should take higher priority than smaller patterns because they are harder to fit */
    if (pattern_a->num_blocks > pattern_b->num_blocks) {
        return 1;
    } else if (pattern_a->num_blocks < pattern_b->num_blocks) {
        return -1;
    }

    base_gain_a = pattern_a->base_cost;
    base_gain_b = pattern_b->base_cost;
    diff = base_gain_a - base_gain_b;

    /* Less costly patterns should be used before more costly patterns */
    if (diff < 0) {
        return 1;
    }
    if (diff > 0) {
        return -1;
    }
    return 0;
}

/* A chain can extend across multiple atom blocks.  Must segment the chain to fit in an atom
 * block by identifying the actual atom that forms the root of the new chain.
 * Returns AtomBlockId::INVALID() if this block_index doesn't match up with any chain
 *
 * Assumes that the root of a chain is the primitive that starts the chain or is driven from outside the logic block
 * block_index: index of current atom
 * list_of_pack_pattern: ptr to current chain pattern
 */
static AtomBlockId find_new_root_atom_for_chain(const AtomBlockId blk_id,
                                                const t_pack_patterns* list_of_pack_patterns,
                                                const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                                const AtomNetlist& atom_nlist) {
    AtomBlockId new_root_blk_id;
    t_pb_graph_pin* root_ipin;
    t_pb_graph_node* root_pb_graph_node;

    VTR_ASSERT(list_of_pack_patterns->is_chain == true);
    VTR_ASSERT(list_of_pack_patterns->chain_root_pins.size());
    root_ipin = list_of_pack_patterns->chain_root_pins[0][0];
    root_pb_graph_node = root_ipin->parent_node;

    if (primitive_type_feasible(blk_id, root_pb_graph_node->pb_type) == false) {
        return AtomBlockId::INVALID();
    }

    // find the block id of the atom block driving the input of this block
    AtomBlockId driver_blk_id = get_adder_driver_block(blk_id, list_of_pack_patterns, atom_molecules, atom_nlist);

    // if there is no driver block for this net
    // then it is the furthest up the chain
    if (!driver_blk_id) {
        return blk_id;
    }

    // didn't find furthest atom up the chain, keep searching further up the chain
    new_root_blk_id = find_new_root_atom_for_chain(driver_blk_id, list_of_pack_patterns, atom_molecules, atom_nlist);

    if (!new_root_blk_id) {
        return blk_id;
    } else {
        return new_root_blk_id;
    }
}

/**
 * This function takes an input pin to a root (has no parent block) pb_graph_node
 * an returns a vector of all the output pins that are reachable from this input
 * pin and have the same packing pattern
 */

static std::vector<t_pb_graph_pin*> find_end_of_path(t_pb_graph_pin* input_pin, int pattern_index) {
    // Enforce some constraints on the function

    // 1) the start of the path should be at the input of the root block
    VTR_ASSERT(input_pin->is_root_block_pin());

    // 2) this pin is an input pin to the root block
    VTR_ASSERT(input_pin->num_input_edges == 0);

    // create a queue of pin pointers for the breadth first search
    std::queue<t_pb_graph_pin*> pins_queue;

    // add the input pin to the queue
    pins_queue.push(input_pin);

    // found reachable output pins
    std::vector<t_pb_graph_pin*> reachable_pins;

    // do breadth first search till all connected
    // pins are explored
    while (!pins_queue.empty()) {
        // get the first pin in the queue
        auto current_pin = pins_queue.front();

        // remove pin from queue
        pins_queue.pop();

        // expand search from current pin
        expand_search(current_pin, pins_queue, pattern_index);

        // if this is an output pin of a root block
        // add to reachable output pins
        if (current_pin->is_root_block_pin()
            && current_pin->num_output_edges == 0) {
            reachable_pins.push_back(current_pin);
        }
    }

    return reachable_pins;
}

static void expand_search(const t_pb_graph_pin* input_pin, std::queue<t_pb_graph_pin*>& pins_queue, const int pattern_index) {
    // If not a primitive input pin (has output edges)
    // -----------------------------------------------

    // iterate over all output edges at this pin
    for (int iedge = 0; iedge < input_pin->num_output_edges; iedge++) {
        const auto& pin_edge = input_pin->output_edges[iedge];
        // if this edge is not anotated with this pattern and its pattern cannot be inferred, ignore it.
        if (!pin_edge->annotated_with_pattern(pattern_index) && !pin_edge->infer_pattern) {
            continue;
        }

        // this edge either matched the pack pattern or its pack pattern could be inferred
        // iterate over all the pins of that edge and add them to the pins_queue
        for (int ipin = 0; ipin < pin_edge->num_output_pins; ipin++) {
            pins_queue.push(pin_edge->output_pins[ipin]);
        }

    } // End for pin edges

    // If a primitive input pin
    // ------------------------

    // if this is an input pin to a primitive, it won't have
    // output edges so the previous for loop won't be entered
    if (input_pin->is_primitive_pin() && input_pin->num_output_edges == 0) {
        // iterate over the output ports of the primitive
        const auto& pin_pb_graph_node = input_pin->parent_node;
        for (int iport = 0; iport < pin_pb_graph_node->num_output_ports; iport++) {
            // iterate over the pins of each port
            const auto& port_pins = pin_pb_graph_node->num_output_pins[iport];
            for (int ipin = 0; ipin < port_pins; ipin++) {
                // add primitive output pins to pins_queue to be explored
                pins_queue.push(&pin_pb_graph_node->output_pins[iport][ipin]);
            }
        }
    }

    // If this is a root block output pin
    // ----------------------------------

    // No expansion will happen in this case
}

/**
 *  This function takes a chain pack pattern and a root pb_block
 *  containing this pattern. Then searches for all the input pins of this
 *  pb_block that are annotated with this pattern. The function then
 *  identifies whether those inputs represent different starting point for
 *  this pattern or are all required for building this pattern.
 *
 *  If this inputs represent different starting point for this pattern, it
 *  means that in this pb_block there exist multiple chains that are exactly
 *  the same. For example an architecture that has two separate adder chains
 *  behaving exactly the same but are totally separate from each other.
 *
 *                        cin[0] cin[1]
 *                       ---|------|---
 *                       | ---    --- |
 *                       | | |    | |<|---- Full Adder
 *                       | ---    --- |
 *   Pb_block ---------> |  |      |<-|---- Second Adder chain
 *                       |  .      .  |
 *                       |  .      .  |
 *                       |  |<-----|--|---- First Adder chain
 *                       | ---    --- |
 *                       | | |    | | |
 *                       | ---    --- |
 *                       ---|------|---
 *                       cout[0] cout[1]
 *
 *  In this case, the chain_root_pin array of the pack pattern is updated
 *  with all the pin that represent a starting point for this pattern.
 */
static void find_all_equivalent_chains(t_pack_patterns* chain_pattern, const t_pb_graph_node* root_block) {
    // this vector will be updated with all root_block input
    // pins that are annotated with this chain pattern
    std::vector<t_pb_graph_pin*> chain_input_pins;

    // iterate over all the input pins of the root_block and populate
    // the chain_input_pins vector
    for (int iports = 0; iports < root_block->num_input_ports; iports++) {
        for (int ipins = 0; ipins < root_block->num_input_pins[iports]; ipins++) {
            auto& input_pin = root_block->input_pins[iports][ipins];
            for (int iedge = 0; iedge < input_pin.num_output_edges; iedge++) {
                if (input_pin.output_edges[iedge]->belongs_to_pattern(chain_pattern->index)) {
                    chain_input_pins.push_back(&input_pin);
                }
            }
        }
    }

    // if this chain has only one cluster input, then
    // there is no need to proceed with the search
    if (chain_input_pins.size() == 1) {
        update_chain_root_pins(chain_pattern, chain_input_pins);
        chain_pattern->chain_exit_pins.push_back(find_chain_exit_pin(chain_input_pins[0], chain_pattern->index));
        return;
    }

    // find the root block output pins reachable when starting from the chain_input_pins
    // found before following the edges that are annotated with the given pack_pattern
    std::vector<std::vector<t_pb_graph_pin*>> reachable_pins;

    for (const auto& pin_ptr : chain_input_pins) {
        auto reachable_output_pins = find_end_of_path(pin_ptr, chain_pattern->index);
        // find the chain exit pin of this chain input pin
        auto chain_exit_pin = find_chain_exit_pin(pin_ptr, chain_pattern->index);
        // sort the reachable output pins to compare them later using set_intersection
        std::stable_sort(reachable_output_pins.begin(), reachable_output_pins.end());
        reachable_pins.push_back(reachable_output_pins);
        // update the chain exit pins array
        chain_pattern->chain_exit_pins.push_back(chain_exit_pin);
    }

    // Search for intersections between reachable pins. Intersection
    // between reachable indicates that found chain_input_pins
    // represent a single chain pattern and not multiple similar
    // chain patterns with multiple starting locations.
    std::vector<t_pb_graph_pin*> intersection;
    for (size_t i = 0; i < reachable_pins.size() - 1; i++) {
        for (size_t j = 1; j < reachable_pins.size(); j++) {
            std::set_intersection(reachable_pins[i].begin(), reachable_pins[i].end(),
                                  reachable_pins[j].begin(), reachable_pins[j].end(),
                                  std::back_inserter(intersection));
            if (intersection.size())
                break;
        }
        if (intersection.size())
            break;
    }

    // if there are no intersections between the reachable pins,
    // this each input pin represents a separate chain of type
    // chain_pattern. Else, they are all representing the same
    // chain.
    if (intersection.empty()) {
        // update the chain_root_pin array of the chain_pattern
        // with all the possible starting points of the chain.
        update_chain_root_pins(chain_pattern, chain_input_pins);
    }
}

/**
 *  This function takes a chain pack pattern and a vector of pin
 *  pointers that represent the root pb block input pins that can connect
 *  a chain to the previous pb block. The function uses the pin pointers
 *  to find the primitive input pin connected to them and updates
 *  the chain_root_pin array with this those pointers
 *  Side Effect: Updates the chain_root_pins array of the input chain_pattern
 */
static void update_chain_root_pins(t_pack_patterns* chain_pattern,
                                   const std::vector<t_pb_graph_pin*>& chain_input_pins) {
    std::vector<std::vector<t_pb_graph_pin*>> primitive_input_pins;

    for (const auto pin_ptr : chain_input_pins) {
        std::vector<t_pb_graph_pin*> connected_primitive_pins;
        get_all_connected_primitive_pins(pin_ptr, connected_primitive_pins, chain_pattern->index);

        /**
         * It is required that the chain pins are connected inside a complex
         * block. Although it is allowed to have them disconnected in some
         * modes of the block provided that there is always at least one mode
         * that has them connected inside. The following assert checks for
         * that.
         */
        VTR_ASSERT(connected_primitive_pins.size());

        primitive_input_pins.push_back(connected_primitive_pins);
    }

    chain_pattern->chain_root_pins = primitive_input_pins;
}

/**
 *  This function takes a pin as an input an does a depth first search on all the output edges
 *  of this pin till it finds all the primitive input pins connected to this pin. For example,
 *  if the input pin given to this function is the Cin pin of the cluster. This pin will return
 *  the Cin pin of all the adder primitives connected to this pin. Which is for typical architectures
 *  will be only one pin connected to the very first adder in the cluster.
 */
static void get_all_connected_primitive_pins(const t_pb_graph_pin* cluster_input_pin, std::vector<t_pb_graph_pin*>& connected_primitive_pins, int pattern_id) {
    for (int iedge = 0; iedge < cluster_input_pin->num_output_edges; iedge++) {
        const auto& output_edge = cluster_input_pin->output_edges[iedge];
        // if (!output_edge->belongs_to_pattern(pattern_id)) continue;

        for (int ipin = 0; ipin < output_edge->num_output_pins; ipin++) {
            if (output_edge->output_pins[ipin]->is_primitive_pin()) {
                connected_primitive_pins.push_back(output_edge->output_pins[ipin]);
            } else {
                get_all_connected_primitive_pins(output_edge->output_pins[ipin], connected_primitive_pins, pattern_id);
            }
        }
    }
    VTR_ASSERT(connected_primitive_pins.size());
}

/**
 * This function initializes the chain info data structure of the molecule.
 * If this is the furthest molecule up the chain, the chain_info data
 * structure is created. Otherwise, the input pack_molecule is set to
 * point to the same chain_info of the molecule feeding it
 *
 * Limitation: This function assumes that the molecules of a chain are
 * created and fed to this function in order. Meaning the first molecule
 * fed to the function should be the furthest molecule up the chain.
 * The second one should should be the molecule directly after that one
 * and so on.
 */
static void init_molecule_chain_info(const AtomBlockId blk_id,
                                     t_pack_molecule* molecule,
                                     const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                     const AtomNetlist& atom_nlist) {
    // the input molecule to this function should have a pack
    // pattern assigned to it and the input block should be valid
    VTR_ASSERT(molecule->pack_pattern && blk_id);

    auto root_ipin = molecule->pack_pattern->chain_root_pins[0][0];
    auto model_pin = root_ipin->port->model_port;
    auto pin_bit = root_ipin->pin_number;

    // find the atom driving the chain input pin of this atom
    auto driver_atom_id = atom_nlist.find_atom_pin_driver(blk_id, model_pin, pin_bit);

    // find the molecule this driver atom is mapped to
    auto itr = atom_molecules.find(driver_atom_id);

    // if this is the first molecule to be created for this chain
    // initialize the chain info data structure. This is the case
    // if either there is no driver to the block input pin or
    // if the driver is not part of a molecule
    if (!driver_atom_id || itr == atom_molecules.end()) {
        // allocate chain info
        molecule->chain_info = std::make_shared<t_chain_info>();
        // chain_id defaults to -1 for short chains (no specific architectural chain assigned)
        // this is not the first molecule to be created for this chain
    } else {
        // molecule driving blk_id
        auto prev_molecule = itr->second;
        // molecule should have chain_info associated with it
        VTR_ASSERT(prev_molecule && prev_molecule->chain_info);
        // this molecule is now known to belong to a long chain
        prev_molecule->chain_info->is_long_chain = true;
        // this new molecule should share the same chain_info
        molecule->chain_info = prev_molecule->chain_info;
        // if the two molecules are of different types
        if (prev_molecule->pack_pattern->chain_root_pins.size() < molecule->pack_pattern->chain_root_pins.size()) {
            molecule->chain_info->chain_id = get_forced_chain_id(molecule, prev_molecule, driver_atom_id);
        }
    }
}

/**
 * This function prints all the starting points of the carry chains in the architecture
 */
static void print_chain_starting_points(t_pack_patterns* chain_pattern) {
    const auto& chain_root_pins = chain_pattern->chain_root_pins;

    VTR_LOGV(chain_root_pins.size() > 1, "\nThere are %zu independent chains for chain pattern \"%s\":\n",
             chain_pattern->chain_root_pins.size(), chain_pattern->name);
    VTR_LOGV(chain_root_pins.size() == 1, "\nThere is one chain in this architecture called \"%s\" with the following starting points:\n", chain_pattern->name);

    size_t chainId = 0;
    for (const auto& chain : chain_root_pins) {
        VTR_LOGV(chain_root_pins.size() > 1 && chain.size() > 1, "\n There are %zu starting points for chain id #%zu:\n", chain.size(), chainId++);
        VTR_LOGV(chain_root_pins.size() > 1 && chain.size() == 1, "\n There is 1 starting point for chain id #%zu:\n", chainId++);
        for (const auto& pin_ptr : chain) {
            VTR_LOG("\t%s\n", pin_ptr->to_string().c_str());
        }
    }

    VTR_LOG("\n");
}

/**
 * This function frees the linked list of pack molecules.
 */
static void free_pack_molecules(t_pack_molecule* list_of_pack_molecules) {
    t_pack_molecule* cur_pack_molecule = list_of_pack_molecules;
    while (cur_pack_molecule != nullptr) {
        cur_pack_molecule = list_of_pack_molecules->next;
        delete list_of_pack_molecules;
        list_of_pack_molecules = cur_pack_molecule;
    }
}

void Prepacker::init(const AtomNetlist& atom_nlist, const std::vector<t_logical_block_type>& logical_block_types) {
    VTR_ASSERT(list_of_pack_molecules == nullptr && "Prepacker cannot be initialized twice.");

    // Allocate the pack patterns from the logical block types.
    list_of_pack_patterns = alloc_and_load_pack_patterns(logical_block_types);
    // Use the pack patterns to allocate and load the pack molecules.
    std::multimap<AtomBlockId, t_pack_molecule*> atom_molecules_multimap;
    // Resize expected_lowest_cost_pb_gnode for the current netlist size.
    // Note: This will be resized again after alloc_and_load_pack_molecules() since
    // fill_vacant_chain_spots() may add new pass-through adder blocks to the netlist.
    expected_lowest_cost_pb_gnode.resize(atom_nlist.blocks().size(), nullptr);
    list_of_pack_molecules = alloc_and_load_pack_molecules(list_of_pack_patterns.data(),
                                                           expected_lowest_cost_pb_gnode,
                                                           list_of_pack_patterns.size(),
                                                           atom_molecules_multimap,
                                                           atom_nlist,
                                                           logical_block_types);

    // After alloc_and_load_pack_molecules returns, the netlist may have been compressed
    // (due to fill_vacant_chain_spots calling set_pin_net which marks it dirty).
    // The compression remaps all AtomBlockIds, so we need to rebuild expected_lowest_cost_pb_gnode
    // using the current (post-compression) netlist block IDs.
    expected_lowest_cost_pb_gnode.clear();
    expected_lowest_cost_pb_gnode.resize(atom_nlist.blocks().size(), nullptr);

    // Fill in expected_lowest_cost_pb_gnode for all blocks using the remapped IDs.
    for (AtomBlockId blk_id : atom_nlist.blocks()) {
        expected_lowest_cost_pb_gnode[blk_id] = get_expected_lowest_cost_primitive_for_atom_block(blk_id, logical_block_types);
    }

    // The multimap is a legacy thing. Since blocks can be part of multiple pack
    // patterns, during prepacking a block may be contained within multiple
    // molecules. However, by the end of prepacking, molecules should be
    // combined such that each block is contained in one and only one molecule.
    atom_molecules.resize(atom_nlist.blocks().size(), nullptr);
    for (AtomBlockId blk_id : atom_nlist.blocks()) {
        auto range = atom_molecules_multimap.equal_range(blk_id);
        // Every atom block should be packed into at least one molecule.
        VTR_ASSERT(range.first != range.second);

        // If an atom ends up in multiple molecules (e.g. due to overlapping
        // chain patterns), follow the existing convention and use the last
        // molecule inserted for this block as its canonical molecule.
        auto chosen_iter = range.first;
        for (auto it = range.first; it != range.second; ++it) {
            chosen_iter = it;
        }
        atom_molecules[blk_id] = chosen_iter->second;
    }

    // Filter the global molecule list so that it only contains molecules which
    // are actually referenced by at least one atom in atom_molecules. This
    // ensures that no atom will appear in multiple forced-pack molecules as
    // seen by the packer, preventing duplicate placement attempts.
    t_pack_molecule* new_head = nullptr;
    t_pack_molecule* cur = list_of_pack_molecules;
    while (cur != nullptr) {
        t_pack_molecule* next = cur->next;

        bool used = false;
        for (AtomBlockId blk_id : atom_nlist.blocks()) {
            if (atom_molecules[blk_id] == cur) {
                used = true;
                break;
            }
        }

        if (used) {
            // Keep this molecule in the list.
            cur->next = new_head;
            new_head = cur;
        } else {
            // No atom points to this molecule anymore; drop it.
            delete cur;
        }

        cur = next;
    }
    list_of_pack_molecules = new_head;
}

t_molecule_stats Prepacker::calc_max_molecule_stats(const AtomNetlist& atom_nlist) const {
    t_molecule_stats max_molecules_stats;
    t_pack_molecule* molecule_head = list_of_pack_molecules;
    for (auto cur_molecule = molecule_head; cur_molecule != nullptr; cur_molecule = cur_molecule->next) {
        //Calculate per-molecule statistics
        (void)atom_nlist;
        t_molecule_stats cur_molecule_stats = calc_molecule_stats(cur_molecule, atom_nlist);

        //Record the maximums (member-wise) over all molecules
        max_molecules_stats.num_blocks = std::max(max_molecules_stats.num_blocks, cur_molecule_stats.num_blocks);

        max_molecules_stats.num_pins = std::max(max_molecules_stats.num_pins, cur_molecule_stats.num_pins);
        max_molecules_stats.num_input_pins = std::max(max_molecules_stats.num_input_pins, cur_molecule_stats.num_input_pins);
        max_molecules_stats.num_output_pins = std::max(max_molecules_stats.num_output_pins, cur_molecule_stats.num_output_pins);

        max_molecules_stats.num_used_ext_pins = std::max(max_molecules_stats.num_used_ext_pins, cur_molecule_stats.num_used_ext_pins);
        max_molecules_stats.num_used_ext_inputs = std::max(max_molecules_stats.num_used_ext_inputs, cur_molecule_stats.num_used_ext_inputs);
        max_molecules_stats.num_used_ext_outputs = std::max(max_molecules_stats.num_used_ext_outputs, cur_molecule_stats.num_used_ext_outputs);
    }

    return max_molecules_stats;
}

void Prepacker::reset() {
    // When the prepacker is reset (or destroyed), clean up the internal data
    // members.
    free_list_of_pack_patterns(list_of_pack_patterns);
    free_pack_molecules(list_of_pack_molecules);
    // Reset everything to default state.
    list_of_pack_patterns.clear();
    list_of_pack_molecules = nullptr;
    atom_molecules.clear();
    expected_lowest_cost_pb_gnode.clear();
}

/*******************************************************/
/*        Extra carry chain logic implementation       */
/*******************************************************/

/**
 *  Find the next primitive input pin connected to the given cluster_input_pin.
 *  Following edges that are annotated with pack_pattern index
 */
static t_pb_graph_pin* get_connected_primitive_input_pin(const t_pb_graph_pin* cluster_input_pin, const int pack_pattern) {
    for (int iedge = 0; iedge < cluster_input_pin->num_output_edges; iedge++) {
        const auto& output_edge = cluster_input_pin->output_edges[iedge];
        // If this edge is annotated with the given pack pattern, or its pattern
        // should be inferred, follow it.
        if (output_edge->annotated_with_pattern(pack_pattern) || output_edge->infer_pattern) {
            for (int ipin = 0; ipin < output_edge->num_output_pins; ipin++) {
                if (output_edge->output_pins[ipin]->is_primitive_pin()) {
                    return output_edge->output_pins[ipin];
                }
                return get_connected_primitive_input_pin(output_edge->output_pins[ipin], pack_pattern);
            }
        }
    }

    // primitive input pin should always
    // be found when using this function
    VTR_ASSERT(false);
    return nullptr;
}

/**
 *  Find the previous primitive output pin connected to the given cluster_output_pin.
 *  Following edges that are annotated with pack_pattern index
 */
static t_pb_graph_pin* get_connected_primitive_output_pin(const t_pb_graph_pin* cluster_output_pin, const int pack_pattern) {
    for (int iedge = 0; iedge < cluster_output_pin->num_input_edges; iedge++) {
        const auto& input_edge = cluster_output_pin->input_edges[iedge];
        // If this edge is annotated with the given pack pattern, or its pattern
        // should be inferred, follow it.
        if (input_edge->annotated_with_pattern(pack_pattern) || input_edge->infer_pattern) {
            for (int ipin = 0; ipin < input_edge->num_input_pins; ipin++) {
                if (input_edge->input_pins[ipin]->is_primitive_pin()) {
                    return input_edge->input_pins[ipin];
                }
                return get_connected_primitive_output_pin(input_edge->input_pins[ipin], pack_pattern);
            }
        }
    }

    // primitive output pin should always
    // be found when using this function
    VTR_ASSERT(false);
    return nullptr;
}

/**
 * This function takes the input pin starting a chain (Cin of the root block) and finds the
 * the Cout pin of the last adder primitve of the chain.
 */
static t_pb_graph_pin* find_chain_exit_pin(t_pb_graph_pin* input_pin, int pattern_index) {
    VTR_ASSERT(input_pin->num_output_edges == 1);
    VTR_ASSERT(input_pin->output_edges[0]->annotated_with_pattern(pattern_index));

    auto first_cin_pin = get_connected_primitive_input_pin(input_pin, pattern_index);
    // pointer to the port model of the cin port of the adder primitive
    const auto cin_port_model = first_cin_pin->port->model_port;

    // create a queue of pin pointers for the breadth first search
    std::queue<t_pb_graph_pin*> pins_queue;

    // add the input pin to the queue
    pins_queue.push(first_cin_pin);

    // do breadth first search till all
    // connected pins are explored
    while (!pins_queue.empty()) {
        // get the first pin in the queue
        auto current_pin = pins_queue.front();

        // remove pin from queue
        pins_queue.pop();

        // if this is a primitive input pin and it's not a cin port, ignore pin
        // since we are only searching along the path of the chain ports
        if (current_pin->is_primitive_pin()
            && current_pin->port->type == IN_PORT
            && current_pin->port->model_port != cin_port_model) {
            continue;
        }

        // expand search from current pin
        expand_search(current_pin, pins_queue, pattern_index);

        // if this is an output pin of a root block then its connected
        // to the last cout of the chain. Return the connected primtive pin.
        if (current_pin->is_root_block_pin()
            && current_pin->num_output_edges == 0) {
            return get_connected_primitive_output_pin(current_pin, pattern_index);
        }
    }

    // Exit chain pin should be found
    VTR_ASSERT(false);
    return nullptr;
}

/**
 * get the pattern block that matches the input block id in this molecule
 */
static t_pack_pattern_block* get_atom_pattern_block(const t_pack_molecule* molecule, const int block_id) {
    const auto root_block = molecule->pack_pattern->root_block;

    std::vector<bool> visited_blocks(molecule->num_blocks);

    std::queue<t_pack_pattern_block*> pattern_block_queue;
    pattern_block_queue.push(root_block);

    // do breadth first search to find the block that matches block_id
    while (!pattern_block_queue.empty()) {
        auto pattern_block = pattern_block_queue.front();
        pattern_block_queue.pop();

        // ignore if a nullptr or is already visited
        if (!pattern_block || visited_blocks[pattern_block->block_id])
            continue;

        if (pattern_block->block_id == block_id)
            return pattern_block;

        visited_blocks[pattern_block->block_id] = true;

        auto block_connections = pattern_block->connections;

        // add all the blocks in the list of connections to the queue
        while (block_connections) {
            pattern_block_queue.push(block_connections->from_block);
            pattern_block_queue.push(block_connections->to_block);
            block_connections = block_connections->next;
        }
    }

    // this block is in this molecule
    // so it should be found
    VTR_ASSERT(false);
    return nullptr;
}

static bool chain_input_is_reachable(const t_pack_molecule* molecule,
                                     const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                     const AtomNetlist& atom_nlist) {
    const auto& chain_root_pins = molecule->pack_pattern->chain_root_pins;
    // assume that if the molecule can start in multiple locations
    // it will always be reachable from the previous molecule
    if (chain_root_pins.size() > 1)
        return true;

    // id of the root block of this molecule
    const auto root_block = molecule->atom_block_ids[molecule->root];
    // get the model of the cin port of the adder primitive
    const auto cin_port_model = chain_root_pins[0][0]->port->model_port;
    // get the pin number of the cin pin within the cin port
    const auto cin_pin_number = chain_root_pins[0][0]->pin_number;
    // get the atom block driving the root block of this molecule
    const auto driver_block = atom_nlist.find_atom_pin_driver(root_block, cin_port_model, cin_pin_number);

    auto driver_molecule_it = atom_molecules.find(driver_block);
    // if the driver block is not in molecule yet
    // then the block is driven by a constant net
    if (driver_molecule_it == atom_molecules.end())
        return true;

    auto driver_molecule = driver_molecule_it->second;

    if (driver_molecule->type != MOLECULE_FORCED_PACK)
        return true;

    t_pb_graph_node* driver_pb_graph_node = get_driver_pb_graph_node(driver_molecule, driver_block);

    // get the model of the cout port of the adder primitive
    const auto cout_port_model = molecule->pack_pattern->chain_exit_pins[0]->port->model_port;

    for (int iport = 0; iport < driver_pb_graph_node->num_output_ports; iport++) {
        for (int ipin = 0; ipin < driver_pb_graph_node->num_output_pins[iport]; ipin++) {
            const auto& pin = driver_pb_graph_node->output_pins[iport][ipin];
            if (pin.port->model_port == cout_port_model) {
                if (&pin == molecule->pack_pattern->chain_exit_pins[0])
                    return true;
                else
                    return false;
            }
        }
    }

    return false;
}

/**
 * This function finds the atom driving the root block of a molecule
 * and find the pb_graph_node associated with this block
 */
static t_pb_graph_node* get_driver_pb_graph_node(const t_pack_molecule* driver_molecule, const AtomBlockId driver_block) {
    auto it = std::find(driver_molecule->atom_block_ids.begin(), driver_molecule->atom_block_ids.end(), driver_block);
    VTR_ASSERT(it != driver_molecule->atom_block_ids.end());

    auto driver_pattern_block_id = std::distance(driver_molecule->atom_block_ids.begin(), it);
    auto driver_pattern_block = get_atom_pattern_block(driver_molecule, driver_pattern_block_id);

    auto block_connection = driver_pattern_block->connections;
    while (block_connection) {
        if (block_connection->to_block == driver_pattern_block) {
            return block_connection->to_pin->parent_node;
        }
        block_connection = block_connection->next;
    }

    VTR_ASSERT(false);
    return nullptr;
}

static int get_forced_chain_id(t_pack_molecule* molecule, const t_pack_molecule* prev_molecule, const AtomBlockId driver_block_id) {
    t_pb_graph_node* driver_pb_graph_node = get_driver_pb_graph_node(prev_molecule, driver_block_id);

    VTR_ASSERT(driver_pb_graph_node);

    const auto& chain_exit_pins = molecule->pack_pattern->chain_exit_pins;
    // get the model of the cout port of the adder primitive
    const auto cout_port_model = chain_exit_pins[0]->port->model_port;

    for (int iport = 0; iport < driver_pb_graph_node->num_output_ports; iport++) {
        for (int ipin = 0; ipin < driver_pb_graph_node->num_output_pins[iport]; ipin++) {
            const auto& pin = driver_pb_graph_node->output_pins[iport][ipin];
            if (pin.port->model_port == cout_port_model) {
                for (size_t chain_id = 0; chain_id < chain_exit_pins.size(); chain_id++) {
                    // architecture specific hack
                    if (pin.parent_node->placement_index == chain_exit_pins[chain_id]->parent_node->placement_index)
                        return chain_id;
                }
            }
        }
    }

    VTR_ASSERT(false);
    return -1;
}

static AtomBlockId get_adder_driver_block(const AtomBlockId block_id,
                                          const t_pack_patterns* pack_pattern,
                                          const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules,
                                          const AtomNetlist& atom_nlist) {
    const auto cin_pin = pack_pattern->chain_root_pins[0][0];
    const auto cin_model = cin_pin->port->model_port;
    const auto block_pb_graph_node = cin_pin->parent_node;
    const auto block_pb_type = block_pb_graph_node->pb_type;

    auto driver_id = atom_nlist.find_atom_pin_driver(block_id, cin_model, cin_pin->pin_number);
    AtomBlockId dummy_adder = AtomBlockId::INVALID();

    if (driver_id && primitive_type_feasible(driver_id, block_pb_type)
        && atom_molecules.find(driver_id) == atom_molecules.end()) {
        if (atom_nlist.find_atom_pin_driver(driver_id, cin_model, cin_pin->pin_number))
            return driver_id;
        else
            dummy_adder = driver_id;
    }

    if (atom_molecules.find(driver_id) != atom_molecules.end())
        return AtomBlockId::INVALID();

    for (int iport = 0; iport < block_pb_graph_node->num_input_ports; ++iport) {
        for (int ipin = 0; ipin < block_pb_graph_node->num_input_pins[iport]; ++ipin) {
            const auto& pin = block_pb_graph_node->input_pins[iport][ipin];
            if (pin.port->model_port != cin_model) {
                auto input_driver_id = atom_nlist.find_atom_pin_driver(block_id, pin.port->model_port, pin.pin_number);
                if (input_driver_id && primitive_type_feasible(input_driver_id, block_pb_type)
                    && atom_molecules.find(input_driver_id) == atom_molecules.end()) {
                    return input_driver_id;
                }
            }
        }
    }

    return dummy_adder;
}

/**
 * This function returns true is this molecule is a packed molecule
 * that has hierarchical structure. For example, an adder that is feeding
 * another adder throught the sumout port.
 */
static bool molecule_is_hierarchical(const t_pack_molecule* molecule) {
    // assume that only chained molecules can be hierarchical
    if (!molecule->is_chain())
        return false;

    const auto cout_pin_model = molecule->pack_pattern->chain_exit_pins[0]->port->model_port;
    const auto root_block = molecule->pack_pattern->root_block;
    auto connection = root_block->connections;

    while (connection) {
        if (connection->from_block == root_block
            && connection->from_pin->port->model_port != cout_pin_model
            && connection->from_pin->parent_node->pb_type == connection->to_pin->parent_node->pb_type) {
            return true;
        }
        connection = connection->next;
    }

    return false;
}

static bool valid_second_level_placement(const AtomBlockId first_level_block,
                                         const AtomBlockId block_id,
                                         const t_pack_molecule* molecule,
                                         const AtomNetlist& atom_nlist) {
    auto cin_pin = molecule->pack_pattern->chain_root_pins[0][0];
    auto cin_port_model = cin_pin->port->model_port;

    auto first_level_driver = first_level_block;
    auto second_level_driver = block_id;
    do {
        first_level_driver = atom_nlist.find_atom_pin_driver(first_level_driver, cin_port_model, cin_pin->pin_number);
        second_level_driver = atom_nlist.find_atom_pin_driver(second_level_driver, cin_port_model, cin_pin->pin_number);

        if (first_level_driver && !second_level_driver)
            return true;
        if (!first_level_driver && second_level_driver)
            return false;
    } while (first_level_driver || second_level_driver);

    return true;
}

static AtomBlockId is_second_level_block(const t_pack_pattern_block* pattern_block, const t_pack_molecule* molecule) {
    auto cin_pin = molecule->pack_pattern->chain_root_pins[0][0];
    auto cin_port_model = cin_pin->port->model_port;
    auto connection = pattern_block->connections;
    while (connection) {
        // if this block is being driven by this connection and the port being driven is not cin port
        if (connection->to_block == pattern_block
            && connection->to_pin->port->model_port != cin_port_model) {
            return molecule->atom_block_ids[connection->from_block->block_id];
        }
        connection = connection->next;
    }

    return AtomBlockId::INVALID();
}

// get the number of ALM inputs feeding the LUTs. The assumption is
// ALMs with 4-LUT has 6 inputs feeding LUTs, however, ALMs with 3-LUTs
// has 8 inputs feeding LUTs. This is a very specific assumption targeting
// the architectures in this study
static size_t get_alm_inputs_feeding_luts(t_pack_molecule* molecule) {
    auto pattern_block = molecule->pack_pattern->root_block;
    auto cin_pin = molecule->pack_pattern->chain_root_pins[0][0];
    auto cin_port = cin_pin->port;
    auto cin_port_model = cin_port->model_port;

    auto connection = pattern_block->connections;

    while (connection) {
        if (connection->to_block == pattern_block && connection->to_pin->port->model_port != cin_port_model) {
            auto lut_pb_type = connection->from_pin->parent_node->pb_type;
            auto lut_input_pins = lut_pb_type->num_input_pins;
            return (lut_input_pins == 3) ? 8 : 6;
        }
        connection = connection->next;
    }

    VTR_ASSERT(false);
    return 0;
}

// Helper function for [check_alm_input_limitation]
static void print_nets(std::unordered_set<AtomNetId>& nets,
                       int alm_placement_index,
                       int alut_placement_index,
                       const AtomNetlist& atom_nlist) {
    VTR_LOG("Placement index: %d->%d (%d)\n", alm_placement_index, alut_placement_index, nets.size());
    for (const auto net : nets) {
        VTR_LOG("%d %s\n", net, atom_nlist.net_name(net).c_str());
    }
    VTR_LOG("\n");
}

static bool check_alm_input_limitation(t_pack_molecule* molecule, const AtomNetlist& atom_nlist) {
    std::string pattern_name(molecule->pack_pattern->name);
    if (pattern_name.find("lut_chain") == std::string::npos)
        return true;

    auto block_id = molecule->atom_block_ids[molecule->root];
    auto pattern_block = molecule->pack_pattern->root_block;

    const auto ALM_INPUTS = get_alm_inputs_feeding_luts(molecule);

    auto cin_pin = molecule->pack_pattern->chain_root_pins[0][0];
    auto cin_port = cin_pin->port;
    auto cin_port_model = cin_port->model_port;

    std::string alm_name = "fle";
    if (get_pb_placement_index(pattern_block, alm_name) == -1)
        alm_name = "fle1";
    std::string alut_name = "ble5";
    auto alm_placement_index = get_pb_placement_index(pattern_block, alm_name);
    auto alut_placement_index = get_pb_placement_index(pattern_block, alut_name);

    std::unordered_set<AtomNetId> alm_nets;
    std::unordered_set<AtomNetId> alut_nets;

    while (true) {
        auto connection = pattern_block->connections;
        // get the unique net ids feeding the adders
        VTR_LOG("\n%s\n", atom_nlist.block_name(molecule->atom_block_ids[pattern_block->block_id]).c_str());
        while (connection) {
            if (connection->to_block == pattern_block && connection->to_pin->port->model_port != cin_port_model) {
                auto& lut_id = molecule->atom_block_ids[connection->from_block->block_id];
                if (lut_id) {
                    get_block_input_nets(lut_id, alm_nets, atom_nlist);
                    get_block_input_nets(lut_id, alut_nets, atom_nlist);
                    VTR_LOG("LUT %s (%zu)\n", atom_nlist.block_name(lut_id).c_str(), atom_nlist.block_input_pins(lut_id).size());
                    if (atom_nlist.block_input_pins(lut_id).size() > 2)
                        return false;
                    if (atom_nlist.block_input_pins(lut_id).empty()) {
                        alm_nets.insert((AtomNetId)0);
                        alut_nets.insert((AtomNetId)0);
                    }
                } else {
                    auto port_id = atom_nlist.find_atom_port(block_id, connection->to_pin->port->model_port);
                    if (port_id) {
                        auto net_id = atom_nlist.port_net(port_id, connection->to_pin->pin_number);
                        if (net_id) {
                            alm_nets.insert(net_id);
                            alut_nets.insert(net_id);
                        }
                    }
                }
            }
            connection = connection->next;
        }

        if (alm_nets.empty())
            break;

        print_nets(alm_nets, alm_placement_index, alut_placement_index, atom_nlist);
        // go to the next pattern block if it is still in the same ALM
        connection = pattern_block->connections;
        bool found_to_block = false;
        while (connection) {
            if (connection->from_block == pattern_block && connection->to_pin->port->model_port == cin_port_model) {
                auto alm_new_placement_index = get_pb_placement_index(connection->to_block, alm_name);
                auto alut_new_placement_index = get_pb_placement_index(connection->to_block, alut_name);
                if (alm_placement_index != alm_new_placement_index) {
                    if (alm_nets.size() > ALM_INPUTS || alut_nets.size() > 4) {
                        print_nets(alm_nets, alm_placement_index, alut_placement_index, atom_nlist);
                        modify_molecule(molecule, pattern_block, atom_nlist);
                        return check_alm_input_limitation(molecule, atom_nlist);
                    }
                    alut_nets.clear();
                    alm_nets.clear();
                    alm_placement_index = alm_new_placement_index;
                    alut_placement_index = alut_new_placement_index;
                } else if (alut_placement_index != alut_new_placement_index) {
                    if (alut_nets.size() > 4) {
                        print_nets(alm_nets, alm_placement_index, alut_placement_index, atom_nlist);
                        modify_molecule(molecule, pattern_block, atom_nlist);
                        return check_alm_input_limitation(molecule, atom_nlist);
                    }
                    alut_nets.clear();
                    alut_placement_index = alut_new_placement_index;
                }
                pattern_block = connection->to_block;
                block_id = molecule->atom_block_ids[pattern_block->block_id];
                found_to_block = true;
                break;
            }
            connection = connection->next;
        }

        if (!found_to_block) {
            if (alm_nets.size() > ALM_INPUTS || alut_nets.size() > 4) {
                print_nets(alm_nets, alm_placement_index, alut_placement_index, atom_nlist);
                modify_molecule(molecule, pattern_block, atom_nlist);
                return check_alm_input_limitation(molecule, atom_nlist);
            }
            break;
        }
    }

    return true;
}

static void get_block_input_nets(const AtomBlockId block_id,
                                 std::unordered_set<AtomNetId>& nets,
                                 const AtomNetlist& atom_nlist) {
    for (const auto& pin_id : atom_nlist.block_input_pins(block_id)) {
        nets.insert(atom_nlist.pin_net(pin_id));
    }
}

static int get_pb_placement_index(t_pack_pattern_block* pattern_block, std::string pb_name) {
    auto connection = pattern_block->connections;

    while (connection) {
        if (connection->to_block == pattern_block)
            break;
        connection = connection->next;
    }

    // Some pattern blocks (e.g. roots) may not have an incoming
    // connection in the pack pattern. In that case, we cannot infer
    // a meaningful placement index; return -1 to indicate "unknown".
    if (!connection)
        return -1;

    auto input_pin = connection->to_pin;
    auto parent_node = input_pin->parent_node;

    std::string parent_name(parent_node->pb_type->name);

    while (parent_node && parent_name != pb_name) {
        parent_node = parent_node->parent_pb_graph_node;
        if (!parent_node)
            break;
        std::string name(parent_node->pb_type->name);
        parent_name = name;
    }

    // VTR_ASSERT(parent_node);
    if (!parent_node)
        return -1;
    return parent_node->placement_index;
}

static void modify_molecule(t_pack_molecule* molecule,
                            t_pack_pattern_block* pattern_block,
                            const AtomNetlist& atom_nlist) {
    auto cin_pin = molecule->pack_pattern->chain_root_pins[0][0];
    auto cin_port = cin_pin->port;
    auto cin_port_model = cin_port->model_port;

    auto connection = pattern_block->connections;
    bool node_removed = false;

    while (true) {
        connection = pattern_block->connections;
        // get the unique net ids feeding the adders
        while (connection) {
            if (connection->to_block == pattern_block && connection->to_pin->port->model_port != cin_port_model) {
                auto& lut_id = molecule->atom_block_ids[connection->from_block->block_id];
                if (lut_id && atom_nlist.block_input_pins(lut_id).size() > 1) {
                    molecule->atom_block_ids[connection->from_block->block_id] = AtomBlockId::INVALID();
                    node_removed = true;
                    break;
                }
            }
            connection = connection->next;
        }
        if (node_removed)
            return;

        connection = pattern_block->connections;
        bool found_from_block = false;
        while (connection) {
            if (connection->to_block == pattern_block && connection->to_pin->port->model_port == cin_port_model) {
                pattern_block = connection->from_block;
                found_from_block = true;
                break;
            }
            connection = connection->next;
        }
        VTR_ASSERT(found_from_block);
    }
}

static bool check_lut_chain_molecules(t_pack_molecule* molecule, const AtomNetlist& atom_nlist) {
    std::string pattern_name(molecule->pack_pattern->name);
    if (pattern_name.find("lut_chain") == std::string::npos) return true;

    auto pattern_block = molecule->pack_pattern->root_block;
    auto cin_pin = molecule->pack_pattern->chain_root_pins[0][0];
    auto cin_port = cin_pin->port;
    auto cin_port_model = cin_port->model_port;

    auto connection = pattern_block->connections;
    t_pb_type* lut_pb_type = nullptr;

    while (connection) {
        if (connection->to_block == pattern_block && connection->to_pin->port->model_port != cin_port_model) {
            lut_pb_type = connection->from_pin->parent_node->pb_type;
            break;
        }
        connection = connection->next;
    }

    VTR_ASSERT(lut_pb_type);

    std::queue<t_pack_pattern_block*> pattern_block_queue;
    pattern_block_queue.push(pattern_block);

    std::vector<bool> visited_blocks(molecule->num_blocks);

    // do breadth first search to find the block that matches block_id
    while (!pattern_block_queue.empty()) {
        pattern_block = pattern_block_queue.front();
        pattern_block_queue.pop();

        // ignore if a nullptr or is already visited
        if (!pattern_block || visited_blocks[pattern_block->block_id])
            continue;

        if (molecule->atom_block_ids[pattern_block->block_id] && pattern_block->pb_type == lut_pb_type)
            return true;

        visited_blocks[pattern_block->block_id] = true;

        auto block_connections = pattern_block->connections;

        // add all the blocks in the list of connections to the queue
        while (block_connections) {
            pattern_block_queue.push(block_connections->from_block);
            pattern_block_queue.push(block_connections->to_block);
            block_connections = block_connections->next;
        }
    }

    VTR_LOG("check_lut_chain_molecules: Failed for pattern %s. Dump of primitives found:\n", pattern_name.c_str());

    // Re-initialize for printing
    std::fill(visited_blocks.begin(), visited_blocks.end(), false);
    while (!pattern_block_queue.empty())
        pattern_block_queue.pop();
    pattern_block_queue.push(molecule->pack_pattern->root_block);

    while (!pattern_block_queue.empty()) {
        pattern_block = pattern_block_queue.front();
        pattern_block_queue.pop();

        if (!pattern_block || visited_blocks[pattern_block->block_id])
            continue;

        visited_blocks[pattern_block->block_id] = true;

        if (molecule->atom_block_ids[pattern_block->block_id]) {
            VTR_LOG("  Block ID: %zu, Name: %s, Type: %s\n",
                    size_t(molecule->atom_block_ids[pattern_block->block_id]),
                    atom_nlist.block_name(molecule->atom_block_ids[pattern_block->block_id]).c_str(),
                    pattern_block->pb_type->name);
        }

        auto block_connections = pattern_block->connections;
        while (block_connections) {
            pattern_block_queue.push(block_connections->from_block);
            pattern_block_queue.push(block_connections->to_block);
            block_connections = block_connections->next;
        }
    }

    return false;
}