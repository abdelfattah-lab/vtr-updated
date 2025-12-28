/*
 * Intra-logic block router determines if a candidate packing solution (or intermediate solution) can route.
 *
 * Author: Jason Luu
 * Date: July 22, 2013
 */
#ifndef CLUSTER_ROUTER_H
#define CLUSTER_ROUTER_H
#include <string>
#include <vector>
#include "atom_netlist_fwd.h"
#include "pack_types.h"

/* Constructors/Destructors */
t_lb_router_data* alloc_and_load_router_data(std::vector<t_lb_type_rr_node>* lb_type_graph, t_logical_block_type_ptr type);
void free_router_data(t_lb_router_data* router_data);
void free_intra_lb_nets(std::vector<t_intra_lb_net>* intra_lb_nets);

/* Routing Functions */
void add_atom_as_target(t_lb_router_data* router_data, const AtomBlockId blk_id);
void remove_atom_from_target(t_lb_router_data* router_data, const AtomBlockId blk_id);
void set_reset_pb_modes(t_lb_router_data* router_data, const t_pb* pb, const bool set);
bool try_intra_lb_route(t_lb_router_data* router_data, int verbosity, t_mode_selection_status* mode_status);
void reset_intra_lb_route(t_lb_router_data* router_data);

/* Accessor Functions */
t_pb_routes alloc_and_load_pb_route(const std::vector<t_intra_lb_net>* intra_lb_nets, t_pb_graph_node* pb_graph_head);
void free_pb_route(t_pb_route* free_pb_route);

/**
 * @brief Get a description of why routing failed.
 *
 * @param router_data     The router data containing routing state
 * @param mode_status     The mode selection status from the routing attempt
 * @return A human-readable description of the routing failure
 */
std::string describe_routing_failure(const t_lb_router_data* router_data,
                                     const t_mode_selection_status& mode_status);

#endif
