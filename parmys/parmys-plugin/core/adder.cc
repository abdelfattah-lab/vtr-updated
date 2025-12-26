/*
 * Copyright 2022 CAS—Atlantic (University of New Brunswick, CASA)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#include "adder.h"
#include "multiplier.h"
#include "netlist_utils.h"
#include "node_utils.h"
#include "odin_globals.h"
#include "odin_types.h"
#include "odin_util.h"
#include "subtractor.h"
#include <string.h>
#include <vector>
#include <set>

#include "vtr_memory.h"
#include "vtr_util.h"

#include "parmys_utils.h"

using vtr::t_linked_vptr;

USING_YOSYS_NAMESPACE

t_model *hard_adders = NULL;
t_linked_vptr *add_list = NULL;
t_linked_vptr *processed_adder_list = NULL;
t_linked_vptr *chain_list = NULL;
int total = 0;
int *adder = NULL;
int min_add = 0;
int min_threshold_adder = 0;

void init_split_adder(nnode_t *node, nnode_t *ptr, int a, int sizea, int b, int sizeb, int cin, int cout, int index, int flag, netlist_t *netlist);
static void cleanup_add_old_node(nnode_t *nodeo, netlist_t *netlist);

/*---------------------------------------------------------------------------
 * (function: init_add_distribution)
 *  For adder, the output will only be the maxim input size + 1
 *-------------------------------------------------------------------------*/
void init_add_distribution()
{
    oassert(hard_adders != NULL);

    int len = hard_adders->inputs->size + hard_adders->inputs->next->size + 1;
    adder = (int *)vtr::calloc(len, sizeof(int));
}

/* These values are collected during the unused logic removal sweep */
extern long adder_chain_count;
extern long longest_adder_chain;
extern long total_adders;

extern double geomean_addsub_length;
extern double sum_of_addsub_logs;

void report_add_distribution()
{
    if (hard_adders == NULL)
        return;

    log("\nHard adder Distribution\n");
    log("============================\n");
    log("\n");
    log("\nTotal # of chains = %ld\n", adder_chain_count);

    log("\nHard adder chain Details\n");
    log("============================\n");

    log("\n");
    log("\nThe Number of Hard Block adders in the Longest Chain: %ld\n", longest_adder_chain);

    log("\n");
    log("\nThe Total Number of Hard Block adders: %ld\n", total_adders);

    log("\n");
    log("\nGeometric mean adder/subtractor chain length: %.2f\n", geomean_addsub_length);

    // vtr::free(adder);
}

/*---------------------------------------------------------------------------
 * (function: find_hard_adders)
 *-------------------------------------------------------------------------*/
void find_hard_adders()
{
    hard_adders = Arch.models;
    // Disable the size in configuration file.(The threshold for the extra bits).
    // min_add = configuration.min_hard_adder;
    min_threshold_adder = configuration.min_threshold_adder;

    while (hard_adders != NULL) {
        if (strcmp(hard_adders->name, "adder") == 0) {
            init_add_distribution();
            return;
        } else {
            hard_adders = hard_adders->next;
        }
    }

    return;
}

/*---------------------------------------------------------------------------
 * (function: declare_hard_adder)
 *-------------------------------------------------------------------------*/
void declare_hard_adder(nnode_t *node)
{
    t_adder *tmp;
    int width_a, width_b, width_sumout;

    /* See if this size instance of adder exists? */
    if (hard_adders == NULL)
        warning_message(NETLIST, node->loc, "%s\n", "Instantiating adder where adders do not exist");

    tmp = (t_adder *)hard_adders->instances;
    width_a = node->input_port_sizes[0];
    width_b = node->input_port_sizes[1];
    width_sumout = node->output_port_sizes[1];

    while (tmp != NULL) {
        if ((tmp->size_a == width_a) && (tmp->size_b == width_b) && (tmp->size_sumout == width_sumout))
            return;
        else
            tmp = tmp->next;
    }

    /* Does not exist - must create an instance */
    tmp = (t_adder *)vtr::malloc(sizeof(t_adder));
    tmp->next = (t_adder *)hard_adders->instances;
    hard_adders->instances = tmp;
    tmp->size_a = width_a;
    tmp->size_b = width_b;
    tmp->size_cin = 1;
    tmp->size_cout = 1;
    tmp->size_sumout = width_sumout;
    return;
}

/*---------------------------------------------------------------------------
 * (function: instantiate_hard_addier )
 *-------------------------------------------------------------------------*/
void instantiate_hard_adder(nnode_t *node, short mark, netlist_t * /*netlist*/)
{
    char *new_name;
    int len, sanity;

    declare_hard_adder(node);

    /* Need to give node proper name */
    len = strlen(node->name);
    len = len + 20; /* 20 chars should hold mul specs */
    new_name = (char *)vtr::malloc(len);

    /* wide input first :) identical branches! */
    // if (node->input_port_sizes[0] > node->input_port_sizes[1])
    // 	sanity = odin_sprintf(new_name, "%s", node->name);
    // else
    sanity = odin_sprintf(new_name, "%s", node->name);

    if (new_name)
        vtr::free(new_name);

    if (len <= sanity) /* buffer not large enough */
        oassert(false);

    /* Give names to the output pins */
    for (int i = 0; i < node->num_output_pins; i++) {
        if (node->output_pins[i]->name == NULL) {
            len = strlen(node->name) + 20; /* 6 chars for pin idx */
            new_name = (char *)vtr::malloc(len);
            odin_sprintf(new_name, "%s[%d]", node->name, node->output_pins[i]->pin_node_idx);
            node->output_pins[i]->name = new_name;
        }
    }

    node->traverse_visited = mark;
    return;
}

/*----------------------------------------------------------------------------
 * function: add_the_blackbox_for_adds()
 *--------------------------------------------------------------------------*/
void add_the_blackbox_for_adds_yosys(Yosys::Design *design)
{

    int hard_add_inputs, hard_add_outputs;
    t_adder *adds;
    t_model_ports *ports;
    char *pa, *pb, *psumout, *pcin, *pcout;

    /* Check to make sure this target architecture has hard adders */
    if (hard_adders == NULL)
        return;

    /* Get the names of the ports for the adder */
    ports = hard_adders->inputs;
    pcin = ports->name;
    ports = ports->next;
    pb = ports->name;
    ports = ports->next;
    pa = ports->name;

    ports = hard_adders->outputs;
    psumout = ports->name;
    ports = ports->next;
    pcout = ports->name;

    /* find the adder devices in the tech library */
    adds = (t_adder *)(hard_adders->instances);
    if (adds == NULL) /* No adders instantiated */
        return;

    /* simplified way of getting the multsize, but fine for quick example */
    while (adds != NULL) {

        Yosys::RTLIL::Module *module = nullptr;

        Yosys::hashlib::dict<Yosys::RTLIL::IdString, std::pair<int, bool>> wideports_cache;

        module = new Yosys::RTLIL::Module;
        module->name = Yosys::RTLIL::escape_id("adder");

        if (design->module(module->name))
            Yosys::log_error("Duplicate definition of module %s!\n", Yosys::log_id(module->name));
        design->add(module);

        /* add the inputs */
        hard_add_inputs = adds->size_a + adds->size_b + adds->size_cin;
        for (int i = 0; i < hard_add_inputs; i++) {
            std::string w_name;
            if (i < adds->size_a) {
                w_name = Yosys::stringf("%s[%d]", pa, i);
            } else if (i < hard_add_inputs - adds->size_cin && i >= adds->size_a) {
                w_name = Yosys::stringf("%s[%d]", pb, i - adds->size_a);
            } else {
                w_name = Yosys::stringf("%s[%d]", pcin, i - adds->size_a - adds->size_b);
            }

            Yosys::RTLIL::Wire *wire = to_wire(w_name, module);
            wire->port_input = true;

            std::pair<Yosys::RTLIL::IdString, int> wp = wideports_split(w_name);
            if (!wp.first.empty() && wp.second >= 0) {
                wideports_cache[wp.first].first = std::max(wideports_cache[wp.first].first, wp.second + 1);
                wideports_cache[wp.first].second = true;
            }
        }

        /* add the outputs */
        hard_add_outputs = adds->size_cout + adds->size_sumout;
        for (int i = 0; i < hard_add_outputs; i++) {
            std::string w_name;
            if (i < adds->size_cout) {
                w_name = Yosys::stringf("%s[%d]", pcout, i);
            } else {
                w_name = Yosys::stringf("%s[%d]", psumout, i - adds->size_cout);
            }

            Yosys::RTLIL::Wire *wire = to_wire(w_name, module);
            wire->port_output = true;

            std::pair<Yosys::RTLIL::IdString, int> wp = wideports_split(w_name);
            if (!wp.first.empty() && wp.second >= 0) {
                wideports_cache[wp.first].first = std::max(wideports_cache[wp.first].first, wp.second + 1);
                wideports_cache[wp.first].second = false;
            }
        }

        handle_wideports_cache(&wideports_cache, module);

        module->fixup_ports();
        wideports_cache.clear();

        module->attributes[Yosys::ID::blackbox] = Yosys::RTLIL::Const(1);

        adds = adds->next;
    }
}

void define_add_function_yosys(nnode_t *node, Yosys::Module *module, Yosys::Design *design)
{

    oassert(node->input_port_sizes[0] > 0);
    oassert(node->input_port_sizes[1] > 0);
    oassert(node->input_port_sizes[2] > 0);
    oassert(node->output_port_sizes[0] > 0);
    oassert(node->output_port_sizes[1] > 0);

    std::string cell_type_name = "adder";

    Yosys::IdString celltype = Yosys::RTLIL::escape_id(cell_type_name);
    Yosys::RTLIL::Cell *cell = module->addCell(NEW_ID, celltype);

    Yosys::hashlib::dict<Yosys::RTLIL::IdString, Yosys::hashlib::dict<int, Yosys::SigBit>> cell_wideports_cache;

    /* Write the input pins*/
    for (int i = 0; i < node->num_input_pins; i++) {
        std::string p, q;

        oassert(node->input_pins[i] != NULL && "input_pins[i] is NULL");
        oassert(node->input_pins[i]->net != NULL && "input_pins[i]->net is NULL");
        oassert(node->input_pins[i]->net->num_driver_pins == 1);
        npin_t *driver_pin = node->input_pins[i]->net->driver_pins[0];

        if (i < node->input_port_sizes[0]) {
            p = Yosys::stringf("%s[%d]", hard_adders->inputs->next->next->name, i);
            if (!driver_pin->name)
                q = driver_pin->node->name;
            else
                q = driver_pin->name;
        } else if (i >= node->input_port_sizes[0] && i < node->input_port_sizes[1] + node->input_port_sizes[0]) {
            p = Yosys::stringf("%s[%d]", hard_adders->inputs->next->name, i - node->input_port_sizes[0]);
            if (!driver_pin->name)
                q = driver_pin->node->name;
            else
                q = driver_pin->name;
        } else {
            p = Yosys::stringf("%s[%d]", hard_adders->inputs->name, i - (node->input_port_sizes[0] + node->input_port_sizes[1]));
            if (!driver_pin->name)
                q = driver_pin->node->name;
            else
                q = driver_pin->name;
        }

        std::pair<Yosys::RTLIL::IdString, int> wp = wideports_split(p);
        if (wp.first.empty())
            cell->setPort(Yosys::RTLIL::escape_id(p), to_wire(q, module));
        else
            cell_wideports_cache[wp.first][wp.second] = to_wire(q, module);
    }

    /* Write the output pins*/
    for (int i = 0; i < node->num_output_pins; i++) {
        std::string p, q;
        if (i < node->output_port_sizes[0]) {
            p = Yosys::stringf("%s[%d]", hard_adders->outputs->next->name, i);
            q = node->output_pins[i]->name;
        } else {
            p = Yosys::stringf("%s[%d]", hard_adders->outputs->name, i - node->output_port_sizes[0]);
            q = node->output_pins[i]->name;
        }

        std::pair<Yosys::RTLIL::IdString, int> wp = wideports_split(p);
        if (wp.first.empty())
            cell->setPort(Yosys::RTLIL::escape_id(p), to_wire(q, module));
        else
            cell_wideports_cache[wp.first][wp.second] = to_wire(q, module);
    }

    handle_cell_wideports_cache(&cell_wideports_cache, design, module, cell);

    return;
}

/*-----------------------------------------------------------------------
 * (function: init_split_adder)
 *###################################
 *	Create a carry chain adder when splitting. Inputs are connected
 *	to original pins, output pins are set to NULL for later connecting
 *	flag = 0: all adders are hard logic block; flag = 1: the last adder in the chain is soft logic block
 *---------------------------------------------------------------------*/
void init_split_adder(nnode_t *node, nnode_t *ptr, int a, int sizea, int b, int sizeb, int cin, int cout, int index, int flag, netlist_t *netlist)
{
    int flaga = 0, flagb = 0;
    int current_sizea, current_sizeb;
    int aa = 0, bb = 0, num = 0;

    // if the input of the first cin is generated by a dummy adder added
    // to the start of the chain, then an offset is needed to compensate
    // for that in various positions in the code, otherwise the offset is 0
    const int offset = (configuration.adder_cin_global) ? 0 : 1;

    /* Copy properties from original node */
    ptr->type = node->type;
    ptr->bit_width = node->bit_width;
    ptr->related_ast_node = node->related_ast_node;
    ptr->traverse_visited = node->traverse_visited;
    ptr->node_data = NULL;

    /* decide the current size of input a and b */
    if (flag == 0) {
        // increase input sizes by one if a dummy adder is
        // added to feed the first cin in the chain
        current_sizea = (a + offset) - sizea * index;
        current_sizeb = (b + offset) - sizeb * index;

        if (current_sizea >= sizea)
            current_sizea = sizea;
        else if (current_sizea <= 0) {
            current_sizea = sizea;
            flaga = 1;
        } else {
            aa = current_sizea;
            current_sizea = sizea;
            flaga = 2;
        }

        if (current_sizeb >= sizeb)
            current_sizeb = sizeb;
        else if (current_sizeb <= 0) {
            current_sizeb = sizeb;
            flagb = 1;
        } else {
            bb = current_sizeb;
            current_sizeb = sizeb;
            flagb = 2;
        }
    } else {
        if (sizea != 0)
            current_sizea = sizea;
        else
            current_sizea = 1;
        if (sizeb != 0)
            current_sizeb = sizeb;
        else
            current_sizeb = 1;
    }

    /* Set new port sizes and parameters */
    ptr->num_input_port_sizes = 3;
    ptr->input_port_sizes = (int *)vtr::malloc(3 * sizeof(int));
    ptr->input_port_sizes[0] = current_sizea;
    ptr->input_port_sizes[1] = current_sizeb;
    ptr->input_port_sizes[2] = cin;
    ptr->num_output_port_sizes = 2;
    ptr->output_port_sizes = (int *)vtr::malloc(2 * sizeof(int));
    ptr->output_port_sizes[0] = cout;

    /* The size of output port sumout equals the maxim size of sizea and sizeb  */
    if (current_sizea > current_sizeb)
        ptr->output_port_sizes[1] = current_sizea;
    else
        ptr->output_port_sizes[1] = current_sizeb;

    /* Set the number of pins and re-locate previous pin entries */
    ptr->num_input_pins = current_sizea + current_sizeb + cin;
    ptr->input_pins = (npin_t **)vtr::malloc(sizeof(void *) * (current_sizea + current_sizeb + cin));
    // if flaga or flagb = 1, the input pins should be empty.
    if (flaga == 1) {
        for (int i = 0; i < current_sizea; i++)
            ptr->input_pins[i] = NULL;
    } else if (flaga == 2) {
        if (index == 0) {
            ptr->input_pins[0] = NULL;
            if (sizea > 1) {
                for (int i = 1; i < aa; i++) {
                    ptr->input_pins[i] = node->input_pins[i + index * sizea - 1];
                    ptr->input_pins[i]->node = ptr;
                    ptr->input_pins[i]->pin_node_idx = i;
                }
                for (int i = 0; i < (sizea - aa); i++)
                    ptr->input_pins[i + aa] = NULL;
            }
        } else {
            for (int i = 0; i < aa; i++) {
                ptr->input_pins[i] = node->input_pins[i + index * sizea - 1];
                ptr->input_pins[i]->node = ptr;
                ptr->input_pins[i]->pin_node_idx = i;
            }
            for (int i = 0; i < (sizea - aa); i++)
                ptr->input_pins[i + aa] = NULL;
        }
    } else {
        if (index == 0 && !configuration.adder_cin_global) {
            if (flag == 0) {
                ptr->input_pins[0] = NULL;
                if (current_sizea > 1) {
                    for (int i = 1; i < current_sizea; i++) {
                        ptr->input_pins[i] = node->input_pins[i - 1];
                        ptr->input_pins[i]->node = ptr;
                        ptr->input_pins[i]->pin_node_idx = i;
                    }
                }
            } else {
                for (int i = 0; i < current_sizea; i++) {
                    ptr->input_pins[i] = node->input_pins[i];
                    ptr->input_pins[i]->node = ptr;
                    ptr->input_pins[i]->pin_node_idx = i;
                }
            }
        } else {
            if (flag == 0) {
                for (int i = 0; i < current_sizea; i++) {
                    // use the offset to compensate for the dummy adder added at start of the chain
                    ptr->input_pins[i] = node->input_pins[i + index * sizea - offset];
                    ptr->input_pins[i]->node = ptr;
                    ptr->input_pins[i]->pin_node_idx = i;
                }
            } else {
                if (sizea == 0)
                    connect_nodes(netlist->gnd_node, 0, ptr, 0);
                else {
                    num = node->input_port_sizes[0];
                    for (int i = 0; i < current_sizea; i++) {
                        ptr->input_pins[i] = node->input_pins[i + num - current_sizea];
                        ptr->input_pins[i]->node = ptr;
                        ptr->input_pins[i]->pin_node_idx = i;
                    }
                }
            }
        }
    }

    if (flagb == 1) {
        for (int i = 0; i < current_sizeb; i++)
            ptr->input_pins[i + current_sizeb] = NULL;
    } else if (flagb == 2) {
        if (index == 0) {
            ptr->input_pins[sizea] = NULL;
            if (current_sizeb > 1) {
                for (int i = 1; i < bb; i++) {
                    ptr->input_pins[i + current_sizea] = node->input_pins[i + a + index * sizeb - 1];
                    ptr->input_pins[i + current_sizea]->node = ptr;
                    ptr->input_pins[i + current_sizea]->pin_node_idx = i + current_sizea;
                }
                for (int i = 0; i < (sizeb - bb); i++)
                    ptr->input_pins[i + current_sizea + bb] = NULL;
            }
        } else {
            for (int i = 0; i < bb; i++) {
                ptr->input_pins[i + current_sizea] = node->input_pins[i + a + index * sizeb - 1];
                ptr->input_pins[i + current_sizea]->node = ptr;
                ptr->input_pins[i + current_sizea]->pin_node_idx = i + current_sizea;
            }
            for (int i = 0; i < (sizeb - bb); i++)
                ptr->input_pins[i + current_sizea + bb] = NULL;
        }
    } else {
        if (index == 0 && !configuration.adder_cin_global) {
            if (flag == 0) {
                ptr->input_pins[sizea] = NULL;
                if (current_sizeb > 1) {
                    for (int i = 1; i < current_sizeb; i++) {
                        ptr->input_pins[i + current_sizea] = node->input_pins[i + a + index * sizeb - 1];
                        ptr->input_pins[i + current_sizea]->node = ptr;
                        ptr->input_pins[i + current_sizea]->pin_node_idx = i + current_sizea;
                    }
                }
            } else {
                for (int i = 0; i < current_sizeb; i++) {
                    ptr->input_pins[i + current_sizea] = node->input_pins[i + a];
                    ptr->input_pins[i + current_sizea]->node = ptr;
                    ptr->input_pins[i + current_sizea]->pin_node_idx = i + current_sizea;
                }
            }
        } else {
            if (flag == 0) {
                for (int i = 0; i < current_sizeb; i++) {
                    ptr->input_pins[i + current_sizea] = node->input_pins[i + a + index * sizeb - offset];
                    ptr->input_pins[i + current_sizea]->node = ptr;
                    ptr->input_pins[i + current_sizea]->pin_node_idx = i + current_sizea;
                }
            } else {
                if (sizeb == 0)
                    connect_nodes(netlist->gnd_node, 0, ptr, current_sizea);
                else {
                    num = node->input_port_sizes[0] + node->input_port_sizes[1];
                    for (int i = 0; i < current_sizeb; i++) {
                        ptr->input_pins[i + current_sizea] = node->input_pins[i + num - current_sizeb];
                        ptr->input_pins[i + current_sizea]->node = ptr;
                        ptr->input_pins[i + current_sizea]->pin_node_idx = i + current_sizea;
                    }
                }
            }
        }
    }

    /* Carry_in should be NULL*/
    for (int i = 0; i < cin; i++) {
        ptr->input_pins[i + current_sizea + current_sizeb] = NULL;
    }

    /* output pins */
    int output;
    if (current_sizea > current_sizeb)
        output = current_sizea + cout;
    else
        output = current_sizeb + cout;

    ptr->num_output_pins = output;
    ptr->output_pins = (npin_t **)vtr::malloc(sizeof(void *) * output);
    for (int i = 0; i < output; i++)
        ptr->output_pins[i] = NULL;

    return;
}

/*-------------------------------------------------------------------------
 * (function: split_adder)
 *
 * This function works to split a adder into several smaller
 *  adders to better "fit" with the available resources in a
 *  targeted FPGA architecture.
 *
 * This function is at the lowest level since it simply receives
 *  a adder and is told how to split it.
 *
 * Note: In this function, we can do padding(default -1), fix the size of hard block adder.
 *-----------------------------------------------------------------------*/

void split_adder(nnode_t *nodeo, int a, int b, int sizea, int sizeb, int cin, int cout, int count, netlist_t *netlist)
{
    nnode_t **node;
    int num, lefta = 0, leftb = 0;
    int max_num = 0;
    int flag = 0;

    // if the input of the first cin is generated by a dummy adder added
    // to the start of the chain, then an offset is needed to compensate
    // for that in various positions in the code, otherwise the offset is 0
    const int offset = (configuration.adder_cin_global) ? 0 : 1;

    /* Check for a legitimate split */
    oassert(nodeo->input_port_sizes[0] == a);
    oassert(nodeo->input_port_sizes[1] == b);

    node = (nnode_t **)vtr::malloc(sizeof(nnode_t *) * (count));

    for (int i = 0; i < count; i++) {
        node[i] = allocate_nnode(nodeo->loc);
        node[i]->name = (char *)vtr::malloc(strlen(nodeo->name) + 20);
        odin_sprintf(node[i]->name, "%s-%d", nodeo->name, i);
        if (i == count - 1) {
            // fixed_hard_adder = 1 then adder need to be exact size;
            if (configuration.fixed_hard_adder == 1)
                init_split_adder(nodeo, node[i], a, sizea, b, sizeb, cin, cout, i, flag, netlist);
            else {
                if (count == 1) {
                    lefta = a;
                    leftb = b;
                } else {
                    lefta = (a + 1) % sizea;
                    leftb = (b + 1) % sizeb;
                }

                max_num = (lefta >= leftb) ? lefta : leftb;
                // if fixed_hard_adder = 0, and the left of a and b is more than min_add, then adder need to be remain the same size.
                if (max_num >= min_add)
                    init_split_adder(nodeo, node[i], a, sizea, b, sizeb, cin, cout, i, flag, netlist);
                else {
                    // Using soft logic to do the addition, No need to pad as the same size
                    flag = 1;
                    init_split_adder(nodeo, node[i], a, lefta, b, leftb, cin, cout, i, flag, netlist);
                }
            }
        } else
            init_split_adder(nodeo, node[i], a, sizea, b, sizeb, cin, cout, i, flag, netlist);

        // store the processed hard adder node for optimization
        processed_adder_list = insert_in_vptr_list(processed_adder_list, node[i]);
    }

    chain_information_t *adder_chain = allocate_chain_info();
    // if flag = 0, the last adder use soft logic, so the count of the chain should be one less
    if (flag == 0)
        adder_chain->count = count;
    else
        adder_chain->count = count - 1;
    adder_chain->num_bits = a + b;
    adder_chain->name = nodeo->name;
    chain_list = insert_in_vptr_list(chain_list, adder_chain);

    // don't add a dummy adder in the beginning of the chain if the first cin will be connected to a global gnd
    if ((flag == 0 || count > 1) && !configuration.adder_cin_global) {
        // connect the a[0] and b[0] of first adder node to ground
        connect_nodes(netlist->vcc_node, 0, node[0], 0);
        connect_nodes(netlist->gnd_node, 0, node[0], sizea);
        // hang the first sumout
        node[0]->output_pins[1] = allocate_npin();
        node[0]->output_pins[1]->name = append_string("", "%s~dummy_output~%d~%d", node[0]->name, 0, 1);
    }

    if (nodeo->num_input_port_sizes == 2) {
        // connect the first cin pin to unconn
        connect_nodes(netlist->pad_node, 0, node[0], node[0]->num_input_pins - 1);
    } else if (nodeo->num_input_port_sizes == 3) {
        // remap the first cin pins)
        remap_pin_to_new_node(nodeo->input_pins[nodeo->num_input_pins - 1], node[0], (node[0]->num_input_pins - 1));
    }
    // if (a + 1) % sizea == 0, the a[0] and b[0] of node[count-1] should connect to gound
    if ((a + 1) % sizea == 0 && (b + 1) % sizeb == 0) {
        if (flag == 0) {
            connect_nodes(netlist->gnd_node, 0, node[count - 1], 0);
            connect_nodes(netlist->gnd_node, 0, node[count - 1], sizea);
        }
    }

    // if any input pins beside first cin pins are NULL, connect those pins to unconn
    for (int i = 0; i < count; i++) {
        num = node[i]->num_input_pins;
        for (int j = 0; j < num - 1; j++) {
            if (node[i]->input_pins[j] == NULL)
                connect_nodes(netlist->pad_node, 0, node[i], j);
        }
    }

    if (configuration.adder_cin_global) {
        // connect first cin to gnd
        connect_nodes(netlist->gnd_node, 0, node[0], (node[0]->num_input_pins - 1));
    }

    // connect cout to next cin
    for (int i = 1; i < count; i++)
        connect_nodes(node[i - 1], 0, node[i], (node[i]->num_input_pins - 1));

    // remap the output pins of each adder to nodeo
    if (count == 1) {
        if (flag == 0) {
            for (int j = 0; j < node[0]->num_output_pins - 2; j++) {
                if (j < nodeo->num_output_pins)
                    remap_pin_to_new_node(nodeo->output_pins[j], node[0], j + 2);
                else {
                    node[0]->output_pins[j + 2] = allocate_npin();
                    node[0]->output_pins[j + 2]->name = append_string("", "%s~dummy_output~%d~%d", node[0]->name, 0, j + 2);
                }
                // hang the first cout
                node[0]->output_pins[0] = allocate_npin();
                node[0]->output_pins[0]->name = append_string("", "%s~dummy_output~%d~%d", node[0]->name, 0, 0);
            }
        } else {
            for (int j = 0; j < node[0]->num_output_pins - 1; j++)
                remap_pin_to_new_node(nodeo->output_pins[j], node[0], j + 1);
            remap_pin_to_new_node(nodeo->output_pins[nodeo->num_output_pins - 1], node[0], 0);
        }
    } else {
        // First adder
        for (int j = 0; j < node[0]->num_output_pins - 2; j++)
            remap_pin_to_new_node(nodeo->output_pins[j], node[0], j + 2);
        // if a dummy adder is added (offset = 1) start from the second adder)
        for (int i = offset; i < count - 1; i++) {
            for (int j = 0; j < node[i]->num_output_pins - 1; j++)
                remap_pin_to_new_node(nodeo->output_pins[i * sizea + j - offset], node[i], j + 1);
        }
        // Last adder
        if (flag == 0) {
            for (int j = 0; j < node[count - 1]->num_output_pins - 1; j++) {
                // if a dummy adder is added to this chain (offset = 1), adjust the index of the adder using the offset constant
                if (((count - 1) * sizea + j - offset) < nodeo->num_output_pins)
                    remap_pin_to_new_node(nodeo->output_pins[(count - 1) * sizea + j - offset], node[count - 1], j + 1);
                else {
                    node[count - 1]->output_pins[j + 1] = allocate_npin();
                    // Pad outputs with a unique and descriptive name to avoid collisions.
                    node[count - 1]->output_pins[j + 1]->name = append_string("", "%s~dummy_output~%d~%d", node[count - 1]->name, count - 1, j + 1);
                }
            }
            // Hang the last cout
            node[count - 1]->output_pins[0] = allocate_npin();
            // Pad outputs with a unique and descriptive name to avoid collisions.
            node[count - 1]->output_pins[0]->name = append_string("", "%s~dummy_output~%d~%d", node[count - 1]->name, count - 1, 0);
        } else {
            for (int j = 0; j < node[count - 1]->num_output_pins - 1; j++)
                // if(((count - 1) * sizea + j - 1) < nodeo->num_output_pins)
                remap_pin_to_new_node(nodeo->output_pins[(count - 1) * sizea + j - 1], node[count - 1], j + 1);
            if (nodeo->output_pins[nodeo->num_output_pins - 1] != NULL)
                remap_pin_to_new_node(nodeo->output_pins[nodeo->num_output_pins - 1], node[count - 1], 0);
            else {
                node[count - 1]->output_pins[0] = allocate_npin();
                // Pad outputs with a unique and descriptive name to avoid collisions.
                node[count - 1]->output_pins[0]->name = append_string("", "%s~dummy_output~%d~%d", node[count - 1]->name, count - 1, 0);
            }
        }
    }

    for (int i = offset; configuration.coarsen && i < count - 1; i++) {
        for (int j = 0; j < node[i]->num_output_pins - 1; j++) {
            char *new_output_pin_name = (char *)vtr::malloc((strlen(node[i]->name) + 20) * sizeof(char)); /* 6 chars for pin idx */
            odin_sprintf(new_output_pin_name, "%s[1]", node[i]->name);
            node[i]->output_pins[1]->name = new_output_pin_name;
        }
    }

    /* Freeing the old node! */
    cleanup_add_old_node(nodeo, netlist);

    vtr::free(node);
    return;
}

/*-------------------------------------------------------------------------
 * (function: is_pin_driven_by_adder_sumout)
 *
 * Check if a pin is driven by another adder's sumout.
 * This is used to detect ternary chain patterns where one adder's
 * result feeds another adder's input.
 *-----------------------------------------------------------------------*/
static bool is_pin_driven_by_adder_sumout(npin_t *pin)
{
    if (pin == NULL || pin->net == NULL)
        return false;

    nnet_t *net = pin->net;
    if (net->num_driver_pins == 0 || net->driver_pins == NULL)
        return false;

    // Check the first driver pin
    npin_t *driver_pin = net->driver_pins[0];
    if (driver_pin == NULL || driver_pin->node == NULL)
        return false;

    nnode_t *driver = driver_pin->node;

    // Check if driver is an ADD node and the pin is from sumout port
    if (driver->type == ADD) {
        // For ADD nodes, output port 0 is cout, port 1 is sumout
        // The sumout pins start after the cout pins
        int cout_size = driver->output_port_sizes[0];
        int pin_idx = driver_pin->pin_node_idx;
        if (pin_idx >= cout_size) {
            return true;  // This is a sumout pin
        }
    }

    return false;
}

/*-------------------------------------------------------------------------
 * (function: swap_adder_ports_for_chain_pattern)
 *
 * For the DCC3 architecture, the pack pattern expects sumout to feed
 * port B (adder[0].sumout -> adder[1].b). This function checks if
 * port A is fed by another adder's sumout and swaps ports A and B
 * so the pack pattern can recognize it.
 *-----------------------------------------------------------------------*/
static void swap_adder_ports_for_chain_pattern(nnode_t *node)
{
    if (node == NULL || node->type != ADD)
        return;

    int size_a = node->input_port_sizes[0];
    int size_b = node->input_port_sizes[1];

    // Count how many pins in port A are driven by adder sumout
    int sumout_count_a = 0;
    for (int i = 0; i < size_a; i++) {
        if (is_pin_driven_by_adder_sumout(node->input_pins[i])) {
            sumout_count_a++;
        }
    }

    // Count how many pins in port B are driven by adder sumout
    int sumout_count_b = 0;
    for (int i = 0; i < size_b; i++) {
        if (is_pin_driven_by_adder_sumout(node->input_pins[size_a + i])) {
            sumout_count_b++;
        }
    }

    // If port A has more sumout-driven pins than port B, swap the ports
    // This ensures sumout feeds port B for the pack pattern
    if (sumout_count_a > sumout_count_b) {
        // Swap port sizes
        node->input_port_sizes[0] = size_b;
        node->input_port_sizes[1] = size_a;

        // Swap pins - need temporary storage
        int total_pins = size_a + size_b;
        npin_t **temp_pins = (npin_t **)vtr::malloc(sizeof(npin_t *) * total_pins);

        // Copy port B pins to temp (will become port A)
        for (int i = 0; i < size_b; i++) {
            temp_pins[i] = node->input_pins[size_a + i];
            if (temp_pins[i] != NULL) {
                temp_pins[i]->pin_node_idx = i;
            }
        }
        // Copy port A pins to temp after port B (will become port B)
        for (int i = 0; i < size_a; i++) {
            temp_pins[size_b + i] = node->input_pins[i];
            if (temp_pins[size_b + i] != NULL) {
                temp_pins[size_b + i]->pin_node_idx = size_b + i;
            }
        }

        // Copy back to node
        for (int i = 0; i < total_pins; i++) {
            node->input_pins[i] = temp_pins[i];
        }

        vtr::free(temp_pins);
    }
}

/*-------------------------------------------------------------------------
 * Structures and functions for padding cascaded adder chains.
 *
 * When two $add operations are chained (one's sumout feeds another's b input),
 * they should form a ternary "chain" molecule. However, if the target chain
 * has DFF-only adders at the start that don't connect to the source chain's
 * sumout, VPR creates two independent "simple_chain" molecules instead.
 *
 * These functions detect such cascade-able pairs and add padding adders to
 * the source chain so the pattern matches what VPR expects.
 *-----------------------------------------------------------------------*/

// Structure to track a pair of cascaded adder chains
struct cascaded_adder_pair_t {
    nnode_t* source_chain_head;      // First adder in source chain (dummy or first real)
    nnode_t* target_chain_head;      // First adder in target chain (dummy or first real)
    int first_connection_pos;        // Target position where sumout connection starts (1-indexed)
    int num_dff_only;                // Number of DFF-only positions needing padding (= first_connection_pos - 1)
};

/*-------------------------------------------------------------------------
 * (function: get_adder_cin_driver)
 *
 * Returns the driver node of an adder's cin pin, or NULL if not found.
 *-----------------------------------------------------------------------*/
static nnode_t* get_adder_cin_driver(nnode_t* adder)
{
    if (adder == NULL || adder->type != ADD)
        return NULL;

    // cin is the last input pin
    int cin_idx = adder->num_input_pins - 1;
    npin_t* cin_pin = adder->input_pins[cin_idx];

    if (cin_pin == NULL || cin_pin->net == NULL)
        return NULL;

    nnet_t* net = cin_pin->net;
    if (net->num_driver_pins == 0 || net->driver_pins == NULL || net->driver_pins[0] == NULL)
        return NULL;

    return net->driver_pins[0]->node;
}

/*-------------------------------------------------------------------------
 * (function: get_adder_cout_fanout)
 *
 * Returns the first fanout node of an adder's cout pin that is an ADD,
 * or NULL if not found.
 *-----------------------------------------------------------------------*/
static nnode_t* get_adder_cout_fanout(nnode_t* adder)
{
    if (adder == NULL || adder->type != ADD)
        return NULL;

    // cout is output pin 0
    npin_t* cout_pin = adder->output_pins[0];
    if (cout_pin == NULL || cout_pin->net == NULL)
        return NULL;

    nnet_t* net = cout_pin->net;
    for (int i = 0; i < net->num_fanout_pins; i++) {
        if (net->fanout_pins[i] != NULL && net->fanout_pins[i]->node != NULL) {
            nnode_t* fanout = net->fanout_pins[i]->node;
            if (fanout->type == ADD) {
                return fanout;
            }
        }
    }
    return NULL;
}

/*-------------------------------------------------------------------------
 * (function: find_chain_head)
 *
 * Given an adder node, walks back through cin->cout connections to find
 * the chain head (the adder whose cin is not driven by another adder's cout).
 *-----------------------------------------------------------------------*/
static nnode_t* find_chain_head(nnode_t* adder)
{
    if (adder == NULL || adder->type != ADD)
        return NULL;

    nnode_t* current = adder;
    while (true) {
        nnode_t* driver = get_adder_cin_driver(current);
        // If driver is not an ADD, or is the same node, we've found the head
        if (driver == NULL || driver->type != ADD || driver == current) {
            return current;
        }
        current = driver;
    }
}

/*-------------------------------------------------------------------------
 * (function: get_chain_adder_at_position)
 *
 * Given a chain head, returns the adder at the specified position (0-indexed).
 * Position 0 is the head itself, position 1 is the next in the chain, etc.
 *-----------------------------------------------------------------------*/
static nnode_t* get_chain_adder_at_position(nnode_t* chain_head, int position)
{
    if (chain_head == NULL || position < 0)
        return NULL;

    nnode_t* current = chain_head;
    for (int i = 0; i < position; i++) {
        current = get_adder_cout_fanout(current);
        if (current == NULL)
            return NULL;
    }
    return current;
}

/*-------------------------------------------------------------------------
 * (function: get_chain_length)
 *
 * Returns the length of an adder chain starting from the given head.
 *-----------------------------------------------------------------------*/
static int get_chain_length(nnode_t* chain_head)
{
    if (chain_head == NULL || chain_head->type != ADD)
        return 0;

    int length = 1;
    nnode_t* current = chain_head;
    while (true) {
        nnode_t* next = get_adder_cout_fanout(current);
        if (next == NULL)
            break;
        length++;
        current = next;
    }
    return length;
}

/*-------------------------------------------------------------------------
 * (function: is_adder_b_port_dff_only)
 *
 * Returns true if all pins in the adder's b port are NOT driven by
 * any adder's sumout (i.e., they come from DFFs or other sources).
 *-----------------------------------------------------------------------*/
static bool is_adder_b_port_dff_only(nnode_t* adder)
{
    if (adder == NULL || adder->type != ADD)
        return false;

    int size_a = adder->input_port_sizes[0];
    int size_b = adder->input_port_sizes[1];

    // Check each pin in port B
    for (int i = 0; i < size_b; i++) {
        npin_t* pin = adder->input_pins[size_a + i];
        if (is_pin_driven_by_adder_sumout(pin)) {
            return false;  // At least one pin is driven by sumout
        }
    }
    return true;  // No pins driven by sumout
}

/*-------------------------------------------------------------------------
 * (function: get_b_port_sumout_driver_chain_head)
 *
 * If any pin in the adder's b port is driven by an adder's sumout,
 * returns the chain head of that source adder. Otherwise returns NULL.
 *-----------------------------------------------------------------------*/
static nnode_t* get_b_port_sumout_driver_chain_head(nnode_t* adder)
{
    if (adder == NULL || adder->type != ADD)
        return NULL;

    int size_a = adder->input_port_sizes[0];
    int size_b = adder->input_port_sizes[1];

    // Check each pin in port B for sumout connection
    for (int i = 0; i < size_b; i++) {
        npin_t* pin = adder->input_pins[size_a + i];
        if (pin != NULL && pin->net != NULL && pin->net->num_driver_pins > 0) {
            npin_t* driver_pin = pin->net->driver_pins[0];
            if (driver_pin != NULL && driver_pin->node != NULL) {
                nnode_t* driver = driver_pin->node;
                if (driver->type == ADD) {
                    // Check if driver pin is from sumout (not cout)
                    int cout_size = driver->output_port_sizes[0];
                    if (driver_pin->pin_node_idx >= cout_size) {
                        // This is a sumout pin, find the chain head
                        return find_chain_head(driver);
                    }
                }
            }
        }
    }
    return NULL;
}

/*-------------------------------------------------------------------------
 * (function: detect_cascaded_adder_pairs)
 *
 * Detects pairs of adder chains where one chain's sumout feeds another's
 * b port, but there are DFF-only adders at the start of the target chain.
 *-----------------------------------------------------------------------*/
static std::vector<cascaded_adder_pair_t> detect_cascaded_adder_pairs(netlist_t* netlist)
{
    std::vector<cascaded_adder_pair_t> pairs;
    std::set<nnode_t*> processed_targets;

    // Iterate through all processed adders
    t_linked_vptr* current = processed_adder_list;
    while (current != NULL) {
        nnode_t* adder = (nnode_t*)current->data_vptr;
        current = current->next;

        if (adder == NULL || adder->type != ADD)
            continue;

        // Find chain head for this adder
        nnode_t* target_head = find_chain_head(adder);
        if (target_head == NULL || processed_targets.count(target_head) > 0)
            continue;

        // Check if this target chain has any adder with b port fed by sumout
        nnode_t* source_head = NULL;
        int first_connection_pos = -1;

        int target_length = get_chain_length(target_head);
        for (int pos = 1; pos < target_length; pos++) {  // Skip position 0 (dummy)
            nnode_t* target_adder = get_chain_adder_at_position(target_head, pos);
            if (target_adder == NULL)
                break;

            nnode_t* driver_head = get_b_port_sumout_driver_chain_head(target_adder);
            if (driver_head != NULL && driver_head != target_head) {
                source_head = driver_head;
                first_connection_pos = pos;
                break;
            }
        }

        // If we found a connection and there are DFF-only positions before it
        if (source_head != NULL && first_connection_pos > 1) {
            // Verify positions 1 to (first_connection_pos - 1) are DFF-only
            bool all_dff_only = true;
            for (int pos = 1; pos < first_connection_pos; pos++) {
                nnode_t* target_adder = get_chain_adder_at_position(target_head, pos);
                if (target_adder == NULL || !is_adder_b_port_dff_only(target_adder)) {
                    all_dff_only = false;
                    break;
                }
            }

            if (all_dff_only) {
                cascaded_adder_pair_t pair;
                pair.source_chain_head = source_head;
                pair.target_chain_head = target_head;
                pair.first_connection_pos = first_connection_pos;
                pair.num_dff_only = first_connection_pos - 1;
                pairs.push_back(pair);
                processed_targets.insert(target_head);
            }
        }
    }

    return pairs;
}

/*-------------------------------------------------------------------------
 * (function: create_padding_adder)
 *
 * Creates a padding adder that passes through one input unchanged.
 * The padding adder computes: a + 0 + cin = a (passthrough).
 *-----------------------------------------------------------------------*/
static nnode_t* create_padding_adder(
    nnode_t* target_adder,
    int position,
    nnode_t* prev_adder,
    netlist_t* netlist)
{
    if (target_adder == NULL || netlist == NULL)
        return NULL;

    // Allocate new node
    nnode_t* padding = allocate_nnode(target_adder->loc);

    // Set name
    padding->name = (char*)vtr::malloc(strlen(target_adder->name) + 30);
    odin_sprintf(padding->name, "%s_padding_%d", target_adder->name, position);

    padding->type = ADD;

    // Get port sizes from target adder
    int size_a = target_adder->input_port_sizes[0];
    int size_b = target_adder->input_port_sizes[1];

    // For padding adder, we use size_b for port A (since we're taking target's b input)
    // and 1 for port B (will be connected to gnd)
    int padding_size = size_b;

    if (padding_size <= 0) {
        return NULL;
    }

    // Set bit_width to ensure this padding adder is treated as a hard adder
    // in partial_map_node() and not decomposed into soft logic (which would
    // leave input_pins NULL and cause crashes in define_add_function_yosys).
    // Use max of padding_size and min_threshold_adder to guarantee it passes
    // the threshold check in partial_map_node().
    padding->bit_width = (padding_size >= min_threshold_adder) ? padding_size : min_threshold_adder;

    // Allocate input ports: A (size_b pins) + B (size_b pins) + cin (1 pin)
    padding->num_input_port_sizes = 3;
    padding->input_port_sizes = (int*)vtr::malloc(sizeof(int) * 3);
    padding->input_port_sizes[0] = padding_size;  // port A
    padding->input_port_sizes[1] = padding_size;  // port B
    padding->input_port_sizes[2] = 1;             // cin

    padding->num_input_pins = padding_size + padding_size + 1;
    padding->input_pins = (npin_t**)vtr::calloc(padding->num_input_pins, sizeof(npin_t*));

    // Allocate output ports: cout (1 pin) + sumout (size_b pins)
    padding->num_output_port_sizes = 2;
    padding->output_port_sizes = (int*)vtr::malloc(sizeof(int) * 2);
    padding->output_port_sizes[0] = 1;            // cout
    padding->output_port_sizes[1] = padding_size; // sumout

    padding->num_output_pins = 1 + padding_size;
    padding->output_pins = (npin_t**)vtr::calloc(padding->num_output_pins, sizeof(npin_t*));

    // Connect port A to GND initially (will be reconnected in rewire_target_to_padding)
    // This ensures no input pins are left NULL even if rewiring fails
    for (int i = 0; i < padding_size; i++) {
        connect_nodes(netlist->gnd_node, 0, padding, i);
    }

    // Connect port B to GND (all pins)
    for (int i = 0; i < padding_size; i++) {
        connect_nodes(netlist->gnd_node, 0, padding, padding_size + i);
    }

    // Connect cin from previous adder's cout (or gnd if first)
    if (prev_adder != NULL && prev_adder->type == ADD) {
        connect_nodes(prev_adder, 0, padding, padding->num_input_pins - 1);
    } else {
        connect_nodes(netlist->gnd_node, 0, padding, padding->num_input_pins - 1);
    }

    // Allocate cout pin (output index 0) - MUST always exist for BLIF output
    // Even if nothing connects to this cout, define_add_function_yosys requires
    // all output pins to be valid with names.
    npin_t* cout_pin = allocate_npin();
    cout_pin->name = append_string("", "%s~cout~0", padding->name);
    cout_pin->type = OUTPUT;
    cout_pin->node = padding;
    cout_pin->pin_node_idx = 0;
    padding->output_pins[0] = cout_pin;

    // Create a net for cout (required even if unused)
    nnet_t* cout_net = allocate_nnet();
    cout_net->name = append_string("", "%s~cout", padding->name);
    add_driver_pin_to_net(cout_net, cout_pin);

    // Allocate sumout pins with nets (rewire step will add fanouts to these nets)
    for (int i = 0; i < padding_size; i++) {
        npin_t* sumout_pin = allocate_npin();
        sumout_pin->name = append_string("", "%s~sumout~%d", padding->name, i);
        sumout_pin->type = OUTPUT;
        sumout_pin->node = padding;
        sumout_pin->pin_node_idx = 1 + i;
        padding->output_pins[1 + i] = sumout_pin;

        // Create a net for this sumout - required for proper BLIF output
        nnet_t* sumout_net = allocate_nnet();
        sumout_net->name = vtr::strdup(sumout_pin->name);
        add_driver_pin_to_net(sumout_net, sumout_pin);
    }

    return padding;
}

/*-------------------------------------------------------------------------
 * (function: rewire_target_to_padding)
 *
 * Rewires the target adder's b port to receive from the padding adder's
 * sumout. The original driver of target's b port is moved to padding's a port.
 *
 * Note: create_padding_adder pre-connects port A to GND and creates nets
 * for sumout pins. This function rewires those connections.
 *-----------------------------------------------------------------------*/
static void rewire_target_to_padding(
    nnode_t* padding,
    nnode_t* target,
    netlist_t* netlist)
{
    if (padding == NULL || target == NULL || netlist == NULL)
        return;

    int size_a = target->input_port_sizes[0];
    int size_b = target->input_port_sizes[1];
    int padding_size = padding->input_port_sizes[0];

    // For each pin in target's b port
    for (int i = 0; i < size_b && i < padding_size; i++) {
        int target_b_idx = size_a + i;
        npin_t* target_pin = target->input_pins[target_b_idx];

        if (target_pin == NULL) {
            // If target has no b pin at this index, padding's a stays connected to GND
            continue;
        }

        nnet_t* original_net = target_pin->net;

        // Get padding's current port A pin (connected to GND from create_padding_adder)
        npin_t* old_padding_a_pin = padding->input_pins[i];

        // Move the original driver to padding's a port
        if (original_net != NULL && original_net->num_driver_pins > 0) {
            // Remove old padding A pin from GND net's fanout
            if (old_padding_a_pin != NULL && old_padding_a_pin->net != NULL) {
                nnet_t* gnd_net = old_padding_a_pin->net;
                for (int j = 0; j < gnd_net->num_fanout_pins; j++) {
                    if (gnd_net->fanout_pins[j] == old_padding_a_pin) {
                        for (int k = j; k < gnd_net->num_fanout_pins - 1; k++) {
                            gnd_net->fanout_pins[k] = gnd_net->fanout_pins[k + 1];
                            // Update pin_net_idx for shifted pins
                            if (gnd_net->fanout_pins[k] != NULL) {
                                gnd_net->fanout_pins[k]->pin_net_idx = k;
                            }
                        }
                        gnd_net->num_fanout_pins--;
                        break;
                    }
                }
            }

            // Create new pin for padding's a input
            npin_t* padding_a_pin = allocate_npin();
            padding_a_pin->type = INPUT;
            padding_a_pin->node = padding;
            padding_a_pin->pin_node_idx = i;
            padding->input_pins[i] = padding_a_pin;

            // Connect padding's a pin to the original net
            add_fanout_pin_to_net(original_net, padding_a_pin);

            // Remove target's b pin from the original net's fanout
            for (int j = 0; j < original_net->num_fanout_pins; j++) {
                if (original_net->fanout_pins[j] == target_pin) {
                    for (int k = j; k < original_net->num_fanout_pins - 1; k++) {
                        original_net->fanout_pins[k] = original_net->fanout_pins[k + 1];
                        // Update pin_net_idx for shifted pins
                        if (original_net->fanout_pins[k] != NULL) {
                            original_net->fanout_pins[k]->pin_net_idx = k;
                        }
                    }
                    original_net->num_fanout_pins--;
                    break;
                }
            }
        }

        // Get the sumout net (created in create_padding_adder)
        npin_t* sumout_pin = padding->output_pins[1 + i];
        nnet_t* sumout_net = sumout_pin->net;

        if (sumout_net == NULL) {
            // Shouldn't happen, but create one just in case
            sumout_net = allocate_nnet();
            sumout_net->name = append_string("", "%s~sumout~%d", padding->name, i);
            add_driver_pin_to_net(sumout_net, sumout_pin);
        }

        // Remove target's b pin from its old net (if different from sumout_net)
        if (target_pin->net != NULL && target_pin->net != sumout_net) {
            nnet_t* old_net = target_pin->net;
            for (int j = 0; j < old_net->num_fanout_pins; j++) {
                if (old_net->fanout_pins[j] == target_pin) {
                    for (int k = j; k < old_net->num_fanout_pins - 1; k++) {
                        old_net->fanout_pins[k] = old_net->fanout_pins[k + 1];
                        // Update pin_net_idx for shifted pins
                        if (old_net->fanout_pins[k] != NULL) {
                            old_net->fanout_pins[k]->pin_net_idx = k;
                        }
                    }
                    old_net->num_fanout_pins--;
                    break;
                }
            }
        }

        // Connect target's b pin to padding's sumout net
        target_pin->net = sumout_net;
        add_fanout_pin_to_net(sumout_net, target_pin);
    }
}

/*-------------------------------------------------------------------------
 * (function: transform_dummies_for_chain)
 *
 * Transforms the dummy adders for proper chain pattern:
 * - Source dummy: a=gnd (outputs 0 instead of 1)
 * - Target dummy: b=source_dummy.sumout (receives 0, same as original gnd)
 *-----------------------------------------------------------------------*/
static void transform_dummies_for_chain(
    nnode_t* source_dummy,
    nnode_t* target_dummy,
    netlist_t* netlist)
{
    if (source_dummy == NULL || target_dummy == NULL || netlist == NULL)
        return;

    // Source dummy: change a[0] from vcc to gnd so sumout = 0
    // First, disconnect the current a[0] input
    npin_t* source_a_pin = source_dummy->input_pins[0];
    if (source_a_pin != NULL && source_a_pin->net != NULL) {
        nnet_t* old_net = source_a_pin->net;
        // Remove from fanout list
        for (int i = 0; i < old_net->num_fanout_pins; i++) {
            if (old_net->fanout_pins[i] == source_a_pin) {
                for (int j = i; j < old_net->num_fanout_pins - 1; j++) {
                    old_net->fanout_pins[j] = old_net->fanout_pins[j + 1];
                    // Update pin_net_idx for shifted pins
                    if (old_net->fanout_pins[j] != NULL) {
                        old_net->fanout_pins[j]->pin_net_idx = j;
                    }
                }
                old_net->num_fanout_pins--;
                break;
            }
        }
        source_a_pin->net = NULL;
    }
    // Connect to gnd
    connect_nodes(netlist->gnd_node, 0, source_dummy, 0);

    // Target dummy: change b[0] to receive from source_dummy's sumout
    int target_size_a = target_dummy->input_port_sizes[0];
    npin_t* target_b_pin = target_dummy->input_pins[target_size_a];  // b[0]

    if (target_b_pin != NULL && target_b_pin->net != NULL) {
        nnet_t* old_net = target_b_pin->net;
        // Remove from fanout list
        for (int i = 0; i < old_net->num_fanout_pins; i++) {
            if (old_net->fanout_pins[i] == target_b_pin) {
                for (int j = i; j < old_net->num_fanout_pins - 1; j++) {
                    old_net->fanout_pins[j] = old_net->fanout_pins[j + 1];
                    // Update pin_net_idx for shifted pins
                    if (old_net->fanout_pins[j] != NULL) {
                        old_net->fanout_pins[j]->pin_net_idx = j;
                    }
                }
                old_net->num_fanout_pins--;
                break;
            }
        }
        target_b_pin->net = NULL;
    }

    // Create net from source_dummy sumout[0] to target_dummy b[0]
    npin_t* source_sumout_pin = source_dummy->output_pins[1];  // sumout is at index 1
    if (source_sumout_pin == NULL) {
        // Allocate if not present
        source_sumout_pin = allocate_npin();
        source_sumout_pin->name = append_string("", "%s~sumout~0", source_dummy->name);
        source_sumout_pin->type = OUTPUT;
        source_sumout_pin->node = source_dummy;
        source_sumout_pin->pin_node_idx = 1;
        source_dummy->output_pins[1] = source_sumout_pin;
    }

    nnet_t* chain_net = allocate_nnet();
    chain_net->name = append_string("", "%s~dummy_chain~%s", source_dummy->name, target_dummy->name);

    // Connect source sumout to net as driver
    if (source_sumout_pin->net == NULL) {
        add_driver_pin_to_net(chain_net, source_sumout_pin);
    } else {
        // Use existing net
        chain_net = source_sumout_pin->net;
    }

    // Connect target b[0] to net as fanout
    if (target_b_pin == NULL) {
        target_b_pin = allocate_npin();
        target_b_pin->type = INPUT;
        target_b_pin->node = target_dummy;
        target_b_pin->pin_node_idx = target_size_a;
        target_dummy->input_pins[target_size_a] = target_b_pin;
    }
    target_b_pin->net = chain_net;
    add_fanout_pin_to_net(chain_net, target_b_pin);
}

/*-------------------------------------------------------------------------
 * (function: rewire_source_cin_to_padding)
 *
 * Rewires the source chain adder's cin to receive from the padding adder's
 * cout instead of from its original driver (previous source adder or dummy).
 * This ensures the source chain flows through the padding adder.
 *-----------------------------------------------------------------------*/
static void rewire_source_cin_to_padding(
    nnode_t* source_adder,
    nnode_t* padding,
    netlist_t* netlist)
{
    if (source_adder == NULL || padding == NULL || netlist == NULL)
        return;

    // cin is the last input pin
    int cin_idx = source_adder->num_input_pins - 1;
    npin_t* cin_pin = source_adder->input_pins[cin_idx];

    if (cin_pin == NULL)
        return;

    // Disconnect cin from its current net
    nnet_t* old_net = cin_pin->net;
    if (old_net != NULL) {
        // Remove from fanout list
        for (int i = 0; i < old_net->num_fanout_pins; i++) {
            if (old_net->fanout_pins[i] == cin_pin) {
                for (int j = i; j < old_net->num_fanout_pins - 1; j++) {
                    old_net->fanout_pins[j] = old_net->fanout_pins[j + 1];
                    // Update pin_net_idx for shifted pins
                    if (old_net->fanout_pins[j] != NULL) {
                        old_net->fanout_pins[j]->pin_net_idx = j;
                    }
                }
                old_net->num_fanout_pins--;
                break;
            }
        }
        cin_pin->net = NULL;
    }

    // Get padding's cout net (cout is output pin 0)
    npin_t* padding_cout_pin = padding->output_pins[0];
    if (padding_cout_pin == NULL || padding_cout_pin->net == NULL)
        return;

    nnet_t* cout_net = padding_cout_pin->net;

    // Connect source's cin to padding's cout net
    cin_pin->net = cout_net;
    add_fanout_pin_to_net(cout_net, cin_pin);
}

/*-------------------------------------------------------------------------
 * (function: pad_cascaded_adder_chains)
 *
 * Main entry point for padding cascaded adder chains.
 * Detects cascade-able pairs and adds padding adders to create proper
 * ternary chain patterns that VPR can recognize.
 *-----------------------------------------------------------------------*/
void pad_cascaded_adder_chains(netlist_t* netlist)
{
    if (netlist == NULL || hard_adders == NULL)
        return;

    // Detect cascaded pairs that need padding
    std::vector<cascaded_adder_pair_t> pairs = detect_cascaded_adder_pairs(netlist);

    if (pairs.empty())
        return;

    log("Padding %zu cascaded adder chain pair(s) for ternary pattern\n", pairs.size());

    for (auto& pair : pairs) {
        if (pair.num_dff_only == 0)
            continue;

        log("  Padding pair: source=%s, target=%s, dff_only_positions=%d\n",
            pair.source_chain_head->name,
            pair.target_chain_head->name,
            pair.num_dff_only);

        // Transform dummy adders
        transform_dummies_for_chain(
            pair.source_chain_head,
            pair.target_chain_head,
            netlist);

        // Create padding adders for each DFF-only position
        nnode_t* prev_padding = pair.source_chain_head;  // Start from source dummy

        for (int pos = 1; pos <= pair.num_dff_only; pos++) {
            nnode_t* target_adder = get_chain_adder_at_position(pair.target_chain_head, pos);
            if (target_adder == NULL) {
                log("    Warning: Could not find target adder at position %d\n", pos);
                continue;
            }

            // Create padding adder
            nnode_t* padding = create_padding_adder(target_adder, pos, prev_padding, netlist);
            if (padding == NULL) {
                log("    Warning: Failed to create padding adder at position %d\n", pos);
                continue;
            }

            // Rewire target's b port to receive from padding's sumout
            rewire_target_to_padding(padding, target_adder, netlist);

            // Rewire source chain's adder at this position to take cin from padding's cout
            nnode_t* source_adder = get_chain_adder_at_position(pair.source_chain_head, pos);
            if (source_adder != NULL) {
                rewire_source_cin_to_padding(source_adder, padding, netlist);
            } else {
                log("    Warning: Could not find source adder at position %d for cin rewiring\n", pos);
            }

            // Add padding to processed list
            processed_adder_list = insert_in_vptr_list(processed_adder_list, padding);

            // Add padding to netlist's internal_nodes so it's properly tracked
            netlist->internal_nodes = (nnode_t**)vtr::realloc(netlist->internal_nodes,
                sizeof(nnode_t*) * (netlist->num_internal_nodes + 1));
            netlist->internal_nodes[netlist->num_internal_nodes++] = padding;

            prev_padding = padding;
            log("    Created padding adder: %s\n", padding->name);
        }
    }
}

/*-------------------------------------------------------------------------
 * (function: iterate_adders)
 *
 * This function will iterate over all of the add operations that
 *	exist in the netlist and perform a splitting so that they can
 *	fit into a basic hard adder block that exists on the FPGA.
 *	If the proper option is set, then it will be expanded as well
 *	to just use a fixed size hard adder.
 *-----------------------------------------------------------------------*/
void iterate_adders(netlist_t *netlist)
{
    int sizea, sizeb, sizecin; // the size of
    int a, b;
    int count, counta, countb;
    int num = 0;
    nnode_t *node;

    // offset to the adder size in case a dummy adder is added to
    // start of the adder chain to feed the first cin with gnd
    const int offset = (configuration.adder_cin_global) ? 0 : 1;

    /* Can only perform the optimization if hard adders exist! */
    if (hard_adders == NULL)
        return;
    // In hard block adder, the summand and addend are same size.
    sizecin = hard_adders->inputs->size;
    sizeb = hard_adders->inputs->next->size;
    sizea = hard_adders->inputs->next->size;

    oassert(sizecin == 1);

    while (add_list != NULL) {
        node = (nnode_t *)add_list->data_vptr;
        add_list = delete_in_vptr_list(add_list);
        oassert(node != NULL);
        if (node->type == HARD_IP)
            node->type = ADD;

        oassert(node->type == ADD);

        // Swap ports if needed to ensure sumout feeds port B for chain pattern
        swap_adder_ports_for_chain_pattern(node);

        a = node->input_port_sizes[0];
        b = node->input_port_sizes[1];
        num = (a >= b) ? a : b;
        node->bit_width = num;
        if (num >= min_threshold_adder && num >= min_add) {
            // if the first cin in a chain is fed by a global input (offset = 0) the adder width is the
            // input width + 1 (to pass the last cout -> sumout) divided by size of the adder input ports
            // otherwise (offset = 1) a dummy adder is added to the chain to feed the first cin with gnd
            // how many adders a can split
            counta = (a + 1) / sizea + offset;
            // how many adders b can split
            countb = (b + 1) / sizeb + offset;
            // how many adders need to be split
            if (counta >= countb)
                count = counta;
            else
                count = countb;
            total++;
            split_adder(node, a, b, sizea, sizeb, 1, 1, count, netlist);
        }
        // Store the node into processed_adder_list if the threshold is bigger than num
        else
            processed_adder_list = insert_in_vptr_list(processed_adder_list, node);
    }

    // After all adders are split, pad cascaded chains for ternary pattern recognition
    pad_cascaded_adder_chains(netlist);

    return;
}

/*-------------------------------------------------------------------------
 * (function: clean_adders)
 *
 * Clean up the memory by deleting the list structure of adders
 *	during optimization
 *-----------------------------------------------------------------------*/
void clean_adders()
{
    while (add_list != NULL)
        add_list = delete_in_vptr_list(add_list);
    return;
}

/*-------------------------------------------------------------------------
 * (function: reduce_operations)
 *
 * reduce the operations that are redundant
 *-----------------------------------------------------------------------*/
void reduce_operations(netlist_t * /*netlist*/, operation_list op)
{
    t_linked_vptr *place = NULL;
    operation_list oper;
    switch (op) {
    case ADD:
        place = add_list;
        oper = ADD;
        break;

    case MULTIPLY:
        place = mult_list;
        oper = MULTIPLY;
        break;

    case MINUS:
        place = sub_list;
        oper = MINUS;
        break;

    default:
        oper = NO_OP;
        break;
    }

    traverse_list(oper, place);
}

/*-------------------------------------------------------------------------
 * (function: traverse_list)
 *
 * traverse the operation lists
 *-----------------------------------------------------------------------*/
void traverse_list(operation_list oper, t_linked_vptr *place)
{
    while (place != NULL && place->next != NULL) {
        match_node(place, oper);
        place = place->next;
    }
}

/*---------------------------------------------------------------------------
 * (function: match_node)
 *-------------------------------------------------------------------------*/
void match_node(t_linked_vptr *place, operation_list oper)
{
    int flag = 0;
    int mark = 0;
    nnode_t *node = NULL;
    nnode_t *next_node = NULL;
    node = (nnode_t *)place->data_vptr;
    t_linked_vptr *pre = place;
    t_linked_vptr *next = NULL;
    if (place->next != NULL)
        next = place->next;
    while (next != NULL) {
        flag = 0;
        mark = 0;
        next_node = (nnode_t *)next->data_vptr;
        if (node->type == next_node->type) {
            if (node->num_input_pins == next_node->num_input_pins) {
                flag = match_ports(node, next_node, oper);
                if (flag == 1) {
                    mark = match_pins(node, next_node);
                    if (mark == 1) {
                        merge_nodes(node, next_node);
                        remove_list_node(pre, next);
                    }
                }
            }
        }
        if (mark == 1)
            next = pre->next;
        else {
            pre = next;
            next = next->next;
        }
    }
}

/*---------------------------------------------------------------------------
 * (function: match_ports)
 *-------------------------------------------------------------------------*/
int match_ports(nnode_t *node, nnode_t *next_node, operation_list oper)
{
    int flag = 0;
    int sign = 0;
    int mark1 = 1;
    int mark2 = 1;
    ast_node_t *ast_node, *ast_node_next;
    char *component_s[2] = {0};
    char *component_o[2] = {0};
    ast_node = node->related_ast_node;
    ast_node_next = next_node->related_ast_node;
    /* in case of coarsen blifs, there is no related ast node, so we skip this part */
    if (ast_node && ast_node->types.operation.op == oper) {
        traverse_operation_node(ast_node, component_s, oper, &sign);
        if (sign != 1) {
            traverse_operation_node(ast_node_next, component_o, oper, &sign);
            if (sign != 1) {
                oassert(component_s[0] && component_o[0] && "missing children on operation");
                switch (oper) {
                case ADD:
                case MULTIPLY: {
                    mark1 = strcmp(component_s[0], component_o[0]);
                    if (component_s[1] && component_o[1]) {
                        if (mark1 == 0) {
                            mark2 = strcmp(component_s[1], component_o[1]);
                        } else {
                            mark1 = strcmp(component_s[0], component_o[1]);
                            mark2 = strcmp(component_s[1], component_o[0]);
                        }
                    }
                } break;

                case MINUS: {
                    mark1 = strcmp(component_s[0], component_o[0]);
                    if (mark1 == 0 && component_s[1] && component_o[1]) {
                        mark2 = strcmp(component_s[1], component_o[1]);
                    }
                } break;

                default:

                    break;
                }
                if (mark1 == 0 && mark2 == 0) {
                    flag = 1;
                }
            }
        }
        for (int i = 0; i < ast_node->num_children; i++) {
            if (ast_node->children[i]->type != IDENTIFIERS) {
                vtr::free(component_s[i]);
            }
        }
        for (int i = 0; i < ast_node_next->num_children; i++) {
            if (ast_node_next->children[i]->type != IDENTIFIERS) {
                vtr::free(component_o[i]);
            }
        }
    }

    return flag;
}

/*-------------------------------------------------------------------------
 * (function: traverse_operation_node)
 *
 * search the ast find the couple of components
 *-----------------------------------------------------------------------*/
void traverse_operation_node(ast_node_t *node, char *component[], operation_list op, int *mark)
{
    if (node == NULL)
        return;

    if (node->types.operation.op == op) {
        for (long i = 0; i < node->num_children; i++) {
            *mark = 0;
            if (node->children[i]->type != IDENTIFIERS && node->children[i]->type != NUMBERS) {
                *mark = 1;
                break;
            } else {
                if (node->children[i]->type == IDENTIFIERS) {
                    component[i] = node->children[i]->types.identifier;
                } else if (node->children[i]->type == NUMBERS) {
                    long value = node->children[i]->types.vnumber->get_value();
                    long len = snprintf(NULL, 0, "%ld", value);
                    component[i] = (char *)vtr::calloc(len + 1, sizeof(char));
                    odin_sprintf(component[i], "%ld", value);
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------
 * (function: merge_node)
 *-------------------------------------------------------------------------*/
void merge_nodes(nnode_t *node, nnode_t *next_node)
{
    remove_fanout_pins(next_node);
    reallocate_pins(node, next_node);
    free_op_nodes(next_node);
}

/*---------------------------------------------------------------------------
 * (function: remove_list_node)
 *-------------------------------------------------------------------------*/
void remove_list_node(t_linked_vptr *pre, t_linked_vptr *next)
{
    if (next->next != NULL)
        pre->next = next->next;
    else
        pre->next = NULL;
    vtr::free(next);
}

/*---------------------------------------------------------------------------
 * (function: remove_fanout_pins)
 *-------------------------------------------------------------------------*/
void remove_fanout_pins(nnode_t *node)
{
    for (int i = 0; i < node->num_input_pins; i++) {
        int j, k;
        int idx = node->input_pins[i]->unique_id;
        for (j = 0; j < node->input_pins[i]->net->num_fanout_pins; j++) {
            if (node->input_pins[i]->net->fanout_pins[j]->unique_id == idx)
                break;
        }
        for (k = j; k < node->input_pins[i]->net->num_fanout_pins - 1; k++) {
            node->input_pins[i]->net->fanout_pins[k] = node->input_pins[i]->net->fanout_pins[k + 1];
            node->input_pins[i]->net->fanout_pins[k]->pin_net_idx = k;
        }
        node->input_pins[i]->net->fanout_pins[k] = NULL;
        node->input_pins[i]->net->num_fanout_pins--;
    }
}

/*---------------------------------------------------------------------------
 * (function: reallocate_pins)
 *-------------------------------------------------------------------------*/
void reallocate_pins(nnode_t *node, nnode_t *next_node)
{
    int pin_idx;
    nnode_t *input_node = NULL;
    nnet_t *net = NULL;
    npin_t *pin = NULL;
    for (int i = 0; i < next_node->num_output_pins; i++) {
        for (int j = 0; j < next_node->output_pins[i]->net->num_fanout_pins; j++) {
            if (next_node->output_pins[i]->net->fanout_pins[j]->node != NULL) {
                input_node = next_node->output_pins[i]->net->fanout_pins[j]->node;
                net = node->output_pins[i]->net;
                pin_idx = next_node->output_pins[i]->net->fanout_pins[j]->pin_node_idx;
                pin = input_node->input_pins[pin_idx];
                add_fanout_pin_to_net(net, pin);
            } else {
                free_npin(next_node->output_pins[i]->net->fanout_pins[j]);
            }
        }
    }
}

/*---------------------------------------------------------------------------
 * (function: free_op_nodes)
 *-------------------------------------------------------------------------*/
void free_op_nodes(nnode_t *node)
{
    for (int i = 0; i < node->num_output_pins; i++) {
        if (node->output_pins[i]->net != NULL) {
            free_nnet(node->output_pins[i]->net);
        }
    }
    free_nnode(node);
}

/*---------------------------------------------------------------------------
 * (function: match_pins)
 *-------------------------------------------------------------------------*/
int match_pins(nnode_t *node, nnode_t *next_node)
{
    for (int i = 0; i < node->num_input_pins; i++) {
        for (int j = 0; j < node->input_pins[i]->net->num_driver_pins; j++) {
            bool found = false;
            long id = node->input_pins[i]->net->driver_pins[j]->unique_id;
            for (int k = 0; k < next_node->num_input_pins && !found; k++) {
                for (int l = 0; l < next_node->input_pins[k]->net->num_driver_pins; l++) {
                    if (id == next_node->input_pins[k]->net->driver_pins[l]->unique_id) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found)
                return -1;
        }
    }

    return 1;
}

/*---------------------------------------------------------------------------------------------
 * connect adder type output pin to a node
 *-------------------------------------------------------------------------------------------*/
static void connect_output_pin_to_node(int *width, int current_pin, int output_pin_id, nnode_t *node, nnode_t *current_adder, short subtraction)
{
    // output
    if (subtraction) {
        remap_pin_to_new_node(node->output_pins[current_pin], current_adder, output_pin_id);
    } else {
        npin_t *node_pin_select =
          node->output_pins[(node->num_input_port_sizes == 2) ? current_pin : (current_pin < width[output_pin_id] - 1) ? current_pin + 1 : 0];
        if (node_pin_select) {
            if (node_pin_select->type != NO_ID || (node->num_input_port_sizes == 2)) {
                remap_pin_to_new_node(node_pin_select, current_adder, output_pin_id);
            } else {
                current_adder->output_pins[output_pin_id] = allocate_npin();
                current_adder->output_pins[output_pin_id]->name = append_string("", "%s~dummy_output~%d", current_adder->name, output_pin_id);
            }
        }
    }
}

/*---------------------------------------------------------------------------------------------
 * make a single half-adder (can do unary subtraction, binary subtraction and addition)
 *-------------------------------------------------------------------------------------------*/
static nnode_t *make_adder(operation_list funct, nnode_t *current_adder, nnode_t *previous_carry, int *width, int current_pin, netlist_t *netlist,
                           nnode_t *node, short subtraction, short mark)
{
    // make a 2 bit 0r 3 bit sum or carry based on previous carry
    nnode_t *new_funct = NULL;
    short is_three_port_gate = 0;

    if (previous_carry == netlist->gnd_node) {
        if (funct == ADDER_FUNC)
            new_funct = make_2port_gate(LOGICAL_XOR, 1, 1, 1, node, mark);
        else if (funct == CARRY_FUNC)
            new_funct = make_2port_gate(LOGICAL_AND, 1, 1, 1, node, mark);
    } else if (previous_carry == netlist->vcc_node) {
        if (funct == ADDER_FUNC)
            new_funct = make_2port_gate(LOGICAL_XNOR, 1, 1, 1, node, mark);
        else if (funct == CARRY_FUNC)
            new_funct = make_2port_gate(LOGICAL_OR, 1, 1, 1, node, mark);
    } else {
        new_funct = make_3port_gate(funct, 1, 1, 1, 1, node, mark);
        connect_nodes(previous_carry, 0, new_funct, 0);
        is_three_port_gate = 1;
    }

    // copy the input pin of a half-adder to another function (CARRY or ADDER)
    if (current_adder != NULL) {
        add_input_pin_to_node(new_funct, copy_input_npin(current_adder->input_pins[0 + is_three_port_gate]), 0 + is_three_port_gate);
        add_input_pin_to_node(new_funct, copy_input_npin(current_adder->input_pins[1 + is_three_port_gate]), 1 + is_three_port_gate);
    }
    // create one from scratch
    else {
        // connect input a
        if (current_pin < width[1]) {
            npin_t *temp_pin = node->input_pins[current_pin];
            oassert(temp_pin->net->num_driver_pins <= 1);
            if (!temp_pin->net->num_driver_pins || temp_pin->net->driver_pins[0]->node->type == GND_NODE) {
                connect_nodes(netlist->gnd_node, 0, new_funct, 0 + is_three_port_gate);
                remove_fanout_pins_from_net(temp_pin->net, temp_pin, temp_pin->pin_net_idx);
            } else if (temp_pin->net->driver_pins[0]->node->type == VCC_NODE) {
                connect_nodes(netlist->vcc_node, 0, new_funct, 0 + is_three_port_gate);
                remove_fanout_pins_from_net(temp_pin->net, temp_pin, temp_pin->pin_net_idx);
            } else {
                remap_pin_to_new_node(temp_pin, new_funct, 0 + is_three_port_gate);
            }
        } else {
            connect_nodes(netlist->gnd_node, 0, new_funct, 0 + is_three_port_gate);
        }

        // connect input b
        if (current_pin < width[2]) {
            // pin a is neighbor to pin b
            npin_t *temp_pin = node->input_pins[current_pin + width[1]];
            oassert(temp_pin->net->num_driver_pins <= 1);
            if (temp_pin->net->num_driver_pins == 0 || temp_pin->net->driver_pins[0]->node->type == GND_NODE) {
                nnode_t *attach_to = (subtraction) ? netlist->vcc_node : netlist->gnd_node;
                connect_nodes(attach_to, 0, new_funct, 1 + is_three_port_gate);
                remove_fanout_pins_from_net(temp_pin->net, temp_pin, temp_pin->pin_net_idx);
            } else if (temp_pin->net->driver_pins[0]->node->type == VCC_NODE) {
                nnode_t *attach_to = (subtraction) ? netlist->gnd_node : netlist->vcc_node;
                connect_nodes(attach_to, 0, new_funct, 1 + is_three_port_gate);
                remove_fanout_pins_from_net(temp_pin->net, temp_pin, temp_pin->pin_net_idx);
            } else {
                if (subtraction) {
                    nnode_t *new_not_cells = make_not_gate(node, mark);
                    remap_pin_to_new_node(temp_pin, new_not_cells, 0);
                    connect_nodes(new_not_cells, 0, new_funct, 1 + is_three_port_gate);
                } else {
                    remap_pin_to_new_node(temp_pin, new_funct, 1 + is_three_port_gate);
                }
            }
        } else {
            nnode_t *attach_to = (subtraction) ? netlist->vcc_node : netlist->gnd_node;
            connect_nodes(attach_to, 0, new_funct, 1 + is_three_port_gate);
        }
    }
    return new_funct;
}

void instantiate_add_w_carry_block(int *width, nnode_t *node, short mark, netlist_t *netlist, short subtraction)
{
    nnode_t *previous_carry = (subtraction) ? netlist->vcc_node : netlist->gnd_node;

    for (int i = 0; i < width[0]; i++) {
        /* set of flags for building purposes */
        short construct_last_carry_flag = (i != width[0] - 1 || !subtraction) ? 1 : 0;

        // build Ripple Carry Adder
        nnode_t *current_adder = make_adder(ADDER_FUNC, NULL, previous_carry, width, i, netlist, node, subtraction, mark);
        if (construct_last_carry_flag)
            previous_carry = make_adder(CARRY_FUNC, current_adder, previous_carry, width, i, netlist, node, subtraction, mark);

        connect_output_pin_to_node(width, i, 0, node, current_adder, subtraction);
    }
}

/**
 * -------------------------------------------------------------------------
 * (function: cleanup_add_old_node)
 *
 * @brief <clean up nodeo, a high level ADD node>
 * In split_adder function, nodeo is splitted to small adders,
 * while because of the complexity of input pin connections they have not been
 * remapped to new nodes, they just copied and added to new nodes. This function
 * will detach input pins from the nodeo. Moreover, it will connect the net of
 * unconnected output signals to the GND node, detach the pin from nodeo and
 * free the output pins to avoid memory leak.
 *
 * @param nodeo representing the old adder node
 * @param netlist representing the current netlist
 *-----------------------------------------------------------------------*/
static void cleanup_add_old_node(nnode_t *nodeo, netlist_t *netlist)
{
    /* Disconnecting input pins from the old node side */
    for (int i = 0; i < nodeo->num_input_pins; i++) {
        nodeo->input_pins[i] = NULL;
    }

    /* connecting the extra output pins to the gnd node */
    for (int i = 0; i < nodeo->num_output_pins; i++) {
        npin_t *output_pin = nodeo->output_pins[i];

        if (output_pin && output_pin->node) {
            /* for now we just pass the signals directly through */
            npin_t *zero_pin = get_zero_pin(netlist);
            int idx_2_buffer = zero_pin->pin_net_idx;

            // Dont eliminate the buffer if there are multiple drivers or the AST included it
            if (output_pin->net->num_driver_pins <= 1) {
                /* join all fanouts of the output net with the input pins net */
                join_nets(zero_pin->net, output_pin->net);

                /* erase the pointer to this buffer */
                zero_pin->net->fanout_pins[idx_2_buffer] = NULL;
            }

            free_npin(zero_pin);
            free_npin(output_pin);

            /* Disconnecting output pins from the old node side */
            nodeo->output_pins[i] = NULL;
        }
    }

    // CLEAN UP
    free_nnode(nodeo);
}

/**
 *-------------------------------------------------------------------------------------------
 * (function: check_missing_ports )
 *
 * @brief check for missing ports such as carry-in/out in case of
 * dealing with generated netlist from Yosys blif file.
 *
 * @param node pointing to the netlist node
 * @param traverse_mark_number unique traversal mark for blif elaboration pass
 * @param netlist pointer to the current netlist file
 *-----------------------------------------------------------------------------------------*/
nnode_t *check_missing_ports(nnode_t *node, uintptr_t traverse_mark_number, netlist_t *netlist)
{
    nnode_t *new_node = NULL;
    int num_input_port = node->num_input_port_sizes;

    /* check for operations that has 2 operands */
    if (num_input_port == 2) {
        int in_port1_size = node->input_port_sizes[0];
        int in_port2_size = node->input_port_sizes[1];
        int out_port_size = (in_port1_size >= in_port2_size) ? in_port1_size + 1 : in_port2_size + 1;

        new_node = make_3port_gate(node->type, in_port1_size, in_port2_size, 1, out_port_size, node, traverse_mark_number);

        /* copy attributes */
        copy_attribute(new_node->attributes, node->attributes);

        for (int i = 0; i < in_port1_size; i++) {
            remap_pin_to_new_node(node->input_pins[i], new_node, i);
        }

        for (int i = 0; i < in_port2_size; i++) {
            remap_pin_to_new_node(node->input_pins[i + in_port1_size], new_node, i + in_port1_size);
        }

        /* adding a cin connected to GND */
        npin_t *cin_pin = get_zero_pin(netlist);
        cin_pin->type = INPUT;
        cin_pin->mapping = vtr::strdup("cin");

        add_input_pin_to_node(new_node, cin_pin, new_node->num_input_pins - 1);

        // moving the output pins to the new node
        for (int i = 0; i < out_port_size; i++) {
            if (i < node->num_output_pins) {
                remap_pin_to_new_node(node->output_pins[i], new_node, i);
            } else {
                npin_t *new_pin1 = allocate_npin();
                npin_t *new_pin2 = allocate_npin();
                nnet_t *new_net = allocate_nnet();
                new_net->name = make_full_ref_name(NULL, NULL, NULL, new_node->name, i);
                /* hook the output pin into the node */
                add_output_pin_to_node(new_node, new_pin1, i);
                /* hook up new pin 1 into the new net */
                add_driver_pin_to_net(new_net, new_pin1);
                /* hook up the new pin 2 to this new net */
                add_fanout_pin_to_net(new_net, new_pin2);
            }
        }

        /**
         * if number of output pins is greater than the max of input pins,
         * here we connect the exceeded pins to the GND
         */
        for (int i = out_port_size; i < node->num_output_pins; i++) {
            /* creating a buf node */
            nnode_t *buf_node = make_1port_gate(BUF_NODE, 1, 1, node, traverse_mark_number);
            /* adding the GND input pin to the buf node */
            add_input_pin_to_node(buf_node, get_zero_pin(netlist), 0);
            /* remapping the outpin to buf node */
            remap_pin_to_new_node(node->output_pins[i], buf_node, 0);
        }

        // CLEAN UP
        free_nnode(node);
    }
    /* otherwise there is unary minus, like -A. no need for any change */
    else if (num_input_port == 1) {
        new_node = node;
    }

    return new_node;
}