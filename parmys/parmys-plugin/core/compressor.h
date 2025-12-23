/*
* Added by Junius Pun (juniuspun00@gmail.com)
* 
* Implementation of compressor trees to shrink multi-level additions into Boolean logic that can be packed into LUTs.
*/

#ifndef _COMPRESSOR_H_
#define _COMPRESSOR_H_

#include "odin_types.h"

enum class compressor_tree_type_e {
    WALLACE,             // Wallace tree - reduces to max rank height 2, then binary adder
    WALLACE_TERNARY,     // Wallace tree for ternary adders - reduces to max rank height 3, then ternary adder chain
    WALLACE_TERNARY_EXP, // Experimental: Wallace ternary with HA preference to preserve height 3
    DADDA,               // Dadda tree
    CASCADE,             // Cascade-friendly sequential accumulation for double-carry-chain architectures
    TERNARY_TREE         // Ternary adder tree for DCC3 chain topology - groups 3 inputs using sumout→input chaining
};

// Implement the compressor tree according to the specified type.
extern signal_list_t *implement_compressor_tree(compressor_tree_type_e tree_type, nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);

#endif // _COMPRESSOR_H_