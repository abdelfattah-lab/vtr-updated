/*
* Added by Junius Pun (juniuspun00@gmail.com)
* 
* Implementation of compressor trees to shrink multi-level additions into Boolean logic that can be packed into LUTs.
*/

#ifndef _COMPRESSOR_H_
#define _COMPRESSOR_H_

#include "odin_types.h"

// Wallace compressor tree implementation.
extern signal_list_t *implement_compressor_tree_wallace(nnode_t *node, short mark, netlist_t *netlist, std::vector<std::vector<npin_t *>> ranks);

#endif // _COMPRESSOR_H_