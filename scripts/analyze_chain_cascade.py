#!/usr/bin/env python3
"""
Analyze pre_packing_molecules_and_patterns.echo files to identify chain molecules
that have invalid cascade connections (row 1's a/b not connected to row 0's sumout).

For the "chain" pattern in DCC3 architecture:
- Row 0: pattern indices 0-19
- Row 1: pattern indices 20-39 (corresponds to row 0 in reverse: 20<->19, 21<->18, ..., 39<->0)

The architecture requires: row0[i].sumout -> row1[39-i].a or row1[39-i].b
If both a and b are connected to gnd/vcc instead, this is a problem.
"""

import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class PatternBlock:
    """Represents a block at a pattern index in a molecule."""
    index: int
    atom_name: str
    atom_id: int
    row: int  # 0 or 1, -1 for unknown/non-adder
    is_empty: bool = False
    is_adder: bool = False  # True if this is an adder primitive

    # Pin connections: pin_name -> net_name
    pins: Dict[str, str] = field(default_factory=dict)

    # Sumout net name (for row 0 blocks)
    sumout_net: Optional[str] = None


@dataclass
class ChainMolecule:
    """Represents a chain molecule."""
    root_atom: str
    blocks: Dict[int, PatternBlock] = field(default_factory=dict)

    # Mapping from row 1 index to row 0 index (for correspondence)
    row1_to_row0: Dict[int, int] = field(default_factory=dict)

    def build_row_correspondence(self):
        """Build the mapping between row 0 and row 1 adder positions."""
        row0_adders = sorted([idx for idx, b in self.blocks.items()
                              if b.is_adder and b.row == 0 and not b.is_empty])
        row1_adders = sorted([idx for idx, b in self.blocks.items()
                              if b.is_adder and b.row == 1 and not b.is_empty], reverse=True)

        # Row 0 goes forward (0,1,2,...), Row 1 goes backward (39,38,37,...)
        # So row0[0] <-> row1[max], row0[1] <-> row1[max-1], etc.
        for i, (r0_idx, r1_idx) in enumerate(zip(row0_adders, row1_adders)):
            self.row1_to_row0[r1_idx] = r0_idx

    def get_row0_block(self, row1_idx: int) -> Optional[PatternBlock]:
        """Get the corresponding row 0 block for a row 1 index."""
        row0_idx = self.row1_to_row0.get(row1_idx)
        if row0_idx is not None:
            return self.blocks.get(row0_idx)
        return None

    def get_all_row0_sumouts(self) -> Set[str]:
        """Get all sumout net names from row 0 adders."""
        sumouts = set()
        for block in self.blocks.values():
            if block.row == 0 and block.is_adder and not block.is_empty and block.sumout_net:
                sumouts.add(block.sumout_net)
        return sumouts

    def check_cascade_issues(self) -> List[Tuple[int, int, str]]:
        """
        Check for cascade connection issues.
        For each row 1 adder, check if its 'a' or 'b' is connected to ANY row 0 sumout.
        Returns list of (row1_idx, row0_idx, issue_description) tuples.
        """
        issues = []
        row0_sumouts = self.get_all_row0_sumouts()

        for idx, block in self.blocks.items():
            # Only check row 1 adder blocks
            if block.row != 1 or block.is_empty or not block.is_adder:
                continue

            # Skip pass-through blocks (they are created by fill_vacant_chain_spots)
            if 'pass_through' in block.atom_name:
                continue

            # Check if row 1's a or b is connected to ANY row 0's sumout
            a_net = block.pins.get('a[0]', '')
            b_net = block.pins.get('b[0]', '')

            a_connected_to_row0 = (a_net in row0_sumouts)
            b_connected_to_row0 = (b_net in row0_sumouts)

            if not a_connected_to_row0 and not b_connected_to_row0:
                # Check if both are gnd/vcc
                a_is_const = 'gnd' in a_net.lower() or 'vcc' in a_net.lower() or a_net == ''
                b_is_const = 'gnd' in b_net.lower() or 'vcc' in b_net.lower() or b_net == ''

                # Try to find the expected row0 block
                row0_idx = self.row1_to_row0.get(idx, -1)
                expected_sumout = ""
                if row0_idx >= 0 and row0_idx in self.blocks:
                    expected_sumout = self.blocks[row0_idx].sumout_net or ""

                issue = f"row1[{idx}] ({block.atom_name}) a={a_net}, b={b_net}"
                if expected_sumout:
                    issue += f"; expected from row0[{row0_idx}].sumout={expected_sumout}"

                if a_is_const and b_is_const:
                    issue += " [BOTH CONST - CRITICAL]"
                elif a_is_const or b_is_const:
                    issue += " [ONE CONST]"
                else:
                    issue += " [NO ROW0 SUMOUT CONNECTION]"

                issues.append((idx, row0_idx, issue))

        return issues


def parse_echo_file(filepath: str) -> List[ChainMolecule]:
    """Parse the pre_packing_molecules_and_patterns.echo file."""
    molecules = []
    current_molecule = None
    current_block = None

    with open(filepath, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()

        # Start of a new chain molecule
        if line.startswith('molecule type: chain'):
            if current_molecule:
                current_molecule.build_row_correspondence()
                molecules.append(current_molecule)
            current_molecule = ChainMolecule(root_atom="")
            current_block = None

        # Pattern index line
        elif current_molecule and 'pattern index' in line:
            match = re.match(r'\s*pattern index (\d+):\s*(.*)', line)
            if match:
                idx = int(match.group(1))
                rest = match.group(2)

                if 'empty' in rest.lower():
                    current_block = PatternBlock(
                        index=idx,
                        atom_name="",
                        atom_id=-1,
                        row=-1,
                        is_empty=True,
                        is_adder=False
                    )
                else:
                    # Parse atom block name and ID
                    # Format: atom block NAME (ID: NUM) row X
                    atom_match = re.match(r'atom block ([^\s]+)\s+\(ID:\s*(\d+)\)\s*(.*)', rest)
                    if atom_match:
                        atom_name = atom_match.group(1)
                        atom_id = int(atom_match.group(2))
                        rest_info = atom_match.group(3)

                        # Detect row from annotation
                        row = -1
                        if 'row 0' in rest_info:
                            row = 0
                        elif 'row 1' in rest_info:
                            row = 1
                        elif 'root' in rest_info.lower():
                            row = 0  # Root is typically row 0

                        # Detect if this is an adder (has adder-like pins or name)
                        is_adder = ('adder' in rest_info.lower() or
                                   'ADD' in atom_name or
                                   '$add' in atom_name or
                                   'primitive: adder' in rest_info)

                        current_block = PatternBlock(
                            index=idx,
                            atom_name=atom_name,
                            atom_id=atom_id,
                            row=row,
                            is_adder=is_adder
                        )

                        if idx == 0 or 'root' in rest_info.lower():
                            current_molecule.root_atom = atom_name

                if current_block:
                    current_molecule.blocks[idx] = current_block

        # Pin connection line
        elif current_block and not current_block.is_empty and '-> pin' in line:
            # Parse: -> pin name[idx]: net net_name (...)
            pin_match = re.match(r'\s*->\s*pin\s+(\w+\[\d+\]):\s*net\s+(\S+)', line)
            if pin_match:
                pin_name = pin_match.group(1)
                net_name = pin_match.group(2)
                current_block.pins[pin_name] = net_name

                if pin_name == 'sumout[0]':
                    current_block.sumout_net = net_name

                # If we see adder pins, mark as adder
                if pin_name in ('cin[0]', 'cout[0]', 'sumout[0]'):
                    current_block.is_adder = True

        # End of molecule (next molecule type or end of relevant section)
        elif current_molecule and line.startswith('molecule type:') and 'chain' not in line:
            current_molecule.build_row_correspondence()
            molecules.append(current_molecule)
            current_molecule = None
            current_block = None

        i += 1

    # Don't forget the last molecule
    if current_molecule:
        current_molecule.build_row_correspondence()
        molecules.append(current_molecule)

    return molecules


def analyze_molecules(molecules: List[ChainMolecule], verbose: bool = False) -> Dict:
    """Analyze molecules for cascade issues."""
    stats = {
        'total_chain_molecules': len(molecules),
        'molecules_with_issues': 0,
        'total_issues': 0,
        'critical_issues': 0,  # Both a and b are const
        'issues_by_position': defaultdict(int),
        'problematic_molecules': []
    }

    for mol in molecules:
        issues = mol.check_cascade_issues()

        if issues:
            stats['molecules_with_issues'] += 1
            stats['total_issues'] += len(issues)

            critical_count = sum(1 for _, _, desc in issues if 'CRITICAL' in desc)
            stats['critical_issues'] += critical_count

            for row1_idx, row0_idx, desc in issues:
                stats['issues_by_position'][row1_idx] += 1

            stats['problematic_molecules'].append({
                'root': mol.root_atom,
                'num_issues': len(issues),
                'critical': critical_count,
                'details': issues if verbose else issues[:3]  # Limit details unless verbose
            })

    return stats


def print_report(stats: Dict, verbose: bool = False):
    """Print analysis report."""
    print("=" * 80)
    print("CHAIN MOLECULE CASCADE CONNECTION ANALYSIS")
    print("=" * 80)
    print()
    print(f"Total chain molecules analyzed: {stats['total_chain_molecules']}")
    print(f"Molecules with cascade issues: {stats['molecules_with_issues']}")
    print(f"Total cascade issues found: {stats['total_issues']}")
    print(f"Critical issues (both a,b const): {stats['critical_issues']}")
    print()

    if stats['issues_by_position']:
        print("Issues by row 1 position (index 20-39):")
        print("-" * 40)
        for pos in sorted(stats['issues_by_position'].keys()):
            count = stats['issues_by_position'][pos]
            row0_pos = 39 - pos
            print(f"  Row1[{pos}] <-> Row0[{row0_pos}]: {count} issues")
        print()

    if stats['problematic_molecules']:
        print("Problematic molecules:")
        print("-" * 40)
        for i, mol_info in enumerate(stats['problematic_molecules'][:20]):  # Limit to first 20
            print(f"\n[{i+1}] Root: {mol_info['root']}")
            print(f"    Issues: {mol_info['num_issues']} (critical: {mol_info['critical']})")
            if verbose:
                for row1_idx, row0_idx, desc in mol_info['details']:
                    print(f"      - {desc}")

        if len(stats['problematic_molecules']) > 20:
            print(f"\n... and {len(stats['problematic_molecules']) - 20} more molecules")

    print()
    print("=" * 80)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze chain molecule cascade connections in VPR echo files'
    )
    parser.add_argument('filepath', help='Path to pre_packing_molecules_and_patterns.echo file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed issue descriptions')
    parser.add_argument('--json', action='store_true',
                        help='Output as JSON instead of text report')

    args = parser.parse_args()

    print(f"Parsing: {args.filepath}")
    molecules = parse_echo_file(args.filepath)
    print(f"Found {len(molecules)} chain molecules")
    print()

    stats = analyze_molecules(molecules, verbose=args.verbose)

    if args.json:
        import json
        # Convert defaultdict to dict for JSON serialization
        stats['issues_by_position'] = dict(stats['issues_by_position'])
        print(json.dumps(stats, indent=2, default=str))
    else:
        print_report(stats, verbose=args.verbose)


if __name__ == '__main__':
    main()
