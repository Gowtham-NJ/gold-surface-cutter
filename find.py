#!/usr/bin/python3

import sys
import math
import argparse

def read_ndx(ndx_file):
    """Read GROMACS index file and return dict of group_name -> list of atom indices."""
    groups = {}
    current = None
    with open(ndx_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('['):
                current = line.strip('[] ')
                groups[current] = []
            elif current is not None:
                groups[current].extend(int(x) for x in line.split())
    return groups

def read_gro(gro_file):
    """Read GRO file and return dict of atom_index -> (x, y, z) in nm."""
    coords = {}
    with open(gro_file, 'r') as f:
        title = f.readline()
        natoms = int(f.readline().strip())
        for i in range(natoms):
            line = f.readline()
            atom_id = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            coords[atom_id] = (x, y, z)
    return coords

def distance(c1, c2):
    return math.sqrt(sum((a - b)**2 for a, b in zip(c1, c2)))

def main():
    parser = argparse.ArgumentParser(description='Find nearest 6 neighbours to a target index group.')
    parser.add_argument('-n', '--ndx', required=True, help='GROMACS index file (.ndx)')
    parser.add_argument('-f', '--gro', required=True, help='GROMACS coordinate file (.gro)')
    args = parser.parse_args()

    # Read files
    groups = read_ndx(args.ndx)
    coords = read_gro(args.gro)

    # List available groups
    print("\nAvailable index groups:")
    print("-" * 40)
    for i, (name, atoms) in enumerate(groups.items()):
        print(f"  {i:3d} : {name:20s} ({len(atoms)} atoms)")
    print()

    # Get target group
    target_name = input("Enter the TARGET index group name: ").strip()
    if target_name not in groups:
        print(f"Error: group '{target_name}' not found.")
        sys.exit(1)

    # Get neighbour search group
    search_name = input("Enter the NEIGHBOUR search group name: ").strip()
    if search_name not in groups:
        print(f"Error: group '{search_name}' not found.")
        sys.exit(1)

    target_atoms = groups[target_name]
    search_atoms = groups[search_name]

    print(f"\nFinding 6 nearest neighbours from [{search_name}] for each atom in [{target_name}]")
    print("=" * 80)

    for t_id in target_atoms:
        if t_id not in coords:
            print(f"Warning: atom {t_id} not found in GRO file, skipping.")
            continue
        tc = coords[t_id]

        # Calculate distances to all search atoms (excluding self)
        dists = []
        for s_id in search_atoms:
            if s_id == t_id:
                continue
            if s_id not in coords:
                continue
            d = distance(tc, coords[s_id])
            dists.append((s_id, d))

        # Sort and pick 6 nearest
        dists.sort(key=lambda x: x[1])
        nearest = dists[:6]

        print(f"\nAtom {t_id:6d}  ({tc[0]:8.3f} {tc[1]:8.3f} {tc[2]:8.3f}) nm")
        print(f"  {'Rank':>4s}  {'AtomID':>8s}  {'Distance (nm)':>14s}  {'Distance (A)':>13s}")
        for rank, (s_id, d) in enumerate(nearest, 1):
            sc = coords[s_id]
            print(f"  {rank:4d}  {s_id:8d}  {d:14.4f}  {d*10:13.4f}")

if __name__ == '__main__':
    main()
