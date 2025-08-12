# Gold Surface Cutter

A Python toolkit for extracting proteins from large gold junction systems and placing them on optimally-sized Au(111) surfaces with controlled padding and registry preservation.

## Overview

This tool addresses a common problem in computational chemistry: when you have a large protein-gold junction system and need to create a smaller, more computationally manageable system while maintaining proper surface registry and adequate spacing for simulations.

**Key Features:**
- 🔄 Preserves crystallographic registry (protein shifted by integer lattice vectors only)
- 📏 Guarantees minimum padding around the protein
- 🎯 Optimizes protein placement to center it on the new slab
- ⚡ Maintains slab origin for DFT calculations
- 📊 Reports final gap distances for verification

## Files

- **`cut.py`** - Main script that analyzes protein dimensions, generates new slab, and places protein
- **`create-gold-slab.py`** - Gold slab generator with Au(111) surface structure

## Requirements

- Python 3.6+
- NumPy
- Standard library modules: `math`, `argparse`, `pathlib`, `importlib`

## Usage

### Basic Usage

```bash
python cut.py old_system.gro
```

### With Custom Padding

```bash
python cut.py old_system.gro --pad 1.5
```

### Parameters

- **`old_system.gro`** - Input GROMACS structure file containing protein + gold surface
- **`--pad`** - Padding distance in nanometers around protein (default: 1.0 nm = 10 Å)

### Output Files

- **`gold-slab.gro`** - Regenerated bare gold slab
- **`gold-slab.xyz`** - XYZ format of the slab  
- **`gold-slab.pdb`** - PDB format of the slab
- **`gold-slab.top`** - GROMACS topology file
- **`protein_on_new_slab.gro`** - Final merged system ready for MD/DFT

## How It Works

### 1. System Analysis
- Reads input `.gro` file and separates gold atoms (AU/AUI/AUC) from protein/solvent
- Calculates protein bounding box dimensions

### 2. Slab Sizing
- Determines minimum slab size needed: `protein_size + 2 × padding`
- Rounds to even number of unit cells for symmetry
- Uses Au(111) lattice parameters: a = 0.293 nm, ly = 0.254 nm

### 3. Registry-Preserving Translation
- Finds all possible integer lattice vector shifts `(m,n)` that fit the protein
- Selects the shift that places protein center-of-mass closest to slab center
- Ensures crystallographic registry is maintained

### 4. Quality Control
- Reports final gaps on all four sides
- Verifies minimum padding requirements are met
- Preserves all atom indices and residue information

## Au(111) Surface Model

The gold surface uses a realistic Au(111) model with:
- **Lattice constant**: 2.93 Å (√2 × 4.144 Å / 2)
- **Surface structure**: FCC(111) with proper layer stacking
- **Virtual sites**: AUI atoms for improved force field accuracy
- **Charge sites**: AUC atoms for electrostatic interactions

## Example Output

```
→ building Au slab: 20 × 26 unit cells  (≥ 15 Å margin)
Protein shift : m = 3, n = -2  →  Δx = 0.732 nm, Δy = -0.461 nm
Gaps (nm)     : left 1.234 | right 1.456 | bottom 1.123 | top 1.334
✓  protein_on_new_slab.gro written
   Box: 5.864 × 6.602 × 10.000 nm
```

## Input File Requirements

Your input `.gro` file should contain:
- **Gold atoms**: Named AU, AUI, or AUC (these define the "slab" and are never moved)
- **Protein/solvent**: All other atom types (these get shifted as a rigid body)

## Technical Details

### Lattice Constants
- **a**: 0.293195 nm (√2 × 4.1436457 Å / 2 / 10)
- **ly**: 0.253934 nm (√3/2 × a)
- **Surface**: Au(111) with (√3 × √3)R30° reconstruction capability

### Registry Preservation
The algorithm ensures the protein is shifted by exactly:
```
Δx = m × a + 0.5 × n × a
Δy = n × ly  
```
where `m` and `n` are integers, preserving the crystallographic relationship.

### Error Handling
- Validates that `create-gold-slab.py` exists in the same directory
- Checks that at least one valid translation exists
- Reports if padding requirements cannot be met
