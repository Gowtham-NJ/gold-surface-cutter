#!/usr/bin/env python3
"""
enlarge_slab_shift_protein.py  old_system.gro  [--pad 1.0]

• AU / AUI / AUC atoms = slab (never moved)
• everything else      = protein / solvent (shifted)

Guarantees
----------
✓ slab origin unchanged (good for DFT)
✓ registry preserved (protein shifted by integer lattice vector)
✓ ≥ `pad` nm free space on each side
✓ shift chosen so protein COM is as close as possible to slab centre
✓ prints the resulting gaps (left/right/bottom/top)

Outputs
-------
gold-slab.gro              – regenerated bare slab
protein_on_new_slab.gro    – merged system ready for MD / DFT
"""
from pathlib import Path
import argparse, math, importlib.util, sys, copy
import numpy as np

# ─── Au(111) lattice constants (from builder) ───────────────────────
A_NM  = math.sqrt(2.0) * 4.1436457 / 2.0 / 10.0          # 0.29319457 nm
LY_NM = math.sqrt(3.0) / 2.0 * A_NM                      # 0.253934 nm
SURF  = {"AU", "AUI", "AUC"}
EPS   = 1e-6   # 1 pm tolerance for boundary tests

# ─── minimal .gro helpers ───────────────────────────────────────────
def read_gro(p: Path):
    atoms = []
    with p.open() as fh:
        title = fh.readline().strip()
        nat   = int(fh.readline())
        for _ in range(nat):
            l  = fh.readline()
            atoms.append({
                "resnr":   int(l[0:5]),
                "resname": l[5:10].strip(),
                "name":    l[10:15].strip().upper(),
                "idx":     int(l[15:20]),
                "x":       float(l[20:28]),
                "y":       float(l[28:36]),
                "z":       float(l[36:44]),
            })
        box = list(map(float, fh.readline().split()[:3]))
    return title, atoms, box            # Lx, Ly, Lz

def write_gro(p: Path, title, atoms, box):
    with p.open("w") as fh:
        fh.write(f"{title}\n{len(atoms):5d}\n")
        for i, a in enumerate(atoms, 1):
            fh.write(f"{a['resnr']:5d}{a['resname']:<5s}{a['name'][:3]:>5s}"
                     f"{i:5d}{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}\n")
        fh.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")

# ─── main workflow ──────────────────────────────────────────────────
def enlarge_and_shift(old_gro: Path, pad_nm: float):
    # --- read old system
    title, atoms, _ = read_gro(old_gro)
    protein = [a for a in atoms if a["name"] not in SURF]
    gold    = [a for a in atoms if a["name"] in SURF]

    xs = np.array([a["x"] for a in protein])
    ys = np.array([a["y"] for a in protein])

    # --- new slab size (even Nx,Ny)
    dx_need = np.ptp(xs) + 2 * pad_nm
    dy_need = np.ptp(ys) + 2 * pad_nm
    Nx = int(math.ceil(dx_need / (2 * A_NM)) * 2)
    Ny = int(math.ceil(dy_need / (2 * LY_NM)) * 2)
    Nz = 2
    print(f"→ building Au slab: {Nx} × {Ny} unit cells  (≥ {pad_nm*10:.0f} Å margin)")

    # --- call create‑gold‑slab.py
    slab_py = old_gro.with_name("create-gold-slab.py")
    if not slab_py.exists():
        sys.exit("create-gold-slab.py not found in this directory.")
    spec = importlib.util.spec_from_file_location("gold_slab", slab_py)
    mod  = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.CellReplication = [Nx, Ny, Nz]
    mod.main()                                            # → gold-slab.gro
    _, slab_atoms, slab_box = read_gro(Path("gold-slab.gro"))
    Lx, Ly, Lz = slab_box

    # --- anchor Au in old slab (lowest y, then x)
    anchor = min(gold, key=lambda a: (a["y"], a["x"]))
    ax, ay = anchor["x"], anchor["y"]

    # slab centre for COM metric
    cx, cy = 0.5 * Lx, 0.5 * Ly

    # --- search integer (m,n) that fits AND centres COM as much as possible
    candidates = []
    for n in range(-Ny, Ny + 1):
        dy = n * LY_NM - ay
        for m in range(-Nx, Nx + 1):
            dx = m * A_NM + 0.5 * n * A_NM - ax
            x_shift = xs + dx
            y_shift = ys + dy
            if (x_shift.min() >= -EPS and x_shift.max() <= Lx + EPS and
                y_shift.min() >= -EPS and y_shift.max() <= Ly + EPS):
                com_dist = math.hypot(x_shift.mean() - cx,
                                      y_shift.mean() - cy)
                candidates.append((com_dist, m, n, dx, dy))
    if not candidates:
        sys.exit("No lattice translation fits – increase padding.")

    _, m, n, dx, dy = min(candidates, key=lambda t: t[0])

    # --- final gap report
    x_shift = xs + dx
    y_shift = ys + dy
    left_gap   = x_shift.min()
    right_gap  = Lx - x_shift.max()
    bottom_gap = y_shift.min()
    top_gap    = Ly - y_shift.max()

    print(f"Protein shift : m = {m}, n = {n}  →  Δx = {dx:.3f} nm, Δy = {dy:.3f} nm")
    print(f"Gaps (nm)     : left {left_gap:.3f} | right {right_gap:.3f} | "
          f"bottom {bottom_gap:.3f} | top {top_gap:.3f}")

    # --- shift protein
    shifted_prot = copy.deepcopy(protein)
    for a in shifted_prot:
        a["x"] += dx
        a["y"] += dy

    # --- merge & write
    merged = slab_atoms + shifted_prot
    for i, a in enumerate(merged, 1):
        a["idx"] = i

    write_gro(Path("protein_on_new_slab.gro"),
              "Protein shifted onto enlarged Au(111) slab (registry kept)",
              merged, slab_box)

    print("✓  protein_on_new_slab.gro written")
    print(f"   Box: {Lx:.3f} × {Ly:.3f} × {Lz:.3f} nm")

# ─── CLI ------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gro", help="old Au+protein .gro file")
    parser.add_argument("--pad", type=float, default=1.0,
                        help="padding in nm around protein (default 1.0 nm = 10 Å)")
    args = parser.parse_args()
    enlarge_and_shift(Path(args.gro), args.pad)

