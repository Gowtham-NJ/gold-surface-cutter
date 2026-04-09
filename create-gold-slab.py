#!/usr/bin/python3

import math
import random

# 5 5x6 layers: 14.7 x 15.2 A  (150 atoms)
# 4 7x8 layers: 20.5 x 20.3 A  (224 atoms)
CellReplication = [18, 25, 2 ]

# Au-Au distance 2.93 A
LatticeConst = math.sqrt(2.0)*4.1436457/2.0

# Z-dimension of the slab model
CellZDimension = 100.0

################################################################################

# Coordinates of atoms in one cell - now BUFFERS instead of writing directly
def cell_coord(f_xyz, buf_aus, buf_aui, buf_aub, v, surf):
  hy = math.sqrt(3.0)/2.0
  hz = math.sqrt(2.0/3.0)
  # basis vector
  a = [ 0.0, 0.0, 0.0 ]
  if ((v[2])%2!=0):
    if ((v[1])%2!=0):
      a[0] = v[0]
    else:
      a[0] = v[0]+0.5
    a[1] = (v[1]-1.0/3.0)*hy
  else:
    if ((v[1])%2!=0):
      a[0] = v[0]+0.5
    else:
      a[0] = v[0]
    a[1] = v[1]*hy
  a[2] = v[2]*hz
  # Atoms
  for i in range(3):
    a[i] *= LatticeConst
  # Q-sites
  q = [ 0.0, 0.0, 0.0 ]
  phi = random.uniform(0.0,2*math.pi)
  tht = random.uniform(0.0,math.pi)
  qr = 0.7
  q[0] = a[0] + qr*math.cos(phi)*math.sin(tht)
  q[1] = a[1] + qr*math.sin(phi)*math.sin(tht)
  q[2] = a[2] + qr*math.cos(tht)
  # Save XYZ (unchanged order)
  f_xyz.write("Au {:12.6f}{:12.6f}{:12.6f}\n".format(a[0],a[1],a[2]))
  # Buffer for reordering
  if (surf!=0):
    # AU + AUC go into buf_aus
    buf_aus.append(('AU', 'AUS', a[:], 'atom'))
    buf_aus.append(('AUC','AUS', q[:], 'charge'))
    # Virtual sites go into buf_aui
    r = 2.93
    x1 = [a[0], a[1] + r*math.sqrt(3.0)/3.0, a[2]]
    x2 = [a[0], a[1] - r*math.sqrt(3.0)/3.0, a[2]]
    buf_aui.append(('AUI','AUI', x1[:]))
    buf_aui.append(('AUI','AUI', x2[:]))
  else:
    # Bulk AU + AUC go into buf_aub
    buf_aub.append(('AU', 'AUB', a[:], 'atom'))
    buf_aub.append(('AUC','AUB', q[:], 'charge'))

# Main function
def main():
  ly = math.sqrt(3.0)/2.0*LatticeConst
  lz = math.sqrt(2.0/3.0)*LatticeConst
  # Number of atoms / Q sites / virtual sites
  na = CellReplication[0]*CellReplication[1]*CellReplication[2]
  nq = na
  nx = 4*CellReplication[0]*CellReplication[1]
  # Cell size
  c = [ 0.0, 0.0, 0.0 ]
  c[0] = CellReplication[0]*LatticeConst
  c[1] = CellReplication[1]*ly
  c[2] = CellZDimension
  # Output files
  f_xyz = open('gold-slab.xyz','w')
  f_xyz.write("{:d}\n".format(na))
  f_xyz.write("Box: {:12.6f}{:12.6f}{:12.6f}\n".format(c[0],c[1],c[2]))

  # Buffers for reordering
  buf_aus = []  # Surface AU + AUC
  buf_aui = []  # Virtual sites AUI
  buf_aub = []  # Bulk AU + AUC

  # Cell origin
  orig = [ 0, 0, 0 ]
  # Z direction
  for k in range(CellReplication[2]):
    orig[2] = k
    surf = 0
    if (k==0 or (k+1)==CellReplication[2]):
      surf = 1
    # X-Y layers
    for j in range(CellReplication[1]):
      orig[1] = j
      for i in range(CellReplication[0]):
        orig[0] = i
        cell_coord(f_xyz, buf_aus, buf_aui, buf_aub, orig, surf)

  f_xyz.close()

  # Now write GRO, PDB, TOP in desired order: AUS -> AUI -> AUB
  f_gro = open('gold-slab.gro','w')
  f_gro.write('Gold slab\n')
  f_gro.write("{:d}\n".format(na+nq+nx))
  f_pdb = open('gold-slab.pdb','w')
  f_pdb.write('{:8s}{:s}\n'.format('HEADER','Gold slab'))
  f_top = open('gold-slab.top','w')
  f_top.write('; Force field parameters\n')
  f_top.write('#include "golp-charmm-27.ff/forcefield.itp"\n\n')
  f_top.write('; Gold topologies\n')
  f_top.write('#include "golp-charmm-27.ff/gold_surface.itp"\n')
  f_top.write('#include "golp-charmm-27.ff/gold_bulk.itp"\n')
  f_top.write('#include "golp-charmm-27.ff/gold_virtual.itp"\n\n')
  f_top.write('[ system ]\n')
  f_top.write('; Name\n')
  f_top.write('Gold slab\n\n')
  f_top.write('[ molecules ]\n')
  f_top.write('; Compound        #mols\n')

  atom_id = 1
  res_id = 1

  # 1) AUS: surface AU + AUC (written in pairs)
  n_aus_pairs = len(buf_aus) // 2
  for i in range(0, len(buf_aus), 2):
    au  = buf_aus[i]    # ('AU', 'AUS', coords, 'atom')
    auc = buf_aus[i+1]  # ('AUC','AUS', coords, 'charge')
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
      res_id,'AUS','AU',atom_id,au[2][0]/10.0,au[2][1]/10.0,au[2][2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n".format(
      'HETATM',atom_id,'AU','AUS',res_id,au[2][0],au[2][1],au[2][2],0.0))
    atom_id += 1
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
      res_id,'AUS','AUC',atom_id,auc[2][0]/10.0,auc[2][1]/10.0,auc[2][2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n".format(
      'HETATM',atom_id,'AUC','AUS',res_id,auc[2][0],auc[2][1],auc[2][2],0.0))
    atom_id += 1
    res_id += 1

  # 2) AUI: virtual sites (written in pairs per original surface atom)
  n_aui_pairs = len(buf_aui) // 2
  for v in buf_aui:
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
      res_id,'AUI','AUI',atom_id,v[2][0]/10.0,v[2][1]/10.0,v[2][2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n".format(
      'HETATM',atom_id,'AUI','AUI',res_id,v[2][0],v[2][1],v[2][2],0.0))
    atom_id += 1
    res_id += 1

  # 3) AUB: bulk AU + AUC (written in pairs)
  for i in range(0, len(buf_aub), 2):
    au  = buf_aub[i]
    auc = buf_aub[i+1]
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
      res_id,'AUB','AU',atom_id,au[2][0]/10.0,au[2][1]/10.0,au[2][2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n".format(
      'HETATM',atom_id,'AU','AUB',res_id,au[2][0],au[2][1],au[2][2],0.0))
    atom_id += 1
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
      res_id,'AUB','AUC',atom_id,auc[2][0]/10.0,auc[2][1]/10.0,auc[2][2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n".format(
      'HETATM',atom_id,'AUC','AUB',res_id,auc[2][0],auc[2][1],auc[2][2],0.0))
    atom_id += 1
    res_id += 1

  # Topology: three summary lines
  f_top.write("GoldSurface       {:d}\n".format(len(buf_aus)//2))
  f_top.write("GoldVirtualSite   {:d}\n".format(len(buf_aui)))
  f_top.write("GoldBulk          {:d}\n".format(len(buf_aub)//2))

  # Cell size
  f_gro.write("{:10.5f} {:10.5f} {:10.5f}\n".format(c[0]/10.0,c[1]/10.0,c[2]/10.0))
  # Close files
  f_gro.close()
  f_pdb.close()
  f_top.close()

################################################################################

if __name__=='__main__':
  main()

################################################################################
