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

# Coordinates of atoms in one cell
#
# f_xyz - output XYZ file
# f_gro - output GRO file
# f_pdb - output PDB file
# f_top - output TOP file
# v     - position of the cell
# surf  - surface layer
# ar    - residuum / atom ID
def cell_coord(f_xyz,f_gro,f_pdb,f_top,v,surf,ar):
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
  # Save coordinates
  f_xyz.write("Au {:12.6f}{:12.6f}{:12.6f}\n".format(a[0],a[1],a[2]))
  if (surf!=0):
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n". \
      format(ar[1],'AUS','AU',ar[0],a[0]/10.0,a[1]/10.0,a[2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n". \
      format('HETATM',ar[0],'AU','AUS',ar[1],a[0],a[1],a[2],0.0))
    ar[0] = ar[0] + 1
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n". \
      format(ar[1],'AUS','AUC',ar[0],q[0]/10.0,q[1]/10.0,q[2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n". \
      format('HETATM',ar[0],'AUC','AUS',ar[1],q[0],q[1],q[2],0.0))
    ar[0] = ar[0] + 1
    ar[1] = ar[1] + 1
    f_top.write("GoldSurface       1\n");
    # Virtual sites
    x = [ 0.0, 0.0, 0.0 ]
    r = 2.93
    x[0] = a[0]
    x[1] = a[1] + r*math.sqrt(3.0)/3.0
    x[2] = a[2]
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n". \
      format(ar[1],'AUI','AUI',ar[0],x[0]/10.0,x[1]/10.0,x[2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n". \
      format('HETATM',ar[0],'AUI','AUI',ar[1],x[0],x[1],x[2],0.0))
    ar[0] = ar[0] + 1
    ar[1] = ar[1] + 1
    x[0] = a[0]
    x[1] = a[1] - r*math.sqrt(3.0)/3.0
    x[2] = a[2]
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n". \
      format(ar[1],'AUI','AUI',ar[0],x[0]/10.0,x[1]/10.0,x[2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n". \
      format('HETATM',ar[0],'AUI','AUI',ar[1],x[0],x[1],x[2],0.0))
    ar[0] = ar[0] + 1
    ar[1] = ar[1] + 1
    f_top.write("GoldVirtualSite   2\n");
  else:
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n". \
      format(ar[1],'AUB','AU',ar[0],a[0]/10.0,a[1]/10.0,a[2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n". \
      format('HETATM',ar[0],'AU','AUB',ar[1],a[0],a[1],a[2],0.0))
    ar[0] = ar[0] + 1
    f_gro.write("{:5d}{:7s}{:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n". \
      format(ar[1],'AUB','AUC',ar[0],q[0]/10.0,q[1]/10.0,q[2]/10.0))
    f_pdb.write("{:8s}{:3d} {:4s} {:3s}{:6d}    {:8.3f}{:8.3f}{:8.3f}{:8.4f}\n". \
      format('HETATM',ar[0],'AUC','AUB',ar[1],q[0],q[1],q[2],0.0))
    ar[0] = ar[0] + 1
    ar[1] = ar[1] + 1
    f_top.write("GoldBulk          1\n");

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
#  c[2] = CellReplication[2]*lz
  c[2] = CellZDimension
  # Output files
  f_xyz = open('gold-slab.xyz','w')
  f_xyz.write("{:d}\n".format(na))
  f_xyz.write("Box: {:12.6f}{:12.6f}{:12.6f}\n".format(c[0],c[1],c[2]))
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
  # Cell origin
  orig = [ 0, 0, 0 ]
  ar = [ 1, 1 ]
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
        # Save coordinates
        cell_coord(f_xyz,f_gro,f_pdb,f_top,orig,surf,ar)
  # Cell size
  f_gro.write("{:10.5f} {:10.5f} {:10.5f}\n". \
    format(c[0]/10.0,c[1]/10.0,c[2]/10.0))
  # Close files
  f_xyz.close()
  f_gro.close()
  f_pdb.close()
  f_top.close()

################################################################################

if __name__=='__main__':
  main()

################################################################################
