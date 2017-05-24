import symmetrize
from math import sqrt
from ase import Atoms
from ase.io import write
import sys
import numpy as np

at = Atoms("Si2", cell=[ [ 1, 0, 4.0 ], [-0.5, sqrt(3.0)/2.0, 4.0], [-0.5, -sqrt(3.0)/2.0, 4.0] ], pbc=[True, True, True])
at.positions[0,:] = [0.0, 0.0, 0.5]
at.positions[1,:] = [0.0, 0.0, 2.5]

(r, t, m) = symmetrize.prep(at, 0.01)

forces = np.random.random_sample(at.positions.shape)
forces_sym = symmetrize.forces(at.get_cell(), at.get_reciprocal_cell().T, forces, r, t, m)
print "forces"
for i in range(len(at)):
    print forces[i,:], "   ", forces_sym[i,:]

stress = np.random.random_sample((3,3))
stress_sym = symmetrize.stress(at.get_cell(), at.get_reciprocal_cell().T, stress, r)
print "stress"
for i in range(3):
    print stress[i,:], "   ", stress_sym[i,:]

