from ase.lattice.cubic import Diamond
import numpy as np

import ase.io, sys

# set of utility routines specific this this model/testing framework
from utilities import relax_config

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
fmax = 0.01 # maximum force following relaxtion [eV/A]

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0, directions=[[1,-1,0],[1,0,-1],[1,1,1]])

# specify that we will use model.calculator to compute forces, energies and stresses
bulk.set_calculator(model.calculator)

# use one of the routines from utilities module to relax the initial
# unit cell and atomic positions
bulk = relax_config(bulk, relax_pos=True, relax_cell=False, tol=fmax, traj_file=None)

# set up supercell
bulk *= (1, 1, 5)


def surface_energy(bulk, z_offset, opening):
    Nat = bulk.get_number_of_atoms()

    # shift so cut is through shuffle plane
    bulk.positions[:,2] += z_offset
    bulk.wrap()

    # relax atom positions, holding cell fixed
    # vac = relax_atoms(vac, fmax=fmax)

    # compute surface formation energy as difference of bulk and expanded cell
    ebulk = bulk.get_potential_energy()
    print('bulk cell energy', ebulk)

    bulk.cell[2,:] += [0.0,0.0,opening]
    eexp  = bulk.get_potential_energy()

    ase.io.write(sys.stdout, bulk, format='extxyz')

    print('expanded cell energy', eexp)
    e_form = 0.5*(eexp - ebulk) / np.linalg.norm(np.cross(bulk.cell[0,:],bulk.cell[1,:]))
    print('unrelaxed 111 surface formation energy', e_form)
    return e_form

# dictionary of computed properties - this is output of this test, to
#   be compared with other models

n_steps = 30
max_opening = 3.0

openings = []
es = []
for i in range(n_steps + 1):
    opening = float(i)/float(n_steps)*max_opening
    openings.append(opening)
    bulk_copy = bulk.copy()
    bulk_copy.set_calculator(model.calculator)
    es.append(surface_energy(bulk_copy, 2.0, opening))

print("openings ", openings)
print("es ", es)
from scipy import interpolate
spline = interpolate.splrep(openings, es, s=0)
stresses = [x for x in interpolate.splev(openings, spline, der=1)]

print("stresses ", stresses)
properties = {'surface_decohesion_unrelaxed_opening': openings, 'surface_decohesion_unrelaxed_energy' : es, 'surface_decohesion_unrelaxed_stress' : stresses}
