# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS:
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines
from ase.lattice.cubic import Diamond

# set of utility routines specific this this model/testing framework
from utilities import relax_config

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
fmax = 0.01 # maximum force following relaxtion [eV/A]

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0)

# specify that we will use model.calculator to compute forces, energies and stresses
bulk.set_calculator(model.calculator)

# use one of the routines from utilities module to relax the initial
# unit cell and atomic positions
bulk = relax_config(bulk, relax_pos=True, relax_cell=False, tol=fmax, traj_file=None)

# set up supercell
bulk *= (5, 1, 1)

def surface_energy(bulk, opening):
    Nat = bulk.get_number_of_atoms()

    # relax atom positions, holding cell fixed
    # vac = relax_atoms(vac, fmax=fmax)

    # compute surface formation energy as difference of bulk and expanded cell
    ebulk = bulk.get_potential_energy()
    print('bulk cell energy', ebulk)

    bulk.cell[0,:] += [opening,0.0,0.0]
    eexp  = bulk.get_potential_energy()

    print('expanded cell energy', eexp)
    e_form = (eexp - ebulk) / (bulk.cell[1,1]*bulk.cell[2,2])
    print('unrelaxed 100 surface formation energy', e_form)
    return e_form

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
n_steps = 35
max_opening = 3.5

al = []

openings = []
es = []
for i in range(n_steps + 1):
    opening = float(i)/float(n_steps)*max_opening
    openings.append(opening)
    bulk_copy = bulk.copy()
    bulk_copy.set_calculator(model.calculator)
    al.append(bulk_copy)
    es.append(surface_energy(bulk_copy, opening))

from ase.io import read, write

write("./decoh_traj.xyz", al)

print("openings ", openings)
print("es ", es)
from scipy import interpolate
spline = interpolate.splrep(openings, es, s=0)
stresses = [x for x in interpolate.splev(openings, spline, der=1)]

print("stresses ", stresses)
properties = {'surface_decohesion_unrelaxed_opening': openings, 'surface_decohesion_unrelaxed_energy' : es, 'surface_decohesion_unrelaxed_stress' : stresses}
