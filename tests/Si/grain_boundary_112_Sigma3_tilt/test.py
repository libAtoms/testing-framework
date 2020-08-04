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
import numpy as np

import ase.io, sys, os

# set of utility routines specific this this model/testing framework
from utilities import relax_config

# the current model
import model

np.random.seed(75)

a0 = 5.44
fmax = 0.01 # maximum force following relaxtion [eV/A]

if not hasattr(model, 'bulk_reference'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=a0)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=fmax, traj_file=None)
else:
    bulk = model.bulk_reference.copy()
    bulk.set_calculator(model.bulk_reference.calc)

print("BOB doing bulk pot en")
e0 = bulk.get_potential_energy()/float(len(bulk))
print("got e0 ", e0)


def interface_energy(e0, filename):

    interface = ase.io.read(os.path.join(os.path.dirname(__file__),filename), format='extxyz')

    # adjust lattice constant
    interface.set_calculator(model.calculator)

    # relax atom positions, holding cell fixed
    interface.positions += 0.01*np.random.random_sample((len(interface),3))
    print("pre relax config")
    ase.io.write(sys.stdout, interface, format='extxyz')
    print("BOB doing gb relax")
    interface = relax_config(interface, relax_pos=True, relax_cell=True, tol=fmax, max_steps=300, traj_file=os.path.splitext(filename)[0]+"-relaxed.opt.xyz")
    print("post relax config")
    ase.io.write(sys.stdout, interface, format='extxyz')

    print("BOB doing gb get pot en")
    einterface  = interface.get_potential_energy()
    print('interface relaxed cell energy', einterface)
    e_form = (einterface-len(interface) * e0)
    area = np.linalg.norm(np.cross(interface.get_cell()[0,:],interface.get_cell()[1,:]))
    print("got area",area)
    print('interface relaxed formation energy/area', e_form/(2.0*area))
    return e_form/(2.0*area)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'gb_energy': interface_energy(e0, "gb.0.25_0.0.xyz") }
