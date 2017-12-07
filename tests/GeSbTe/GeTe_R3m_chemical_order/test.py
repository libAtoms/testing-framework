import os.path, ase.io
from utilities import relax_config, run_root
import random

test_dir = os.path.abspath(os.path.dirname(__file__))

bulk = ase.io.read(os.path.join(test_dir,"bulk.xyz"), format="extxyz")
tol=1.0e-3 # hack to work around bad convergence
relaxed_bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, 
    traj_file=None, config_label='bulk', from_base_model=True, save_config=True)
relaxed_bulk_pe = relaxed_bulk.get_potential_energy()/len(bulk)

tol=2.0e-3 # hack to work around convergence issues
bulk *= (2,2,2)

random.seed(10)
unrelaxed_energies = []
relaxed_energies = []
for config_i in range(20):
    Z = bulk.get_atomic_numbers()
    (i1, i2) = random.sample(range(len(bulk)),2)
    while Z[i1] == Z[i2]:
        (i1, i2) = random.sample(range(len(bulk)),2)
    t_Z = Z[i1]
    Z[i1] = Z[i2]
    Z[i2] = t_Z
    bulk.set_atomic_numbers(Z)
    unrelaxed_energies.append(bulk.get_potential_energy()/len(bulk)-relaxed_bulk_pe)
    relaxed_bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, 
        traj_file=None, config_label='supercell_antisite_{}'.format(config_i), from_base_model=True, save_config=True)
    relaxed_energies.append(bulk.get_potential_energy()/len(bulk)-relaxed_bulk_pe)

properties = { 'unrelaxed_energy_per_atom' : unrelaxed_energies, 'relaxed_energy_per_atom' : relaxed_energies }
