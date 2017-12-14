import os.path, ase.io
from utilities import relax_config, run_root
import random, sys, model

test_dir = os.path.abspath(os.path.dirname(__file__))

bulk = ase.io.read(os.path.join(test_dir,"bulk.xyz"), format="extxyz")
tol=1.0e-3
relaxed_bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, 
    traj_file=None, config_label='bulk', from_base_model=True, save_config=True, keep_symmetry=True)
relaxed_bulk_pe = relaxed_bulk.get_potential_energy()/len(relaxed_bulk)
ase.io.write(sys.stdout, relaxed_bulk, format="extxyz")

tol=2.0e-3 # hack to work around convergence issues
antisites = relaxed_bulk * (2,2,2)

random.seed(10)
unrelaxed_energies = []
relaxed_energies = []
for config_i in range(20):
    Z = antisites.get_atomic_numbers()
    (i1, i2) = random.sample(range(len(antisites)),2)
    while Z[i1] == Z[i2]:
        (i1, i2) = random.sample(range(len(antisites)),2)
    t_Z = Z[i1]
    Z[i1] = Z[i2]
    Z[i2] = t_Z
    antisites.set_atomic_numbers(Z)
    ase.io.write(sys.stdout, antisites, format="extxyz")
    antisites.set_calculator(model.calculator)
    unrelaxed_energies.append(antisites.get_potential_energy()/len(antisites)-relaxed_bulk_pe)
    relaxed_antisites = antisites.copy()
    relaxed_antisites = relax_config(relaxed_antisites, relax_pos=True, relax_cell=True, tol=tol, 
        traj_file=None, config_label='supercell_antisite_{}'.format(config_i), from_base_model=True, save_config=True)
    relaxed_energies.append(relaxed_antisites.get_potential_energy()/len(antisites)-relaxed_bulk_pe)
    ase.io.write(sys.stdout, relaxed_antisites, format="extxyz")

properties = { 'unrelaxed_energy_per_atom' : unrelaxed_energies,
               'relaxed_energy_per_atom' : relaxed_energies }
