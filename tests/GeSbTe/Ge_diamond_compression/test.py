import os.path, ase.io
from utilities import relax_config, run_root
import random
import quippy

test_dir = os.path.abspath(os.path.dirname(__file__))

bulk = ase.io.read(os.path.join(test_dir,"bulk.xyz"), format="extxyz")
tol=1.0e-3 # hack to work around bad convergence
relaxed_bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, 
    traj_file=None, config_label='bulk', from_base_model=True, save_config=True)
relaxed_bulk_pe = relaxed_bulk.get_potential_energy()/len(bulk)

orig_cell = relaxed_bulk.get_cell()
energies = []
ns=20
for i in range(0,ns+1):
    l_factor = 1.0-0.5*float(i)/float(ns)
    relaxed_bulk.set_cell(orig_cell*l_factor, True)
    qrb = quippy.Atoms(relaxed_bulk)
    qrb.set_cutoff(3.0)
    qrb.calc_connect()
    try:
        energies.append([ min([nn.distance for nn in qrb.neighbours[1]]), relaxed_bulk.get_potential_energy()/len(relaxed_bulk) ])
    except:
        break

properties={ 'unrelaxed_energies' :  energies }
