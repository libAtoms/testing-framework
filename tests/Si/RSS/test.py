import ase.io, numpy as np, os.path
from utilities import relax_config

np.random.seed(17)

N = 3
tol=1.0e-3

ats = ase.io.read(os.path.join(os.path.abspath(os.path.dirname(__file__)),'initial_configs.xyz'),':')
ats = np.random.choice(ats, N, replace=False)

energies = []
volumes = []
for (i_at, at) in enumerate(ats):
    scaled_bulk = relax_config(at, relax_pos=True, relax_cell=True, tol=tol, traj_file="RSS_traj_%04d.extxyz" % i_at, method='lbfgs', keep_symmetry=True, config_label="RSS_%04d" % i_at)
    energies.append(at.get_potential_energy()/len(at))
    volumes.append(at.get_volume()/len(at))

properties = { 'energies' : energies, 'volumes' : volumes }
