import ase.io, numpy as np, os.path
from utilities import robust_minim
import model
import minim_sd2

np.random.seed(17)

N = 3
tol=2.0e-2

ats = ase.io.read(os.path.join(os.path.abspath(os.path.dirname(__file__)),'initial_configs.xyz'),':')
ats = np.random.choice(ats, N, replace=False)

energies = []
volumes = []
for (i_at, at) in enumerate(ats):
    ## robust_minim(at, tol, "RSS_%04d" % i_at)
    robust_minim(at, tol, label="RSS_%04d" % i_at, max_sd2_iter=50, max_lbfgs_iter=20)
    if hasattr(model, "new_config"):
        model.new_config(at)
    energies.append(at.get_potential_energy()/len(at))
    volumes.append(at.get_volume()/len(at))
    ase.io.write("RSS_relaxed_%04d.extxyz" % i_at, at)

properties = { 'energies' : energies, 'volumes' : volumes }
