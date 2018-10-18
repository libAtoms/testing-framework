import ase.io
from itertools import izip
from utilities import robust_minim_cell_pos

def do_RSS(initial_configs_file, index=':', tol=0.01):
    import model

    ats = ase.io.read(initial_configs_file, index)
    range_slice_args = [ None if i == '' else int(i) for i in index.split(':')]

    print "got index ", index, "range_slice_args",range_slice_args
    print "using i_config",range(len(ats))[slice(*range_slice_args)]

    energies = []
    volumes = []
    for (i_config, at) in izip(range(len(ats))[slice(*range_slice_args)], ats):
        robust_minim_cell_pos(at, tol, "RSS_%04d" % i_config)
        print "RSS completed minimization"
        if hasattr(model, "fix_cell_dependence"):
            model.fix_cell_dependence()
        energies.append(at.get_potential_energy()/len(at))
        volumes.append(at.get_volume()/len(at))
        ase.io.write("RSS_relaxed_%04d.extxyz" % i_config, at)

    return { 'energies' : energies, 'volumes' : volumes }
