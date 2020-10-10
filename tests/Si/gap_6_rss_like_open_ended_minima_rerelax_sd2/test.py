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
from ase.constraints import UnitCellFilter

# set of utility routines specific this this model/testing framework
from utilities_old import relax_atoms, relax_atoms_cell

# the current model
#print "import model"
import model
#print "done import model"

fmax = 0.01 # maximum force following relaxtion [eV/A]

if os.path.basename(os.path.dirname(model.__file__)) == 'CASTEP_file':
    #print "using tpsd for CASTEP_file"
    model.geom_method='tpsd'
    model.geom_force_tol=fmax
    model.geom_stress_tol=fmax
    model.geom_energy_tol=1.0e8
    model.geom_disp_tol=1.0e8

if not hasattr(model, 'bulk_reference'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=5.43)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_atoms_cell(bulk, tol=fmax, traj_file=None)
else:
    bulk = model.bulk_reference.copy()
    bulk.set_calculator(model.calculator)

#print "calling bulk.get_potential_energy()"
e0_per_atom = bulk.get_potential_energy()/len(bulk)
#print "got e0_per_atom ", e0_per_atom

def struct_energy(e0_per_atom, filename):
    V_unrel = []
    V_rel = []
    E_unrel = []
    E_rel = []

    ats = ase.io.read(os.path.dirname(__file__)+"/"+filename, index=':', format='extxyz')
    for (i, at_orig) in enumerate(ats):
        at = at_orig.copy()
        at.set_calculator(model.calculator)
        # # The commented out bit that follows seems to break CASTEP_file - no minimization runs are set up, apparently because call to relax_atoms_cell doesn't actually trigger a calculate()
        # unrelaxed_e_per_atom  = at.get_potential_energy()/len(at)
        # print 'at unrelaxed energy per atom', unrelaxed_e_per_atom
        # E_unrel.append(unrelaxed_e_per_atom-e0_per_atom)

        # relax fully
        #print "BOB: calling relax_atoms_cell"
        if os.path.basename(os.path.dirname(model.__file__)) == 'CASTEP_file':
            at = relax_atoms_cell(at, tol=fmax, traj_file="model-"+model.name+"-RSS-{}-".format(i)+filename+"-relaxed.opt.xyz")
        else:
            n_at = len(at)
            config_minim = UnitCellFilter(at)
            x = config_minim.get_positions()
            converged = False
            for i_minim in range(500):
                config_minim.set_positions(x)
                grad_f = - config_minim.get_forces()
                E = config_minim.get_potential_energy()
                try:
                    pred_err = predictive_error(config_minim.atoms)
                except:
                    pred_err = None
                if isinstance(config_minim, UnitCellFilter):
                    print(i,"SD2: {} {} {} {}".format(i_minim, E, np.max(np.abs(grad_f[:n_at])), np.max(np.abs(grad_f[n_at:]))))
                else:
                    print(i,"SD2: {} {} {}".format(i_minim, E, np.max(np.abs(grad_f[:n_at]))))
                if np.max(np.abs(grad_f[:n_at])) < fmax and np.max(np.abs(grad_f[n_at:])) < fmax:
                    converged = True
                    break
                if i_minim == 0:
                    alpha = 1.0e-6
                else:
                    alpha = np.sum((x-x_old)*(grad_f-grad_f_old)) / np.sum((grad_f-grad_f_old)**2)
                x_old = x.copy()
                grad_f_old = grad_f.copy()
                x -= np.abs(alpha)*grad_f
        #print "BOB: done relax_atoms_cell"
        ase.io.write(sys.stdout, at, format='extxyz')

        #print "BOB: calling get_potential_energy for relaxed"
        relaxed_e_per_atom  = at.get_potential_energy()/len(at)
        #print "BOB: done get_potential_energy for relaxed"
        #print 'at relaxed energy per atom', relaxed_e_per_atom
        E_rel.append(relaxed_e_per_atom-e0_per_atom)

        V_rel.append(at.get_volume()/len(at))

    # This is a workaround for the problem mentioned above
    ats = ase.io.read(os.path.dirname(__file__)+"/"+filename, index=':', format='extxyz')
    for (i, at) in enumerate(ats):
        at.set_calculator(model.calculator)

        # relax fully
        #print "BOB: calling get_potential_energy for unrelaxed"
        unrelaxed_e_per_atom  = at.get_potential_energy()/len(at)
        #print "BOB: done get_potential_energy for unrelaxed"
        #print 'at unrelaxed energy per atom', unrelaxed_e_per_atom
        E_unrel.append(unrelaxed_e_per_atom-e0_per_atom)
        V_unrel.append(at.get_volume()/len(at))
    # end of workaround

    return (V_unrel, E_unrel, V_rel, E_rel)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
(V_unrel, E_unrel, V_rel, E_rel) = struct_energy(e0_per_atom, "model-GAP-6-test-RSS_like_open_ended-unique-minima.xyz")
properties = {'V_unrelaxed' : V_unrel, 'E_unrelaxed' : E_unrel, 'V_relaxed' : V_rel, 'E_relaxed' : E_rel }
