from ase.io import read, write
from ase.calculators.calculator import Calculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import FIRE
from ase.optimize.precon import PreconLBFGS
from ase.constraints import ExpCellFilter, FixAtoms, voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_stress
from ase.units import GPa
import numpy as np
import os.path
import sys
import time
from ase import Atoms

try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
except:
    pass

from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry, refine_symmetry

#from quippy.io import AtomsWriter
#from quippy.cinoutput import CInOutput,OUTPUT
#from quippy.atoms import Atoms

def path_of_file(file):
    return os.path.abspath(os.path.dirname(file))

def name_of_file(file):
    if '/' in file:
        name = os.path.basename(file)
    else:
        name = file
    return name

def model_test_root(u_model_name=None, u_test_name=None, base_model=False):
    if u_model_name is None:
        if base_model:
            u_model_name = base_model_name
        else:
            u_model_name = model_name
    if u_test_name is None:
        u_test_name = test_name
    if system_label != '':
        return '{0}-model-{1}-test-{2}'.format(system_label, u_model_name, u_test_name)
    else:
        return 'model-{0}-test-{1}'.format(model_name, test_name)

def sd2_run(log_prefix, config_minim, tol, converged, max_iter=1000, max_cell_vec_change_ratio=3.0, exception_if_invalid_config = None):
    traj = []
    x = config_minim.get_positions()
    # x_old = 0.0
    # grad_f_old = 0.0
    try:
        initial_cell_mag = np.linalg.norm(config_minim.atoms.get_cell(), axis=1)
    except:
        initial_cell_mag = None

    try:
        underlying_atoms = config_minim.atoms
    except:
        underlying_atoms = config_minim

    return_stat = "unconverged"
    for i_minim in range(max_iter):
        if "n_minim_iter" in underlying_atoms.info:
            underlying_atoms.info["n_minim_iter"] += 1
        if initial_cell_mag is not None:
            cur_cell_mag = np.linalg.norm(config_minim.atoms.get_cell(), axis=1)
            if np.max(cur_cell_mag/initial_cell_mag) > max_cell_vec_change_ratio:
                print("SD2: i {} bad lattice constant".format(i_minim))
                sys.stdout.flush()
                return_stat="failed"
                break

        config_minim.set_positions(x)
        if exception_if_invalid_config is not None:
            exception_if_invalid_config(config_minim)
        grad_f = - config_minim.get_forces()
        E = config_minim.get_potential_energy()
        done = converged(i_minim)

        try: # ExpCellFilter
            traj.append(config_minim.atoms.copy())
        except:
            traj.append(config_minim.copy())

        if done:
            print(log_prefix,"SD2: i {} E {} {} {}".format(i_minim, E, config_minim.log_message, done))
            sys.stdout.flush()
            return_stat="converged"
            break

        if i_minim == 0:
            alpha = 1.0e-4
        else:
            alpha_num = np.sum((x-x_old)*(grad_f-grad_f_old))
            if alpha_num < 0.0:
                alpha = 1.0e-2
            else:
                alpha = alpha_num / np.sum((grad_f-grad_f_old)**2)

        print(log_prefix,"SD2: i {} E {} {} alpha {} {}".format(i_minim, E, config_minim.log_message, alpha, done))
        sys.stdout.flush()

        x_old = x.copy()
        grad_f_old = grad_f.copy()
        x -= alpha*grad_f

    return (traj, return_stat)

def f_conv_crit_sq(forces):
    return (forces**2).sum(axis=1).max()
def s_conv_crit_sq(stress):
    return (stress**2).max()

def sd2_converged(minim_ind, atoms, fmax, smax=None):
    forces = atoms.get_forces()

    if isinstance(atoms, ExpCellFilter):
        fmax_sq = f_conv_crit_sq(forces[:len(atoms)-3])
        smax_sq = s_conv_crit_sq(forces[len(atoms)-3:])
        if smax is None:
            smax = fmax
        atoms.log_message = "f {} s {} ".format(np.sqrt(fmax_sq), np.sqrt(smax_sq))
        f_conv = (fmax_sq < fmax**2 and smax_sq < smax**2)
    else:
        fmax_sq = f_conv_crit_sq(forces)
        atoms.log_message = "f {} ".format(np.sqrt(fmax_sq))
        f_conv = fmax_sq < fmax**2

    return f_conv


def relax_config(atoms, relax_pos, relax_cell, tol=1e-3, method='lbfgs', max_steps=200, save_traj=False, constant_volume=False,
    refine_symmetry_tol=None, keep_symmetry=False, strain_mask = None, config_label=None, from_base_model=False, save_config=False,
    try_restart=False, fix_cell_dependence=False, applied_P=0.0, **kwargs):

    # get from base model if requested
    import model

    if config_label is not None:
        save_file = run_root+'-'+config_label+'-relaxed.xyz'
        traj_file = run_root+'-'+config_label+'-traj.xyz'
    else:
        save_file = None
        traj_file = None

    if try_restart:
        # try from saved final config
        try:
            atoms = read(save_file, format='extxyz')
            print("relax_config read config from final file", save_file)
            return atoms
        except:
            pass
        # try from last config in traj file
        try:
            atoms_in = read(traj_file, -1, format='extxyz')
            # set positions from atoms_in rescaling to match current cell
            saved_cell = atoms.get_cell().copy()
            atoms.set_cell(atoms_in.get_cell())
            atoms.set_positions(atoms_in.get_positions())
            atoms.set_cell(saved_cell, scale_atoms=True)
            print("relax_config read config from traj", traj_file)
        except:
            pass
    else:
        if from_base_model:
            if config_label is None:
                raise ValueError('from_base_model is set but no config_label provided')
            try:
                base_run_file = os.path.join('..',base_run_root,base_run_root+'-'+config_label+'-relaxed.xyz')
                atoms_in = read(base_run_file, format='extxyz')
                # set positions from atoms_in rescaling to match current cell
                saved_cell = atoms.get_cell().copy()
                atoms.set_cell(atoms_in.get_cell())
                atoms.set_positions(atoms_in.get_positions())
                atoms.set_cell(saved_cell, scale_atoms=True)
                print("relax_config read config from base model ",base_run_file)
            except:
                try:
                    print("relax_config failed to read base run config from ",base_run_root+'-'+config_label+'-relaxed.xyz')
                except:
                    print("relax_config failed to determined base_run_root")

    print("relax_config symmetry before refinement at default tol 1.0e-6")
    check_symmetry(atoms, 1.0e-6, verbose=True)
    if refine_symmetry_tol is not None:
        refine_symmetry(atoms, refine_symmetry_tol)
        print("relax_config symmetry after refinement")
        check_symmetry(atoms, refine_symmetry_tol, verbose=True)
    if keep_symmetry:
        print("relax_config trying to maintain symmetry")
        atoms.set_constraint(FixSymmetry(atoms))

    atoms.set_calculator(model.calculator)

    # if needed, fix cell dependence before running
    if fix_cell_dependence and hasattr(model, "fix_cell_dependence"):
        model.fix_cell_dependence(atoms)

    if method == 'lbfgs' or method == 'sd2':
        if 'move_mask' in atoms.arrays:
            atoms.set_constraint(FixAtoms(np.where(atoms.arrays['move_mask'] == 0)[0]))
        if relax_cell:
            atoms_cell = ExpCellFilter(atoms, mask=strain_mask, constant_volume=constant_volume, scalar_pressure=applied_P*GPa)
        else:
            atoms_cell = atoms
        atoms.info["n_minim_iter"] = 0
        if method == 'sd2':
            (traj, run_stat) = sd2_run("", atoms_cell, tol, lambda i : sd2_converged(i, atoms_cell, tol), max_steps)
            if save_traj is not None:
                write(traj_file, traj)
        else:
            # precon="Exp" specified to resolve an error with the lbfgs not optimising
            opt = PreconLBFGS(atoms_cell, use_armijo=False, **kwargs)
            if save_traj:
                traj = open(traj_file, "w")
                def write_trajectory():
                    if "n_minim_iter" in atoms.info:
                        atoms.info["n_minim_iter"] += 1
                    write(traj, atoms, format='extxyz')
                    traj.flush()
                opt.attach(write_trajectory)
    elif method == 'cg_n':
        raise ValueError('minim method cg_n not supported in new python3 quippy')
        # if strain_mask is not None:
            # raise(Exception("strain_mask not supported with method='cg_n'"))
        # atoms.info['Minim_Constant_Volume'] = constant_volume
        # opt = Minim(atoms, relax_positions=relax_pos, relax_cell=relax_cell, method='cg_n')
    else:
        raise ValueError('unknown method %s!' % method)

    if method != 'sd2':
        opt.run(tol, max_steps)

    if refine_symmetry_tol is not None:
        print("symmetry at end of relaxation at desired tol")
        check_symmetry(atoms, refine_symmetry_tol, verbose=True)
    print("symmetry at end of relaxation at default tol 1e-6")
    check_symmetry(atoms, 1.0e-6, verbose=True)

    # in case we had a trajectory saved
    try:
        traj.close()
    except:
        pass

    if save_config:
        if config_label is None:
            raise ValueError('save_config is set but no config_label provided')
        write(save_file, atoms, format='extxyz')

    if keep_symmetry:
        for (i_c, c) in enumerate(atoms.constraints):
            if isinstance(c, FixSymmetry):
                del atoms.constraints[i_c]
                break

    # undo fix cell dependence
    if fix_cell_dependence and hasattr(model, "fix_cell_dependence"):
        sys.stderr.write("WARNING: relax_config undoing fix_cell_dependence, whether or not it was set before it started\n")
        model.fix_cell_dependence()

    return atoms

def evaluate_atoms_list(atoms_list, do_energy=True, do_forces=True, do_stress=True):
    results = []
    for (at_i, at) in enumerate(atoms_list):
        print("evaluation ",at_i,"/",len(atoms_list))
        evaluate(at, do_energy, do_forces, do_stress)
        result = {}
        if do_energy:
            result["energy"] = at.get_potential_energy()
        if do_forces:
            result["forces"] = at.get_forces().tolist()
        if do_stress:
            result["stress"] = at.get_stress().tolist()
        results.append(result)
    return results

def evaluate_file(file, do_energy=True, do_forces=True, do_stress=True, do_predictive_error=None):
    al = read(file, index=":")
    for a in al:
        evaluate(a, do_energy, do_forces, do_stress, do_predictive_error)
    write(model_test_root()+"_"+name_of_file(file), al)
    return al

def evaluate(atoms, do_energy=True, do_forces=True, do_stress=True, do_predictive_error=None):
    import model

    if do_predictive_error:
        try:
            orig_calc_args = model.calculator.get_calc_args()
            new_calc_args = orig_calc_args.copy()
            new_calc_args["local_gap_variance"] = "predictive_error"
            if isinstance(do_predictive_error,float):
                new_calc_args["gap_variance_regularisation"] = do_predictive_error
            elif not isinstance(do_predictive_error,bool):
                raise ValueError("do_predictive_error not float or bool {}".format(str(do_predictive_error)))
            model.calculator.set_calc_args(new_calc_args)
        except AttributeError:
            pass

    atoms.set_calculator(model.calculator)

    stress = None
    if do_stress:
        stress = atoms.get_stress()

    forces = None
    if do_forces:
        forces = atoms.get_forces()

    energy = None
    if do_energy:
        energy = atoms.get_potential_energy()

    spc = SinglePointCalculator(atoms,
                                energy=energy,
                                forces=forces,
                                stress=stress)
    atoms.set_calculator(spc)

    if do_predictive_error:
        atoms.arrays["predictive_error"] = np.sqrt(model.calculator.results["predictive_error"])
        try:
            model.calculator.set_calc_args(orig_calc_args)
        except AttributeError:
            pass

    return atoms

def robust_minim_cell_pos(atoms, final_tol, label="robust_minim", max_sd2_iter=50, sd2_tol=1.0, max_lbfgs_iter=20, max_n_lbfgs=50, keep_symmetry=True):
    import model

    # do each minim at fixed cell-dependent model params (e.g. k-point mesh)
    if hasattr(model, "fix_cell_dependence"):
        model.fix_cell_dependence(atoms)
    relax_config(atoms, relax_pos=True, relax_cell=True, tol=sd2_tol, max_steps=max_sd2_iter,
        save_traj=True, method='sd2', keep_symmetry=keep_symmetry, config_label=label )

    done=False
    i_iter = 0
    while not done and i_iter < max_n_lbfgs:
        try:
            if hasattr(model, "fix_cell_dependence"):
                model.fix_cell_dependence(atoms)
            relax_config(atoms, relax_pos=True, relax_cell=True, tol=final_tol, max_steps=max_lbfgs_iter,
                save_traj=True, method='lbfgs', keep_symmetry=keep_symmetry, config_label=f'{label}.{i_iter}' )
            done = (atoms.info["n_minim_iter"] < max_lbfgs_iter)
            print("robust_minim relax_configs LBFGS finished in ",atoms.info["n_minim_iter"],"iters, max", max_lbfgs_iter)
        except:
            print("robust_minim relax_configs LBFGS failed, trying again")
        i_iter += 1

    # Undo fixed cell dependence. Hopefully no one is using robust_minim as part of a
    # more complex process that is doing its own fix_cell_depdence()
    if hasattr(model, "fix_cell_dependence"):
        model.fix_cell_dependence()

def get_relaxed_bulk(bulk_struct_test, model_name=None):
    bulk_model_test_relaxed = os.path.join('..',model_test_root(u_model_name=model_name, u_test_name=bulk_struct_test)+"-relaxed.xyz")
    try:
        bulk = read(bulk_model_test_relaxed, format='extxyz')
    except:
        sys.stderr.write("Failed to read relaxed bulk '{}', perhaps bulk test hasn't been run yet\n".format(bulk_model_test_relaxed))
        sys.exit(1)

    return bulk

def rescale_to_relaxed_bulk(supercell):
    # read bulk
    bulk = get_relaxed_bulk(supercell.info['bulk_struct_test'])

    # rescale supercell cell
    try:
        supercell_a1_lattice = supercell.info['supercell_a1_in_bulk_lattice_coords']
        bulk_lattice = bulk.get_cell()
        supercell_a1_in_bulk = np.dot(supercell_a1_lattice,bulk_lattice)
        cell_ratio = np.linalg.norm(supercell_a1_in_bulk) / np.linalg.norm(supercell.get_cell()[0,:])
        if 'supercell_a2_in_lattice_coords' in supercell.info:
            raise ValueError('anisotropic rescaling of supercellace cell not implemented')
        if 'supercell_a3_in_lattice_coords' in supercell.info:
            raise ValueError('anisotropic rescaling of supercellace cell not implemented')
    except:
        print("'supercell_a1_in_bulk_lattice_coords' is not in supercell.info (imported from surface.xyz). Assuming a cell_ratio of 1.0")
        cell_ratio = 1.0

    supercell.set_cell(supercell.get_cell()*cell_ratio, scale_atoms=True)

    return bulk

def phonons(model,bulk,supercell,dx,mesh=None,points=None,n_points=50):

    import model

    unitcell = PhonopyAtoms(symbols=bulk.get_chemical_symbols(),
                            cell=bulk.get_cell(),
                            scaled_positions=bulk.get_scaled_positions())
    phonon = Phonopy(unitcell,supercell)
    phonon.generate_displacements(distance=dx)

    sets_of_forces = []

    for s in phonon.get_supercells_with_displacements():
        at = Atoms(cell=s.get_cell(),
                       symbols=s.get_chemical_symbols(),
                       scaled_positions=s.get_scaled_positions(),
                       pbc=3*[True])
        at.set_calculator(model.calculator)
        sets_of_forces.append(at.get_forces())

    phonon.set_forces(sets_of_forces=sets_of_forces)
    phonon.produce_force_constants()

    properties = {}

    if mesh is not None:
        phonon.set_mesh(mesh,is_gamma_center=True)
        qpoints, weights, frequencies, eigvecs = phonon.get_mesh()

        properties["frequencies"] = frequencies.tolist()
        properties["weights"] = weights.tolist()

    if points is not None:
        bands = []
        for i in range(len(points)-1):
            band = []
            for r in np.linspace(0,1,n_points):
                band.append(points[i]+(points[i+1]-points[i])*r)
            bands.append(band)

        phonon.set_band_structure(bands,is_eigenvectors=True,is_band_connection=False)
        band_q_points, band_distances, band_frequencies, band_eigvecs = phonon.get_band_structure()

        band_distance_max = np.max(band_distances)
        band_distances = [(_b / band_distance_max).tolist() for _b in band_distances]

        band_frequencies = [ _b.tolist() for _b in band_frequencies ]

        properties["band_q_points"] = band_q_points
        properties["band_distances"] = band_distances
        properties["band_frequencies"] = band_frequencies
        properties["band_eigvecs"] = band_eigvecs

    properties["phonopy"] = phonon

    return properties
