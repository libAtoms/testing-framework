from ase.io import read, write
from ase.calculators.calculator import Calculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import FIRE
from ase.optimize.precon import Exp, PreconLBFGS
from ase.constraints import UnitCellFilter, FixAtoms, voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_stress
from quippy.potential import Minim
import numpy as np
import os.path

import symmetrize

#from quippy.io import AtomsWriter
#from quippy.cinoutput import CInOutput,OUTPUT
#from quippy.atoms import Atoms

def name_of_file(file):
    if '/' in file:
        name = os.path.basename(os.path.dirname(file))
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

class SymmetrizedCalculator(Calculator):
   implemented_properties = ['energy','forces','stress']
   def __init__(self, calc, atoms, *args, **kwargs):
      Calculator.__init__(self, *args, **kwargs)
      self.calc = calc
      (self.rotations, self.translations, self.symm_map) = symmetrize.prep(atoms)

   def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if system_changes:
            self.results = {}
        if 'energy' in properties and 'energy' not in self.results:
            self.results['energy'] = self.calc.get_potential_energy(atoms)
        if 'forces' in properties and 'forces' not in self.results:
            raw_forces = self.calc.get_forces(atoms)
            self.results['forces'] = symmetrize.forces(atoms.get_cell(), atoms.get_reciprocal_cell().T, raw_forces, 
                self.rotations, self.translations, self.symm_map)
        if 'stress' in properties and 'stress' not in self.results:
            raw_stress = voigt_6_to_full_3x3_stress(self.calc.get_stress(atoms))
            symmetrized_stress = symmetrize.stress(atoms.get_cell(), atoms.get_reciprocal_cell().T, raw_stress, self.rotations)
            self.results['stress'] = full_3x3_to_voigt_6_stress(symmetrized_stress)

def relax_config(atoms, relax_pos, relax_cell, tol=1e-3, method='lbfgs', max_steps=200, traj_file=None, constant_volume=False,
    refine_symmetry_tol=None, keep_symmetry=False, strain_mask = None, config_label=None, from_base_model=False, save_config=False, **kwargs):

    # get from base model if requested
    import model
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
            print "relax_config read config from ",base_run_file
        except:
            try:
                print "relax_config failed to read base run config from ",base_run_root+'-'+config_label+'-relaxed.xyz'
            except:
                print "relax_config failed to determined base_run_root"

    if refine_symmetry_tol is not None:
        symmetrize.refine(atoms, refine_symmetry_tol)
    if keep_symmetry:
        atoms.set_calculator(SymmetrizedCalculator(model.calculator, atoms))
    else:
        atoms.set_calculator(model.calculator)

    if method == 'lbfgs':
        if 'move_mask' in atoms.arrays:
            atoms.set_constraint(FixAtoms(np.where(atoms.arrays['move_mask'] == 0)[0]))
        if relax_cell:
            atoms = UnitCellFilter(atoms, mask=strain_mask, constant_volume=constant_volume)
        precon = Exp(3.0, use_pyamg=False)
        opt = PreconLBFGS(atoms, precon=precon, **kwargs)
        if traj_file is not None:
            traj = open(traj_file, "w")
            def write_trajectory():
                write(traj, atoms, format='extxyz')
            opt.attach(write_trajectory)
    elif method == 'cg_n':
        if strain_mask is not None:
            raise(Exception("strain_mask not supported with method='cg_n'"))
        atoms.info['Minim_Constant_Volume'] = constant_volume
        opt = Minim(atoms, relax_positions=relax_pos, relax_cell=relax_cell, method='cg_n')
    else:
        raise ValueError('unknown method %s!' % method)

    opt.run(tol, max_steps)

    # in case we wrapped in a UnitCellFilter
    try:
        atoms = atoms.atoms
    except:
        pass

    symmetrize.check(atoms, refine_symmetry_tol)
    symmetrize.check(atoms, 1.0e-6)

    # in case we had a trajectory saved
    try:
        traj.close()
    except:
        pass

    if save_config:
        if config_label is None:
            raise ValueError('save_config is set but no config_label provided')
        write(run_root+'-'+config_label+'-relaxed.xyz', atoms, format='extxyz')

    if keep_symmetry:
        atoms.set_calculator(model.calculator)

    return atoms


def evaluate(atoms, do_energy=True, do_forces=True, do_stress=True):
    import model
    atoms.set_calculator(model.calculator)

    energy = None
    if do_energy:
        energy = atoms.get_potential_energy()

    forces = None
    if do_forces:
        forces = atoms.get_forces()

    stress = None
    if do_stress:
        stress = atoms.get_stress()

    spc = SinglePointCalculator(atoms,
                                energy=energy,
                                forces=forces,
                                stress=stress)
    atoms.set_calculator(spc)
    return atoms
