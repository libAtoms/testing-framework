from ase.calculators.lammpsrun import LAMMPS
import os

os.environ["ASE_LAMMPSRUN_COMMAND"]="/Users/Cas/miniconda3/bin/lmp_serial"

parameters = {'pair_style': 'sw',
              'pair_coeff': ['* * Si.sw Si']}

files = ['/Users/Cas/gits/testing-framework/run_dir/Si/models/SW_LAMMPS/Si.sw']

calculator = LAMMPS(parameters=parameters, files=files)
