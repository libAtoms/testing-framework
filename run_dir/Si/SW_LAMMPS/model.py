from ase.calculators.lammpsrun import LAMMPS
import os

os.environ["ASE_LAMMPSRUN_COMMAND"]="/Users/Cas/miniconda3/bin/lmp_serial"

parameters = {'pair_style': 'sw',
              'pair_coeff': ['* * Si.sw Si']}

model_dir = os.path.dirname(os.path.realpath(__file__))
os.path.join(model_dir, "Si.sw")

files = [os.path.join(model_dir, "Si.sw")]

calculator = LAMMPS(parameters=parameters, files=files)
