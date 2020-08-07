import os
from ase.calculators.lammpsrun import LAMMPS

os.environ["ASE_LAMMPSRUN_COMMAND"]="lmp_serial"

parameters = {'pair_style': 'sw',
              'pair_coeff': ['* * Si.sw Si']}

model_dir = os.path.dirname(os.path.realpath(__file__))

files = [os.path.join(model_dir, "Si.sw")]

calculator = LAMMPS(parameters=parameters, files=files)

no_checkpoint = True
