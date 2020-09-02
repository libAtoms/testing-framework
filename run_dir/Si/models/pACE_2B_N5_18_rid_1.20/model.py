import os
from ase.calculators.lammpsrun import LAMMPS

os.environ["ASE_LAMMPSRUN_COMMAND"]="/Users/Cas/gits/lammps-ace/src/lmp_serial"

model_dir = os.path.dirname(os.path.realpath(__file__))

parameters = {'pair_style': 'pace',
              'pair_coeff': ['* * Si_2B_N5_18_072_rid_1.2_rep_2B+ACE.ace Si']}

files = [os.path.join(model_dir, "Si_2B_N5_18_072_rid_1.2_rep_2B+ACE.ace")]

calculator = LAMMPS(parameters=parameters, files=files)

no_checkpoint = True
