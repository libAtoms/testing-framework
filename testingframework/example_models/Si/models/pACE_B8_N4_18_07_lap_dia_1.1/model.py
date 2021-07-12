import os
from ase.calculators.lammpsrun import LAMMPS

os.environ["ASE_LAMMPSRUN_COMMAND"] = "/Users/Cas/gits/lammps-ace/src/lmp_serial"

model_dir = os.path.dirname(os.path.realpath(__file__))

parameters = {
    "pair_style": "pace",
    "pair_coeff": ["* * Si_B8_N4_18_07_lap_dia_1.1_rep_2B+ACE.ace Si"],
}

files = [os.path.join(model_dir, "Si_B8_N4_18_07_lap_dia_1.1_rep_2B+ACE.ace")]

calculator = LAMMPS(parameters=parameters, files=files)

name = "ACE"

no_checkpoint = True
