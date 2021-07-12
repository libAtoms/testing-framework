# from ase.calculators.eam import EAM
# import os
#
# model_dir = os.path.dirname(os.path.realpath(__file__))
#
# file = os.path.join(model_dir, "V_Olsson_CMS2009.eam.alloy")
#
# calculator = EAM(potential=file)
#
# no_checkpoint = True
#
# name = "ACE"


import os
from ase.calculators.lammpsrun import LAMMPS

os.environ["ASE_LAMMPSRUN_COMMAND"] = "lmp_serial"

model_dir = os.path.dirname(os.path.realpath(__file__))

parameters = {
    "pair_style": "eam/alloy",
    "pair_coeff": ["* * V_Olsson_CMS2009.eam.alloy V"],
}

files = [os.path.join(model_dir, "V_Olsson_CMS2009.eam.alloy")]

calculator = LAMMPS(parameters=parameters, files=files)

no_checkpoint = True
