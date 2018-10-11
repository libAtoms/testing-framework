# Model for VASP for Si

from quippy import Potential, Dictionary
import os.path

# define calculator symbol

mydir=os.path.abspath(os.path.dirname(__file__))
calculator=Potential("FilePot command=vasp_driver")
calculator.set_calc_args(Dictionary('vasp=vasp.para INCAR_template="{0}/INCAR.template" kpoints_file=_NONE_ potcar_files={{14 {0}/POTCAR.Si}}'.format(mydir)))
calculator.set_default_properties(['energy','forces','stress'])

# want checkpointing, but broken now because of force_consistent
# no_checkpoint = False
no_checkpoint = True
name = 'VASP'

