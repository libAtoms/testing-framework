import os
import quippy
from ase.calculators.loggingcalc import LoggingCalculator


model_dir = os.path.dirname(__file__)

if 'USE_NP' in os.environ:
    np = int(os.environ['USE_NP'])
else:
    np = 1

if 'VASP_COMMAND' in os.environ:
      vasp = os.environ['VASP_COMMAND']
else:
      if np > 1:
         vasp_cmd = ('vasp.para')
      else:
         vasp_cmd = ('vasp.serial')

name = 'VASP'

no_checkpoint = True

global calculator

if not os.path.isdir(model_dir+"/pot"):
    print "ERROR: Can't find pot/ subdirectory undel model_dir %s" % model_dir
    raise

VASP = quippy.Potential('FilePot command=vasp_driver')
VASP.set(vasp=vasp_cmd, incar_template='%s/INCAR.template' % model_dir, kpoints_file='_NONE_', potcar_files='32 %s/pot/Ge/POTCAR 51 %s/pot/Sb/POTCAR 52 %s/pot/Te/POTCAR' % (model_dir, model_dir, model_dir), persistent=False, clean_up_keep_n=10)
VASP.set_default_properties(['forces', 'energy', 'stress'])
calculator = LoggingCalculator(VASP)
