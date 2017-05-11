import os

from quippy import Potential
import __builtin__

orig_dir = os.getcwd()
model_dir = os.path.dirname(__file__)
if model_dir != '':
    os.chdir(model_dir)
try:
    if hasattr(__builtin__, 'mpi_glob'):
        calculator = Potential(init_args='Potential xml_label="GAP_2017_5_7_60_20_48_12_148"',
                                               param_filename='gp_merged_2bdmbd.xml', mpi_obj=mpi_glob)
    else:
        calculator = Potential(init_args='Potential xml_label="GAP_2017_5_7_60_20_48_12_148"',
                                               param_filename='gp_merged_2bdmbd.xml')
    Potential.__str__ = lambda self: '<GAP Potential>'
finally:
    os.chdir(orig_dir)

no_checkpoint = True

name = 'GAP_2017_5_7_60_20_48_12_148'
