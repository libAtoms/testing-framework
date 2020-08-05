import os
from distutils import spawn

#from matscipy.socketcalc import CastepClient, SocketCalculator
from ase.calculators.castep import Castep

model_abs_dir = os.path.abspath(os.path.dirname(__file__))

mpirun = "mpirun"
mpirun_args = "-n 16"
castep = "castep.mpi"

os.environ['CASTEP_COMMAND'] = '{0} {1} {2}'.format(mpirun, mpirun_args, castep)

print(os.environ['CASTEP_COMMAND'])

no_checkpoint = True

#def start(test_name):
#       global calculator
#       calculator = Castep()
#
#       calculator._directory=test_name
#
#       calculator.param.cut_off_energy=250
#       calculator.param.spin_polarised=False
#       calculator.param.mixing_scheme='Pulay'
#       calculator.param.fine_grid_scale=4
#       calculator.param.smearing_width=0.05
#       calculator.param.finite_basis_corr='automatic'
#       calculator.param.calculate_stress=True
#
#       calculator.cell.kpoints_mp_spacing=0.03

def start(test_name):
    global calculator
    calculator = Castep(directory=test_name,
                        cut_off_energy=250,
                        spin_polarized=False,
                        opt_strategy='speed',
                        xc_functional='PW91',
                        elec_energy_tol='0.0000001',
                        max_scf_cycles=250,
                        fix_occupancy=False,
                        calculate_stress=True,
                        finite_basis_corr='automatic',
                        smearing_width='0.05',
                        fine_grid_scale=4,
                        mixing_scheme='pulay',
                        mix_history_length=20,
                        num_dump_cycles=0,
                        kpoints_mp_spacing='0.030', # note that other values were used for some tests, e.g. 0.015 for bulk E(V), 0.07 for minimization of GAP amorphous structures, and perhaps some other variations
                        species_pot = ("Si","{}/Si_OTF.usp".format(model_abs_dir)),
                        perc_extra_bands=200)

#                        symmetry_generate=False,
