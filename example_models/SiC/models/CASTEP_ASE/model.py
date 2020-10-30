import os

from ase.calculators.castep import Castep

model_abs_dir = os.path.abspath(os.path.dirname(__file__))

mpirun = "mpirun"
mpirun_args = "-n 16"
castep = "castep.mpi"

os.environ['CASTEP_COMMAND'] = '{0} {1} {2}'.format(mpirun, mpirun_args, castep)

no_checkpoint = True


def start(test_name):
    global calculator
    calculator = Castep(directory=test_name,
                        cut_off_energy=500,
                        kpoint_mp_spacing=0.05,
                        xc_functional='pbesol',
                        elec_energy_tol=1E-6,
                        finite_basis_corr='none',
                        spin_polarized=False,
                        perc_extra_bands=100,
                        max_scf_cycles=200,
                        num_dump_cycles=0,
                        mixing_scheme='pulay',
                        mix_history_length=20,
                        smearing_width='0.2',
                        fix_occupancy=False,
                        calculate_stress=True,
                        species_pot=[
                            ("C", "{}/C_C19_PBESOL_OTF.usp".format(model_abs_dir)),
                            ("Si", "{}/Si_C19_PBESOL_OTF.usp".format(model_abs_dir)),
                        ],
                        )
