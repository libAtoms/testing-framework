# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS:
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines

import model
import json

from ase.build import bulk
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


HBAR = 0.6582119514
BOLTZMANN_K = 8.6173303e-05
GPA = 160.21766208

import ase
from ase import Atoms
from testingframework.share.utilities import relax_config, phonons
import numpy as np

supercell = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
dx = 0.03
mesh = [16, 16, 16]

temperature_range = [1.0e-4, 1500]

points = np.array([[0, 0, 0], [0, 0.5, 0.5], [1, 1, 1], [0.5, 0.5, 0.5]])

n_points = 50

# set up cell and bulk modulus
try:
    with open("../model-{}-test-bulk_diamond-properties.json".format(model.name)) as f:
        j = json.load(f)
    a0 = j["diamond_a0"]
    bulk_at = bulk("Si", "diamond", a0)

    bulk_modulus = j["diamond_bulk_modulus"] / GPA
except:
    bulk_at = bulk("Si", "diamond", 5.43)
    bulk_at.set_calculator(model.calculator)
    bulk_at = relax_config(
        bulk_at, relax_pos=True, relax_cell=True, tol=1.0e-4, traj_file=None
    )
    a0 = bulk_at.get_cell_lengths_and_angles()[0] * np.sqrt(2.0)
    bulk_at = bulk("Si", "diamond", a0)

    from matscipy import elasticity
    from ase.optimize import BFGS

    bulk_at.set_calculator(model.calculator)
    opt = BFGS
    b = elasticity.fit_elastic_constants(bulk_at, symmetry="cubic", optimizer=opt)
    bulk_modulus = elasticity.elastic_moduli(b[0])[3]

try:
    with open(
        "../model-{}-test-phonon_diamond-properties.json".format(model.name)
    ) as f:
        j = json.load(f)

    f0 = j["phonon_diamond_frequencies"]
    bf0 = j["phonon_diamond_band_frequencies"]
    band_distances = j["phonon_diamond_band_distances"]
    weights = j["phonon_diamond_weights"]
except:
    phonon_properties = phonons(
        model,
        bulk_at,
        supercell=supercell,
        dx=dx,
        mesh=mesh,
        points=points,
        n_points=n_points,
    )
    f0 = phonon_properties["frequencies"]
    bf0 = phonon_properties["band_frequencies"]
    band_distances = phonon_properties["band_distances"]
    weights = phonon_properties["weights"]

n_frequencies = len(f0)
n_band_frequencies = len(bf0)

f0 = np.array(f0) * 1.0e-3
bf0 = np.array(bf0) * 1.0e-3
weights = 1.0 * np.array(weights)

c = 0.01
fd = []
band_frequencies_d = []
a = []
for dc in [-c, c]:
    atd = bulk_at.copy()
    atd.set_cell(bulk_at.get_cell() * (1 + dc), scale_atoms=True)
    phonon_properties = phonons(
        model,
        atd,
        supercell=supercell,
        dx=dx,
        mesh=mesh,
        points=points,
        n_points=n_points,
    )
    a.append(atd.get_cell_lengths_and_angles()[0] * np.sqrt(2.0))
    f = phonon_properties["frequencies"]
    if n_frequencies != len(f):
        raise ValueError
    fd.append(np.array(f) * 1.0e-3)

    bf = phonon_properties["band_frequencies"]
    if n_band_frequencies != len(bf):
        raise ValueError
    band_frequencies_d.append(np.array(bf) * 1.0e-3)

mode_gruneisens = -(fd[1] - fd[0]) / (a[1] - a[0]) * a0 / f0 / 3.0
where_Nans = np.isclose(f0, 0.0, atol=1.0e-4)
mode_gruneisens[where_Nans] = 0.0

mode_band_gruneisens = (
    -(band_frequencies_d[1] - band_frequencies_d[0]) / (a[1] - a[0]) * a0 / bf0 / 3.0
)
where_Nans = np.isclose(bf0, 0.0, atol=1.0e-4)
# mode_band_gruneisens[where_Nans] = 0.0
for i_direction, direction in enumerate(bf0):
    for i_q, q, in enumerate(direction):
        for i_band, band in enumerate(q):
            if np.isclose(band, 0.0, atol=1.0e-4):
                # print band, i_band, i_q,i_direction, bf0[i_direction,i_q,i_band],mode_band_gruneisens[i_direction,i_q,i_band]
                if i_q == 0:
                    mode_band_gruneisens[
                        i_direction, i_q, i_band
                    ] = mode_band_gruneisens[i_direction, i_q + 1, i_band]
                elif i_q == len(direction) - 1:
                    mode_band_gruneisens[
                        i_direction, i_q, i_band
                    ] = mode_band_gruneisens[i_direction, i_q - 1, i_band]
                else:
                    mode_band_gruneisens[i_direction, i_q, i_band] = (
                        mode_band_gruneisens[i_direction, i_q - 1, i_band]
                        + mode_band_gruneisens[i_direction, i_q + 1, i_band]
                    ) / 2
                # print mode_band_gruneisens[i_direction,i_q,i_band]

temperature = np.logspace(
    np.log10(temperature_range[0]), np.log10(temperature_range[1]), 200
)
heat_capacity = []
thermal_expansion = []
gruneisen = []

for t in temperature:
    c = 2 * np.pi * HBAR * f0 / (2 * BOLTZMANN_K * t)
    c = BOLTZMANN_K * c ** 2 / np.sinh(c) ** 2

    c = c * weights[:, np.newaxis] / np.sum(weights)

    heat_capacity.append(np.sum(c))
    thermal_expansion.append(
        np.sum(c * mode_gruneisens) / (3.0 * bulk_at.get_volume() * bulk_modulus)
    )
    gruneisen.append(np.sum(c * mode_gruneisens) / np.sum(c))

temperature = [_f for _f in temperature]
heat_capacity = [_f for _f in heat_capacity]
thermal_expansion = [_f for _f in thermal_expansion]
gruneisen = [_f for _f in gruneisen]
mode_band_gruneisens = [_f.tolist() for _f in mode_band_gruneisens]

properties = {
    "a0": a0,
    "qha_temperature": temperature,
    "qha_heat_capacity": heat_capacity,
    "qha_thermal_expansion": thermal_expansion,
    "qha_gruneisen": gruneisen,
    "qha_band_distances": band_distances,
    "qha_mode_band_gruneisens": mode_band_gruneisens,
}
