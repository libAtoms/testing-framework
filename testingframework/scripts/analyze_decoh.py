#!/usr/bin/env python

from analyze_utils import *
import numpy as np
from ase.units import _k
import matplotlib.pyplot as plt

(args, models, tests, default_analysis_settings) = analyze_start("surface-decohesion-*")

try:
    (mcc_compositions, mcc_energies) = get_multicomponent_constraints(
        args.test_set, models, default_analysis_settings["multicomponent_constraints"]
    )
except:
    (mcc_compositions, mcc_energies) = (None, None)

# print("multicomponent_constraints_decoh_data ", multicomponent_constraints_decoh_data)

from multicomponent_mu_range import mu_range

ref_linestyles = ["-", "--"]
other_linestyles = [":", "-."]
struct_colors = [
    "black",
    "red",
    "blue",
    "orange",
    "green",
    "brown",
    "grey",
    "magenta",
    "cyan",
]
ref_model_name = default_analysis_settings["ref_model"]

# read and parse all decoh_data
decoh_data = read_properties(models, tests, args.test_set)

for test in tests:
    ref_opening = decoh_data[ref_model_name][test][
        "surface_decohesion_unrelaxed_opening"
    ]
    ref_energy = decoh_data[ref_model_name][test]["surface_decohesion_unrelaxed_energy"]
    ref_stress = decoh_data[ref_model_name][test]["surface_decohesion_unrelaxed_stress"]

    # f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 6))

    ax1.plot(
        ref_opening,
        ref_energy,
        label=ref_model_name,
        linestyle=ref_linestyles[0],
        color="black",
    )
    # ax1.set_xlabel("Unrelaxed opening ($\AA$)")
    ax1.set_ylabel("Energy (eV)")
    ax2.plot(
        ref_opening,
        ref_stress,
        label=ref_model_name,
        linestyle=ref_linestyles[0],
        color="black",
    )
    ax2.set_xlabel("Unrelaxed opening ($\AA$)")
    ax2.set_ylabel(r"Stress (eV/$\AA$)")

    i = 0
    for model in sorted(decoh_data.keys()):
        if bool(decoh_data[model]) is True and model != ref_model_name:
            i += 1
            ax1.plot(
                decoh_data[model][test]["surface_decohesion_unrelaxed_opening"],
                decoh_data[model][test]["surface_decohesion_unrelaxed_energy"],
                label=model,
                color=struct_colors[i],
            )  # , linestyle=other_linestyles[0])
            ax2.plot(
                decoh_data[model][test]["surface_decohesion_unrelaxed_opening"],
                decoh_data[model][test]["surface_decohesion_unrelaxed_stress"],
                label=model,
                color=struct_colors[i],
            )  # , linestyle=other_linestyles[0])

    # plt.title(test)
    ax1.legend()
    plt.tight_layout()
    plt.savefig("{}.pdf".format(test), bbox_inches="tight")  ####
