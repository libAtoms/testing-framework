#!/usr/bin/env python3

from analyze_utils import *
import numpy as np
from ase.units import _k
import matplotlib.pyplot as plt

(args, models, tests, default_analysis_settings) = analyze_start("Bain_path_*")

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

# read and parse all Bain_data
Bain_data = read_properties(models, tests, args.test_set)

for test in tests:
    ref_Bain = Bain_data[ref_model_name][test]["Bain_path_E"]

    f, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        [v[0] for v in ref_Bain],
        [v[1] for v in ref_Bain],
        label=ref_model_name,
        linestyle=ref_linestyles[0],
        color="black",
    )
    # ax1.set_xlabel("Unrelaxed opening ($\AA$)")
    ax.set_xlabel("c/a")
    ax.set_ylabel("Energy (eV/atom)")

    i = 0
    for model in sorted(Bain_data.keys()):
        if bool(Bain_data[model]) is True and model != ref_model_name:
            i += 1
            ax.plot(
                [v[0] for v in Bain_data[model][test]["Bain_path_E"]],
                [v[1] for v in Bain_data[model][test]["Bain_path_E"]],
                label=model,
                color=struct_colors[i],
                linestyle=other_linestyles[0],
            )

    # plt.title(test)
    ax.legend()
    f.tight_layout()
    f.savefig("{}.pdf".format(test), bbox_inches="tight")  ####
