#!/usr/bin/env python

from analyze_utils import *
import numpy as np
from ase.units import _k
import matplotlib.pyplot as plt

(args, models, tests, default_analysis_settings) = analyze_start('qha_diamond')

try:
    (mcc_compositions, mcc_energies) = get_multicomponent_constraints(args.test_set, models, default_analysis_settings["multicomponent_constraints"])
except:
    (mcc_compositions, mcc_energies) = (None, None)

# print("multicomponent_constraints_data ", multicomponent_constraints_data)

from multicomponent_mu_range import mu_range

ref_linestyles=[ "-", "--" ]
other_linestyles=[ ":", "-." ]
struct_colors = [ "black", "red", "blue", "orange", "green", "brown", "grey", "magenta","cyan" ]
ref_model_name = default_analysis_settings["ref_model"]

# read and parse all data
qha_data = read_properties(models, tests, args.test_set)

fig, axs = plt.subplots(3,1, sharex=False, sharey=False,figsize=(6.6,9.9))

axs[0].plot(qha_data[ref_model_name]['qha_diamond']['qha_temperature'], qha_data[ref_model_name]['qha_diamond']['qha_heat_capacity'], label=ref_model_name, linestyle=ref_linestyles[0], color="black")
axs[1].plot(qha_data[ref_model_name]['qha_diamond']['qha_temperature'], qha_data[ref_model_name]['qha_diamond']['qha_gruneisen'], label=ref_model_name, linestyle=ref_linestyles[0], color="black")
axs[2].plot(qha_data[ref_model_name]['qha_diamond']['qha_temperature'], qha_data[ref_model_name]['qha_diamond']['qha_thermal_expansion'], label=ref_model_name, linestyle=ref_linestyles[0], color="black")

i = 0
for model in sorted(qha_data.keys()):
    if bool(qha_data[model]) is True and model != ref_model_name:
        i += 1
        axs[0].plot(qha_data[model]['qha_diamond']['qha_temperature'], qha_data[model]['qha_diamond']['qha_heat_capacity'], label=model, color=struct_colors[i])
        axs[1].plot(qha_data[model]['qha_diamond']['qha_temperature'], qha_data[model]['qha_diamond']['qha_gruneisen'], label=model, color=struct_colors[i])
        axs[2].plot(qha_data[model]['qha_diamond']['qha_temperature'], qha_data[model]['qha_diamond']['qha_thermal_expansion'], label=model, color=struct_colors[i])

axs[0].set_ylabel("Heat Capacity [kb]")
axs[1].set_ylabel("Gruneisen Parameter")
axs[2].set_ylabel("Thermal Expansion [1E-6/K]")
axs[2].set_xlabel("Temperature [K]")
plt.legend()
plt.savefig("./qha.pdf")
