#!/usr/bin/env python

from analyze_utils import *
import numpy as np
from ase.units import _k
import matplotlib.pyplot as plt

(args, models, tests, default_analysis_settings) = analyze_start('RSS')

try:
    (mcc_compositions, mcc_energies) = get_multicomponent_constraints(args.test_set, models, default_analysis_settings["multicomponent_constraints"])
except:
    (mcc_compositions, mcc_energies) = (None, None)

# print("multicomponent_constraints_data ", multicomponent_constraints_data)

from multicomponent_mu_range import mu_range

ref_linestyles=[ "-", "--" ]
other_linestyles=[ ":", "-." ]
struct_colors = [ "red", "blue", "orange", "green", "brown", "grey", "magenta","cyan" ]
ref_model_name = default_analysis_settings["ref_model"]

# read and parse all data
rss_data = read_properties(models, tests, args.test_set)

print(rss_data)

ace_name = "ACE_B8_N4_18_lap_1.1"

print(len(rss_data[ace_name]['RSS']["volumes"]))
print(len(rss_data['GAP']['RSS']["volumes"]))

plt.scatter(rss_data[ace_name]['RSS']["volumes"], rss_data[ace_name]['RSS']["energies"], label="ACE")
plt.scatter(rss_data["GAP"]['RSS']["volumes"], rss_data["GAP"]["RSS"]["energies"], label="GAP")
plt.legend()
plt.xlim(10, 75)
plt.ylim(-163.4, -161.5)
plt.savefig("RSS_random.pdf")

# f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(5,7))
#
# a0.plot(A, B, '_', marker='o', color="grey", markersize=2, label="CASTEP_unrelaxed")
# a0.plot(C, D, '_', marker='o', color="black", markersize=2, label="CASTEP_relaxed")
# L = list(zip(list(zip(A,B)),list(zip(C,D))))
# lc = LineCollection(L, colors=["grey"], linewidth=0.5)
# a0.add_collection(lc)
#
# K = []
# K.append(D)
# histcolors = ["black"]
#
# for (i,model) in enumerate(models):
#     if model != ref_model_name:
#         C, D = rss_data[model]['gap_6_rss_like_open_ended_minima_rerelax_sd2']["V_relaxed"], rss_data[model]['gap_6_rss_like_open_ended_minima_rerelax_sd2']["E_relaxed"]
#         L = list(zip(list(zip(A,B)),list(zip(C,D))))
#         lc = LineCollection(L, colors=["grey"], linewidth=0.5)
#         a0.add_collection(lc)
#         a0.plot(C, D, '_', marker='o', color=struct_colors[i], markersize=2, label=model)
#         histcolors.append(struct_colors[i])
#         K.append(D)
#         #a1.hist(D, color=struct_colors[i], label=model, bins=15, alpha=0.8)
#
# #a1.hist(D, color="black", label="CASTEP_relaxed", bins=15, alpha=0.8)
# a1.hist(K, color=histcolors, bins=18)
# a1.set_xlabel("Energy (ev/Atom)")
# a1.set_ylabel("Occurence")
# a0.set_ylabel("Energy (ev/Atom)")
# a0.set_xlabel("Volume [A^3/Atom]")
#
# a0.legend()
# plt.xlabel("Volumes [A^3]")
# plt.ylabel("Energy [eV/atom]")
# plt.xlim(10, 75)
# plt.ylim(-164.5, -161.5)
