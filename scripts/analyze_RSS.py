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
struct_colors = [ "black", "red", "blue", "orange", "green", "brown", "grey", "magenta","cyan" ]
ref_model_name = default_analysis_settings["ref_model"]

# read and parse all data
rss_data = read_properties(models, tests, args.test_set)

print(rss_data)

plt.hist(rss_data["GAP"]["RSS"]["energies"], bins=20)
plt.hist(rss_data["ACE_B14_N7_18_lc_rrqr_1.0e-10"]["RSS"]["energies"], bins=20)
plt.savefig("RSS.pdf")
