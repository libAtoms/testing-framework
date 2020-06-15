#!/usr/bin/env python

import numpy as np
from analyze_utils import *
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot

(args, models, tests, default_analysis_settings) = analyze_start(['phonon_DOS_bulks'])

data = read_properties(models, tests, args.test_set)

ref_model_name = default_analysis_settings["ref_model"]

n_bulks = len(data[ref_model_name]["phonon_DOS_bulks"])

dE_rms_list = {}
dF_rms_list = {}
dS_rms_list = {}
print("PHONON RESULTS")
for model_name in models:
    if model_name == ref_model_name:
        continue
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    for (bulk_i, bulk_struct_test) in enumerate(data[ref_model_name]["phonon_DOS_bulks"]):
        ref_model_data = data[ref_model_name]["phonon_DOS_bulks"][bulk_struct_test]
        try:
            model_data = data[model_name]["phonon_DOS_bulks"][bulk_struct_test]
        except:
            continue

        ax.plot(ref_model_data['DOS'][0], ref_model_data['DOS'][1], "-", color='C{}'.format(bulk_i), label=bulk_struct_test)
        if model_name != ref_model_name:
            ax.plot(model_data['DOS'][0], model_data['DOS'][1], ":", color='C{}'.format(bulk_i), label=None)

    ax.set_xlabel("freq (cm$^{-1}$)")
    ax.set_ylabel("DOS (arb. units)")
    xlim = ax.get_xlim()
    # ax.set_xlim(0.0, xlim[0])
    ax.legend()

    fig.savefig("phonon_DOS_bulks-"+model_name+".pdf")
    pyplot.clf()

    print("")
