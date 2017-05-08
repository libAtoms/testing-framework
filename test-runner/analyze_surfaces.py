#!/usr/bin/env python

from analyze_utils import *
import numpy as np

(args, models, tests, default_analysis_settings) = analyze_start('*surface*')

(mcc_compositions, mcc_energies) = get_multicomponent_constraints(args.label, models, default_analysis_settings["multicomponent_constraints"])

# print "multicomponent_constraints_data ", multicomponent_constraints_data

from multicomponent_mu_range import mu_range

# read and parse all data
data = read_properties(models, tests, args.label)

n_fig = 0
for model_name in models:
    for test_name in tests:
        print "DO {} {}".format(model_name, test_name)

        (cur_min_EV, cur_composition) = read_ref_bulk_model_struct(args.label, model_name, data[model_name][test_name]["bulk_struct"])
        (stable_mu_extrema, full_mu_range) = mu_range(cur_min_EV, cur_composition, data[model_name][test_name]["bulk_struct"], mcc_compositions, mcc_energies[model_name])
        # print model_name, test_name, "stable_mu_extrema", stable_mu_extrema
        # print "full_mu_range", full_mu_range
        if stable_mu_extrema is not None:
            if len(stable_mu_extrema) > 0:
                print "stable mu range:"
                for pt in stable_mu_extrema:
                    for Z in sorted(pt.keys()):
                        print "mu_{} = {}".format(Z, pt[Z]),
                    print ""
            else:
                print "stable mu range: None"

        # print "data", data[model_name][test_name]
        if data[model_name][test_name]["dmu"] is None: # single component
            print "SURFACE", model_name, test_name, data[model_name][test_name]["Ef"]
        else: # multicomponent
            print "SURFACE", model_name, test_name, data[model_name][test_name]["Ef"],"+ ",
            mu_contrib_min_total = 0.0
            mu_contrib_max_total = 0.0
            for mu_Z in data[model_name][test_name]["dmu"]:
                n_dmu = data[model_name][test_name]["dmu"][mu_Z]
                if n_dmu != 0:
                    print "( {} * mu_{} =".format(n_dmu,mu_Z),
                    mu_min = min([ mu_pt[int(mu_Z)] for mu_pt in stable_mu_extrema] )
                    mu_max = max([ mu_pt[int(mu_Z)] for mu_pt in stable_mu_extrema] )
                    mu_contrib_min = np.finfo(float).max
                    mu_contrib_max = -np.finfo(float).max
                    mu_contrib_min = min(mu_contrib_min, n_dmu*mu_min)
                    mu_contrib_min = min(mu_contrib_min, n_dmu*mu_max)
                    mu_contrib_max = max(mu_contrib_max, n_dmu*mu_min)
                    mu_contrib_max = max(mu_contrib_max, n_dmu*mu_max)
                    print "[",mu_contrib_min,"--",mu_contrib_max,"] )",
                    mu_contrib_min_total += mu_contrib_min
                    mu_contrib_max_total += mu_contrib_max
            print " = [", data[model_name][test_name]["Ef"] + mu_contrib_min_total,"--",data[model_name][test_name]["Ef"] + mu_contrib_max_total,"]"
        print ""
