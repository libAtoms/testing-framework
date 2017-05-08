#!/usr/bin/env python

from analyze_utils import *

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

        print "data", data[model_name][test_name]
        if data[model_name][test_name]["dmu"] is None: # single component
            print "SURFACE", model_name, test_name, data[model_name][test_name]["Ef"]
        else: # multicomponent
#
#        label = data[model_name][test_name][0]
#
#        for defect_label in data[model_name][test_name]["defects"]:
#            # print "defect label",defect_label
#            defect = data[model_name][test_name]["defects"][defect_label]
#
#            ind = defect['atom_ind']
#            Z = defect['Z']
#            Ef = defect['Ef']
#            if 'dmu' in defect:
#                n_Z = defect['dmu'][0]
#                mu_Z = defect['dmu'][1]
#                if len(stable_mu_extrema) == 0:
#                    print "DEFECT", model_name, test_name, ind, Z, "(E_f0 =", Ef,") but no stable mu range exists"
#                else:
#                    mu_min = min([ mu_pt[mu_Z] for mu_pt in stable_mu_extrema] )
#                    mu_max = max([ mu_pt[mu_Z] for mu_pt in stable_mu_extrema] )
#                    print "DEFECT", model_name, test_name, "atom",ind, "Z",Z, "Ef ",Ef," + ( mu_{} = [".format(mu_Z),mu_min,"--",mu_max,"] ) = [",Ef+mu_min,"--",Ef+mu_max,"]"
#            else:
#                print "DEFECT", model_name, test_name, "atom",ind, "Z", Z, "Ef", Ef
#
        print ""
