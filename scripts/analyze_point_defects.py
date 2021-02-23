#!/usr/bin/env python

from analyze_utils import *

(args, models, tests, default_analysis_settings) = analyze_start(['point_defect_*'])

try:
    (mcc_compositions, mcc_energies) = get_multicomponent_constraints(args.test_set, models, default_analysis_settings["multicomponent_constraints"])
except:
    (mcc_compositions, mcc_energies) = (None, None)

# print("multicomponent_constraints_data ", multicomponent_constraints_data)

from multicomponent_mu_range import mu_range

# read and parse all data
data = read_properties(models, tests, args.test_set)

n_fig = 0
defect_table_data = {}
for model_name in models:
    table_entry = []
    for test_name in tests:
        print("DO {} {}".format(model_name, test_name))

        try:
            (cur_min_EV, cur_composition) = read_ref_bulk_model_struct(args.test_set, model_name, data[model_name][test_name]["bulk_struct_test"])
        except:
            print("No data")
            continue
        try:
            (stable_mu_extrema, full_mu_range) = mu_range(cur_min_EV, cur_composition, data[model_name][test_name]["bulk_struct_test"], mcc_compositions, mcc_energies[model_name])
        except:
            (stable_mu_extrema, full_mu_range) = (None,None)
        # print(model_name, test_name, "stable_mu_extrema", stable_mu_extrema)
        # print("full_mu_range", full_mu_range)
        if stable_mu_extrema is not None:
            if len(stable_mu_extrema) > 0:
                print("stable mu range:")
                for pt in stable_mu_extrema:
                    for Z in sorted(pt.keys()):
                        print("mu_{} = {}".format(Z, pt[Z]),end='')
                    print("")
            else:
                print("stable mu range: None")

        for defect_label in data[model_name][test_name]["defects"]:
            # print("defect label",defect_label)
            defect = data[model_name][test_name]["defects"][defect_label]

            ind = defect['atom_ind']
            Z = defect['Z']
            Ef = defect['Ef']
            Ef0 = defect['Ef0']
            if 'dmu' in defect:
                n_Z = defect['dmu'][0]
                mu_Z = defect['dmu'][1]
                if len(stable_mu_extrema) == 0:
                    print(f'DEFECT {model_name}, {test_name}, {ind}, {Z}, (Ef0 = {Ef0:.4f} Ef = {Ef:.4f}) but no stable mu range exists')
                else:
                    mu_min = min([ mu_pt[mu_Z] for mu_pt in stable_mu_extrema] )
                    mu_max = max([ mu_pt[mu_Z] for mu_pt in stable_mu_extrema] )
                    print("DEFECT", model_name, test_name, "atom",ind, "Z",Z)
                    l1 = f'{Ef0:.4f} + ( mu_{mu_Z} = [ {mu_min:.4f} -- {mu_max:.4f} ] ) = [ {Ef0+mu_min:.4f} -- {Ef0+mu_max:.4f} ]'
                    print('DEFECT Ef0', l1)
                    l2 = f'{Ef:.4f} + ( mu_{mu_Z} = [ {mu_min:.4f} -- {mu_max:.4f} ] ) = [ {Ef+mu_min:.4f} -- {Ef+mu_max:.4f} ]'
                    print('DEFECT Ef', l2)
                    table_entry.append((test_name, ind, Z, l1, l2))
            else:
                l1 = f'{Ef0:.4f}'
                l2 = f'{Ef:.4f}'
                print("DEFECT", model_name, test_name, "atom",ind, "Z", Z, "Ef0", l1, "Ef", l2)
                table_entry.append((test_name, ind, Z, l1, l2))

        for (test_name, ind, Z, Ef0, Ef) in table_entry:
            l = test_name+" "+str(Z)+" "+str(ind)
            if l not in defect_table_data:
                defect_table_data[l] = {}
            defect_table_data[l][model_name] = (Ef0,Ef)

        print("")

defects_sorted = sorted(defect_table_data.keys())
defects_sorted_print = [ n.replace("point_defect_","").replace("_"," ") for n in defects_sorted ]
print("\\begin{tabular}{ l & " + " & ".join(["c"]*2*len(defects_sorted)) + " }")
print("model & "+" & ".join( [ "\\multicolumn{{2}}{{c}}{{{}}}".format(n) for n in defects_sorted_print ]) + "\\\\")
print("      & "+" & ".join( [" E$_f^0$ & E$_f$ "]*len(defects_sorted) ),end='')
for model_name in models:
    print("\\\\\n"+model_name.replace("_","\\_"),end='')
    for d in defects_sorted:
        if model_name in defect_table_data[l]:
            print(" & " + "{}".format(defect_table_data[d][model_name][0])+" & "+"{}".format(defect_table_data[d][model_name][1]),end='')
        else:
            print(" & & ",end='')
print("\n\\end{tabular}")
