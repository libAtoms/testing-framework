#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
from matplotlib.pyplot import *

import json
import ase.io
from ase.data import chemical_symbols
from analyze_utils import *
import math

(args, models, bulk_tests, default_analysis_settings) = analyze_start('bulk_*')

ref_linestyles=[ "-", "--" ]
other_linestyles=[ ":", "-." ]
struct_colors = [ "black", "red", "blue", "cyan", "orange", "magenta", "green", "grey", "brown" ]

element_ref_struct_data = get_element_ref_structs(args.label, models, default_analysis_settings["element_ref_struct"])

# read and parse all data
data = {}
struct_data = {}
for model_name in models:
    sys.stderr.write("reading data for model {}\n".format(model_name))
    data[model_name] = {}
    cur_model_data = {}
    for bulk_test_name in bulk_tests:
        sys.stderr.write("   reading data for test {}\n".format(bulk_test_name))

        # read bulk test structure
        struct_filename = "{}model-{}-test-{}-relaxed.xyz".format(args.label, model_name, bulk_test_name)
        try:
            struct = ase.io.read(struct_filename, format="extxyz")
        except:
            sys.stderr.write("No struct file '{}'\n".format(struct_filename))
            continue

        # read bulk test properties
        prop_filename ="{}model-{}-test-{}-properties.json".format(args.label, model_name, bulk_test_name)
        try:
            with open(prop_filename, "r") as model_data_file:
                json_data = json.load(model_data_file)
        except:
            sys.stderr.write("No properties file '{}'\n".format(prop_filename))
            continue

        cur_model_data[bulk_test_name] = json_data.copy()

        # parse elements present
        elements_present = {}
        for Z in set(struct.get_atomic_numbers()):
            elements_present[Z] = sum(struct.get_atomic_numbers() == Z)

        # shift energy by minimum energies of reference element structures
        E0 = 0.0
        for Z in elements_present:
            symb = chemical_symbols[Z]
            print "looking for min E ",symb,model_name
            try:
                E = element_ref_struct_data[symb]["min_Es"][model_name]
                E0 += sum(struct.get_atomic_numbers() == Z)*E
            except:
                pass
        E0 /= len(struct)

        # shift E_vs_V
        E_vs_V_orig = cur_model_data[bulk_test_name]["E_vs_V"]
        cur_model_data[bulk_test_name]["E_vs_V"] = [ [ EV[0], EV[1]-E0] for EV in E_vs_V_orig ]

        if not bulk_test_name in struct_data:
            struct_data[bulk_test_name] = {}
        struct_data[bulk_test_name]["formula_unit"] = formula_unit(struct.get_atomic_numbers())

    # print "got data for ",model_name, cur_model_data.keys()
    data[model_name] = cur_model_data.copy()

ref_model_name = default_analysis_settings["ref_model"]
n_fig = 1
for model_name in models:
    figure_nums = {}
    bulk_inds = {}
    for bulk_test_name in bulk_tests:
        try:
            min_EV  = min(data[model_name][bulk_test_name]["E_vs_V"], key = lambda x : x[1])
        except:
            print "no data for",model_name,bulk_test_name
        continue
        print "BULK_E_V_MIN",model_name,bulk_test_name, min_EV[0], min_EV[1]

    for bulk_test_name in bulk_tests:
        fu = ""
        for (Z, i) in struct_data[bulk_test_name]["formula_unit"]:
            if i == 1:
                fu += format(ase.data.chemical_symbols[Z])
            else:
                fu += "{}{}".format(ase.data.chemical_symbols[Z], i)
        if fu not in figure_nums:
            figure_nums[fu] = n_fig
            n_fig += 1
        figure(figure_nums[fu])
        if fu not in bulk_inds:
            bulk_inds[fu] = 0
        else:
            bulk_inds[fu] += 1
        if bulk_test_name not in data[model_name]:
            sys.stderr.write("skipping struct {} in plotting model {}\n".format(bulk_test_name, model_name))
            continue

        ref_linestyle = ref_linestyles[int(math.floor(bulk_inds[fu]/len(struct_colors)))]
        other_linestyle = other_linestyles[int(math.floor(bulk_inds[fu]/len(struct_colors)))]
        color = struct_colors[bulk_inds[fu] % len(struct_colors)]
        if model_name != ref_model_name:
            line, = plot( [x[0] for x in data[model_name][bulk_test_name]["E_vs_V"]], [x[1] for x in data[model_name][bulk_test_name]["E_vs_V"]], other_linestyle)
            line.set_color(color)
            line.set_marker('o')
            line.set_markersize(2.5)
        try:
            line, = plot( [x[0] for x in data[ref_model_name][bulk_test_name]["E_vs_V"]], [x[1] for x in data[ref_model_name][bulk_test_name]["E_vs_V"]], ref_linestyle, label=bulk_test_name)
        except:
            print "no data for struct",bulk_test_name,"ref model",ref_model_name
            pass
        line.set_color(color)

    for fu in figure_nums:
        figure(figure_nums[fu])
        legend(loc="center left", bbox_to_anchor=[1, 0.5])
        xlabel("V ($A^3$/atom)")
        ylabel("E (eV/atom)")
        savefig("{}_{}_bulk.pdf".format(model_name,fu), bbox_inches='tight')
