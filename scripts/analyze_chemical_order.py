#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
from matplotlib.pyplot import *

import json
import ase.io
from ase.data import chemical_symbols
from analyze_utils import *
import math

(args, models, chem_order_tests, default_analysis_settings) = analyze_start('chemical_order_*')
print "got tests",chem_order_tests
ref_model_name = default_analysis_settings["ref_model"]

# read and parse all data
data = {}
for model_name in models:
    sys.stderr.write("reading data for model {}\n".format(model_name))
    cur_model_data = {}
    for chem_order_test_name in chem_order_tests:
        sys.stderr.write("   reading data for test {}\n".format(chem_order_test_name))

        # read chem_order test properties
        prop_filename ="{}-model-{}-test-{}-properties.json".format(args.test_set, model_name, chem_order_test_name)
        try:
            with open(prop_filename, "r") as model_data_file:
                json_data = json.load(model_data_file)
        except:
            sys.stderr.write("No properties file '{}'\n".format(prop_filename))
            continue

        cur_model_data[chem_order_test_name] = json_data.copy()

    if len(cur_model_data.keys()) > 0:
        print "got data for ",model_name, cur_model_data.keys()
        data[model_name] = cur_model_data.copy()

n_fig = 0
figure_nums = {}
for (model_i, model_name) in enumerate(sorted(data)):
    print "plot model",model_name, data[model_name]

    for chem_order_test_name in chem_order_tests:
        if chem_order_test_name not in figure_nums:
            n_fig += 1
            figure_nums[chem_order_test_name] = n_fig
        figure(figure_nums[chem_order_test_name])
        if chem_order_test_name not in data[model_name]:
            sys.stderr.write("skipping struct {} in plotting model {}\n".format(chem_order_test_name, model_name))
            continue

        if model_name != ref_model_name:
            if 'unrelaxed_energy_per_atom' in data[model_name][chem_order_test_name]:
                line, = plot( data[ref_model_name][chem_order_test_name]["unrelaxed_energy_per_atom"],
                              data[model_name][chem_order_test_name]["unrelaxed_energy_per_atom"],
                              ' ', color='C{}'.format(model_i), label="{} {} unrelaxed".format(model_name,chem_order_test_name) )
                line.set_marker('o')
                line.set_markerfacecolor("None")
                line.set_markersize(4.0)
            if 'relaxed_energy_per_atom' in data[model_name][chem_order_test_name]:
                line, = plot( data[ref_model_name][chem_order_test_name]["relaxed_energy_per_atom"],
                              data[model_name][chem_order_test_name]["relaxed_energy_per_atom"],
                              ' ', color='C{}'.format(model_i), label="{} {} relaxed".format(model_name,chem_order_test_name) )
                line.set_marker('o')
                line.set_markersize(4.0)

for chem_order_test_name in figure_nums:
    figure(figure_nums[chem_order_test_name])
    legend(loc="center left", bbox_to_anchor=[1, 0.5])
    xlabel("ref E (eV/atom)")
    ylabel("E (eV/atom)")
    x_range = ( min(data[ref_model_name][chem_order_test_name]["unrelaxed_energy_per_atom"]), max(data[ref_model_name][chem_order_test_name]["unrelaxed_energy_per_atom"]) ) 
    try:
        x_range = ( min(x_range[0], min(data[ref_model_name][chem_order_test_name]["relaxed_energy_per_atom"])), max(x_range[0], max(data[ref_model_name][chem_order_test_name]["relaxed_energy_per_atom"])))
    except:
        pass
    plot([x_range[0],x_range[1]], [x_range[0],x_range[1]], '--', color='black', label=None)
    savefig("{}_scatterplot.pdf".format(chem_order_test_name), bbox_inches='tight')
