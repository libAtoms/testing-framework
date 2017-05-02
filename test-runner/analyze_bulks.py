#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
from matplotlib.pyplot import *
import json
import argparse
import glob, os
import ase.io
from ase.data import chemical_symbols

parser = argparse.ArgumentParser(description='Analyze bulk lattices')
parser.add_argument('--models_re', '-m', action='store', type=str, help='models to include', default='*')
parser.add_argument('--tests_re', '-t', action='store', type=str, help='tests to include', default='*')
parser.add_argument('--label', '-l', action='store', help='optional label for models/tests directories', default='')
args = parser.parse_args()

with open("default_analysis_settings.json") as default_analysis_file:
    default_analysis_settings = json.load(default_analysis_file)

with open("default_run_opts.json") as default_run_file:
    default_run_opts = json.load(default_run_file)
    if 'label' in default_run_opts and not args.label:
        args.label = default_run_opts['label']

models = glob.glob(os.path.join('..', 'models',args.label,args.models_re))
bulk_tests = glob.glob(os.path.join('..', 'tests',args.label,'bulk_'+args.tests_re))

ref_color = "black"
ref_symbol="-"
other_colors = [ "red", "blue", "cyan", "orange", "magenta" ] 
other_symbol="--"

if args.label is None:
    args.label = ""
else:
    args.label = args.label+"-"

# read and parse all data
element_min_Es = {}
data = {}
for model in models:
    print "reading data for model ",model
    model_name = os.path.basename(model)
    data[model_name] = {}
    cur_model_data = {}
    for bulk_test in bulk_tests:
        print "   reading data for test ", bulk_test
        bulk_test_name = os.path.basename(bulk_test).replace("bulk_","")

        # read bulk test structure
        struct_filename = "{}model-{}-test-{}-relaxed.xyz".format(args.label, model_name, "bulk_"+bulk_test_name)
        try:
            struct = ase.io.read(struct_filename, format="extxyz")
        except:
            sys.stderr.write("No struct file '{}'\n".format(struct_filename))
            continue

        # parse elements present
        elements_present = {}
        for Z in set(struct.get_atomic_numbers()):
            elements_present[Z] = sum(struct.get_atomic_numbers() == Z)

        # read bulk test properties
        prop_filename ="{}model-{}-test-{}-properties.json".format(args.label, model_name, "bulk_"+bulk_test_name)
        try:
            with open(prop_filename, "r") as model_data_file:
                json_data = json.load(model_data_file)
        except:
            sys.stderr.write("No properties file '{}'\n".format(prop_filename))
            continue

        cur_model_data[bulk_test_name] = json_data.copy()

        # shift energy by minimum energies of reference element structures
        E0 = 0.0
        for Z in elements_present:
            symb = chemical_symbols[Z]
            element_bulk_struct = default_analysis_settings["{}_ref_struct".format(symb)]
            try:
                E = element_min_Es[element_bulk_struct]
            except:
                with open("{}model-{}-test-bulk_{}_{}-properties.json".format(args.label, model_name, symb, element_bulk_struct), "r") as f:
                    d = json.load(f)
                E = min(d["E_vs_V"], key = lambda x : x[1])[1]
                element_min_Es[element_bulk_struct] = E
            E0 += sum(struct.get_atomic_numbers() == Z)*E
        E0 /= len(struct)

        # shift E_vs_V
        E_vs_V_orig = cur_model_data[bulk_test_name]["E_vs_V"]
        cur_model_data[bulk_test_name]["E_vs_V"] = [ [ EV[0], EV[1]-E0] for EV in E_vs_V_orig ]

    print "got data for ",model_name, cur_model_data
    data[model_name] = cur_model_data.copy()

ref_model_name = default_analysis_settings["ref_model"]
n_fig = 1
for model in models:
    model_name = os.path.basename(model)
    if model_name == ref_model_name:
        continue

    bulk_ind = -1
    figure(n_fig)
    for bulk_test in bulk_tests:
        bulk_ind += 1
        if bulk_test not in data[model_name]:
            continue

        line, = plot( [x[0] for x in data[model_name]["E_vs_V"]],  [x[1] - E0 for x in data[model_name]["E_vs_V"]], other_symbol)
        line.set_color(other_colors[bulk_ind])
        line, = plot( [x[0] for x in data[ref_model_name]["E_vs_V"]],  [x[1] - E0 for x in data[ref_model_name]["E_vs_V"]], ref_symbol)
        line.set_color(ref_color)


    savefig("{}_bulk.pdf".format(model_name))
    n_fig += 1
