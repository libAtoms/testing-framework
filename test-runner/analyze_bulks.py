#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
from matplotlib.pyplot import *
import json
import argparse
import glob, os
import ase.io
from ase.data import chemical_symbols
from analyze_utils import *

parser = argparse.ArgumentParser(description='Analyze bulk lattices')
parser.add_argument('--models_re', '-m', action='store', type=str, help='models to include', default='*')
parser.add_argument('--tests_re', '-t', action='store', type=str, help='tests to include', default='*')
parser.add_argument('--label', '-l', action='store', help='optional label for models/tests directories', default='')
args = parser.parse_args()

try:
    defaults_label = os.environ["DEFAULTS_LABEL"]+"_"
except:
    defaults_label = ""

with open("{}default_analysis_settings.json".format(defaults_label)) as default_analysis_file:
    default_analysis_settings = json.load(default_analysis_file)

with open("{}default_run_opts.json".format(defaults_label)) as default_run_file:
    default_run_opts = json.load(default_run_file)
    if 'label' in default_run_opts and not args.label:
        args.label = default_run_opts['label']

models = [os.path.basename(f) for f in glob.glob(os.path.join('..', 'models',args.label,args.models_re))]
bulk_tests = [os.path.basename(f).replace("bulk_","") for f in glob.glob(os.path.join('..', 'tests',args.label,'bulk_'+args.tests_re))]

ref_symbol="-"
other_symbol="--"
struct_colors = [ "red", "blue", "cyan", "orange", "magenta" ] 

if args.label is None:
    args.label = ""
else:
    args.label = args.label+"-"

ref_struct_data = get_ref_structs(args.label, models, default_analysis_settings)

# read and parse all data
element_min_Es = {}
data = {}
for model_name in models:
    print "reading data for model ",model_name
    element_min_Es[model_name] = {}
    data[model_name] = {}
    cur_model_data = {}
    for bulk_test_name in bulk_tests:
        print "   reading data for test ", bulk_test_name

        # read bulk test structure
        struct_filename = "{}model-{}-test-{}-relaxed.xyz".format(args.label, model_name, "bulk_"+bulk_test_name)
        try:
            struct = ase.io.read(struct_filename, format="extxyz")
        except:
            sys.stderr.write("No struct file '{}'\n".format(struct_filename))
            continue

        # read bulk test properties
        prop_filename ="{}model-{}-test-{}-properties.json".format(args.label, model_name, "bulk_"+bulk_test_name)
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
            element_bulk_struct = default_analysis_settings["{}_ref_struct".format(symb)]
            E = ref_struct_data["min_Es"][model_name][element_bulk_struct]
            E0 += sum(struct.get_atomic_numbers() == Z)*E
        E0 /= len(struct)

        # shift E_vs_V
        E_vs_V_orig = cur_model_data[bulk_test_name]["E_vs_V"]
        cur_model_data[bulk_test_name]["E_vs_V"] = [ [ EV[0], EV[1]-E0] for EV in E_vs_V_orig ]

    # print "got data for ",model_name, cur_model_data.keys()
    data[model_name] = cur_model_data.copy()

ref_model_name = default_analysis_settings["ref_model"]
n_fig = 1
for model_name in models:
    for bulk_test_name in bulk_tests:
        min_EV  = min(data[model_name][bulk_test_name]["E_vs_V"], key = lambda x : x[1])
        print "BULK_E_V",model_name,bulk_test_name, min_EV[0], min_EV[1]
    if model_name == ref_model_name:
        continue

    bulk_ind = -1
    figure(n_fig)
    for bulk_test_name in bulk_tests:
        bulk_ind += 1
        if bulk_test_name not in data[model_name]:
            sys.stderr.write("skipping struct {} in plotting model {}\n".format(bulk_test_name, model_name))
            continue

        line, = plot( [x[0] for x in data[model_name][bulk_test_name]["E_vs_V"]], [x[1] for x in data[model_name][bulk_test_name]["E_vs_V"]], other_symbol)
        line.set_color(struct_colors[bulk_ind])
        line, = plot( [x[0] for x in data[ref_model_name][bulk_test_name]["E_vs_V"]], [x[1] for x in data[ref_model_name][bulk_test_name]["E_vs_V"]], ref_symbol, label=bulk_test_name)
        line.set_color(struct_colors[bulk_ind])

    legend(loc="center left", bbox_to_anchor=[1, 0.5])
    xlabel("V ($A^3$/atom)")
    ylabel("E (eV/atom)")
    savefig("{}_bulk.pdf".format(model_name), bbox_inches='tight')
    n_fig += 1
