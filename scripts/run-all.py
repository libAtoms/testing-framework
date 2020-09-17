#!/usr/bin/env python

import os
import sys
import glob
import re
import json
import argparse
import itertools

my_path=os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='Run all tests on all models.  Default command line options in dict of lists in "default_run_opts.json", keys are "global" or regexps matching model names')
parser.add_argument('--test_set', '-s', action='store', help='label for tests directories',required=True)
parser.add_argument('--models', '-m', action='store', nargs='+', type=str, help='models to include')
parser.add_argument('--tests', '-t', action='store', nargs='+', type=str, help='tests to include')
parser.add_argument('--omit-tests', '-o', action='store', nargs='+', type=str, help='tests to omit')
parser.add_argument('--force', '-f', action='store_true', help='force rerunning of tests')
parser.add_argument('--bugs', '-b', action='store_true', help='use bugs to generate job scripts')
parser.add_argument('--MPI', '-M', action='store_true', help='use MPI')
parser.add_argument('--OpenMP', '-O', action='store_true', help='use OpenMP')
parser.add_argument('--base_model', '-B', action='store', type=str, help='model to use as initial config for tests where it is enabled')
parser.add_argument('--models_path', '-P', action='store', type=str, help='path to models directory', default=os.path.join(os.getcwd(),'../models'))
parser.add_argument('--n_procs', '-N', action='store', type=int, help='number of processors to round to', default=16)

global_default_run_opts = []
default_run_opts = []
if os.path.exists("default_run_opts.json"):
    with open("default_run_opts.json","r") as f:
        default_run_opts = json.load(f)

    if "global" in default_run_opts:
        global_default_run_opts = default_run_opts["global"]

args = parser.parse_args(sys.argv[1:] + global_default_run_opts)

tests_path = os.path.join(my_path,"..","tests",args.test_set)

models = None
tests = None
omittests = None
if args.models is not None:
    print("Models asked for: ", args.models)
    models = list(itertools.chain.from_iterable([ glob.glob(os.path.join(args.models_path, d, 'model.py')) for d in args.models ]))
else:
    models = glob.glob(os.path.join(args.models_path, '*', 'model.py'))
print("Models path: ", args.models_path)
print("Models found: ", models)
models = [ os.path.split(d)[0] for d in models ]

if args.tests is not None:
    tests = list(itertools.chain.from_iterable([ glob.glob(os.path.join(tests_path, d, 'test.py')) for d in args.tests ]))
else:
    tests = glob.glob(os.path.join(tests_path, '*', 'test.py'))
tests = [ os.path.split(t)[0] for t in tests ]
if args.omit_tests is not None:
    for omit_test in list(itertools.chain.from_iterable([ glob.glob(os.path.join(tests_path, d)) for d in args.omit_tests ])):
        if omit_test in tests:
            tests.remove(omit_test)
        else:
            sys.stderr.write("WARNING: omitted test %s not in list of tests\n" % omit_test)


force = args.force
bugs = args.bugs

run_model_test = os.path.join(my_path,'run-model-test.py')

print('Models:', models)
print('Tests:', tests)

for model in models:
    model_name = os.path.split(model)[-1]

    model_default_run_opts = []
    if default_run_opts is not None: 
        for m in default_run_opts:
            if re.search(m, model_name) is not None:
                model_default_run_opts = default_run_opts[m]
                break

    args = parser.parse_args(sys.argv[1:] + model_default_run_opts + global_default_run_opts)

    for test in tests:
        test_name = os.path.split(test)[-1]
        try:
            fp = open(model+"/COMPUTATIONAL_COST")
            model_cost = float(fp.readline().strip())
        except:
            model_cost = 1.0
        try:
            fp = open(test+"/COMPUTATIONAL_COST")
            test_cost = float(fp.readline().strip())
        except:
            test_cost = 1.0
        np=int(round(int(test_cost*model_cost*args.n_procs)/float(args.n_procs)))*args.n_procs
        if np < args.n_procs:
            np = args.n_procs

        cmd_args = '{0} {1} {2} --test_set {3}'.format(run_model_test, "'"+os.path.join(args.models_path,model_name)+"'", test_name, args.test_set)
        if force:
            cmd_args += ' -force'
        if args.base_model is not None:
            cmd_args += ' --base_model '+args.base_model

        if args.bugs:
            ident_string= '{0}_{1}_{2}'.format(args.test_set, model_name, test_name)
            if args.MPI:
                bugs_script='test.bugs_script_mpi'
                mpi_cmd=''
            else:
                mpi_cmd='-no_mpirun'
                bugs_script='test.bugs_script'
            if args.OpenMP:
                bugs_script+='_openmp'
            bugs_script += '.'+os.environ['HOSTNAME'].split('.')[0]
            cmd=('env REDIRECT_IO="'+cmd_args+'" bugs -exec=python '+mpi_cmd+' -time=96h -np='+('%d' % np)+' -script='+bugs_script+
                 ' -name='+ident_string+' -fileroot='+ident_string+' -in_cwd -output=job.'+ident_string)
        else:
            cmd="python "+cmd_args
        print(cmd)
        # TODO: convert this to using subprocess.run()
        os.system(cmd)
        
