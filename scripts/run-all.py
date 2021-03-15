#!/usr/bin/env python

import os
import sys
import subprocess
import glob
import re
import json
import argparse
import itertools
import shlex

my_path=os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='Run all tests on all models.  Default command line options in dict of lists in "default_run_opts.json", keys are "global" or regexps matching model names')
parser.add_argument('--test_set', '-s', nargs='+', help='label(s) for tests directories', required=True)
parser.add_argument('--models', '-m', action='store', nargs='+', type=str, help='models to include')
parser.add_argument('--tests', '-t', action='store', nargs='+', type=str, help='tests to include')
parser.add_argument('--omit-tests', '-o', action='store', nargs='+', type=str, help='tests to omit')
parser.add_argument('--force', '-f', action='store_true', help='force rerunning of tests')
parser.add_argument('--bugs', '-b', action='store_true', help='use bugs to generate job scripts')
parser.add_argument('--bugs_np', '-n', help='override number of processor when using bugs')
parser.add_argument('--bugs_host', '-H', help='host to force bugs to use',
                    default=re.sub('[0-9]*$', '', os.environ.get('HOSTNAME', 'None.FQDN').split('.')[0]))
parser.add_argument('--MPI', '-M', action='store_true', help='use MPI')
parser.add_argument('--OpenMP', '-O', action='store_true', help='use OpenMP')
parser.add_argument('--base_model', '-B', action='store', type=str, help='model to use as initial config for tests where it is enabled')
parser.add_argument('--models_path', '-MP', action='store', type=str, help='path to models directory', default=os.path.join(os.getcwd(),'../models'))
parser.add_argument('--tests_path', '-TP', action='store', type=str, help='path to tests directory (above test_set)', default=os.path.join(my_path,'../tests'))
parser.add_argument('--n_procs', '-N', action='store', type=int, help='number of processors to round to', default=16)

global_default_run_opts = []
default_run_opts = []
if os.path.exists("default_run_opts.json"):
    with open("default_run_opts.json","r") as f:
        default_run_opts = json.load(f)

    if "global" in default_run_opts:
        global_default_run_opts = default_run_opts["global"]

args = parser.parse_args(sys.argv[1:] + global_default_run_opts)


models = None
tests = None
if args.models is not None:
    print("Models asked for: ", args.models)
    models = list(itertools.chain.from_iterable([ glob.glob(os.path.join(args.models_path, d, 'model.py')) for d in args.models ]))
else:
    models = glob.glob(os.path.join(args.models_path, '*', 'model.py'))
print("Models path: ", args.models_path)
print("Models found: ", models)
models = [ os.path.split(d)[0] for d in models ]

def clean_path(p):
    if p.startswith(os.environ['HOME']):
        p_home = '"$HOME"'
        p_rest = p.replace(os.environ['HOME'], '')
    else:
        p_home = ''
        p_rest = p
    return p_home+shlex.quote(p_rest)

tests = []
for test_set in args.test_set:
    tests_path = os.path.join(args.tests_path,test_set)

    if args.tests is not None:
        tests_t = list(itertools.chain.from_iterable([ glob.glob(os.path.join(tests_path, d, 'test.py')) for d in args.tests ]))
    else:
        tests_t = glob.glob(os.path.join(tests_path, '*', 'test.py'))
    tests_t = [ os.path.split(t)[0] for t in tests_t ]
    if args.omit_tests is not None:
        for omit_test in list(itertools.chain.from_iterable([ glob.glob(os.path.join(tests_path, d)) for d in args.omit_tests ])):
            if omit_test in tests_t:
                tests_t.remove(omit_test)
            else:
                sys.stderr.write("WARNING: omitted test %s not in list of tests\n" % omit_test)
    tests.extend([(test, test_set) for test in tests_t])

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

    for test, test_set in tests:
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

        if args.bugs_np is None:
            np=int(round(int(test_cost*model_cost*args.n_procs)/float(args.n_procs)))*args.n_procs
            if np < args.n_procs:
                np = args.n_procs
        else:
            np = args.bugs_np

        if args.bugs:
            # only clean for bugs, which will feed path into shell
            cmd_args = [ clean_path(run_model_test), clean_path(os.path.join(args.models_path,model_name)),
                         clean_path(test) ]
        else:
            cmd_args = [ run_model_test, os.path.join(args.models_path,model_name), test ]
        cmd_args += [ '--system_label', test_set ]

        if force:
            cmd_args += [ '-force' ]
        if args.base_model is not None:
            cmd_args += [ '--base_model', args.base_model ]

        if args.bugs:
            ident_string= '{0}_{1}_{2}'.format(test_set, model_name, test_name)
            cmd = ['bugs', '-exec=python' ]
            if args.MPI:
                bugs_script='test.bugs_script_mpi'
            else:
                cmd += ['-no_mpirun']
                bugs_script='test.bugs_script'
            if args.OpenMP:
                bugs_script+='_openmp'
            bugs_script += '.'+args.bugs_host
            cmd += [ '-time=96h', f'-np={np}', f'-script={bugs_script}', f'-host={args.bugs_host}', f'-name={ident_string}',
                     f'-fileroot={ident_string}', '-in_cwd', f'-output=job.{ident_string}' ]

            cmd_env = os.environ.copy()
            cmd_env['REDIRECT_IO'] = ' '.join(cmd_args)
        else:
            cmd = ['python'] + cmd_args
            cmd_env = None
        print(cmd)
        # TODO: convert this to using subprocess.run()
        subprocess.run(cmd, env=cmd_env)

