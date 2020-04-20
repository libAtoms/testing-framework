Testing Framework for atomistic models
======================================

The purpose of this repository is to aid the testing of a large number
of interatomic potentials for a variety of systems (materials or
molecules). It uses the [Atomic Simulation Environment (ASE)](https://gitlab.com/ase/ase) to glue
things together. Although the relevant tests might differ between
different systems, there are a lot of commonalities, e.g. calculating
response of various crystal phases to deformation, or the evaluation
of a series of test configurations. Coding up the tests in a
model-agnostic way ensures that exactly the same tests are being
performed for different models.


Structure of the framework
--------------------------

The following subdirectories comprise the testing framework:

- `share` : routines that are frequently used in many different
  tests, e.g. relaxation of geometry, calculation of E-V curves,
  energy and force evaluations for a series of configurations
- `scripts` : top level scripts that are used to initiate tests
- `tests` : this is where the actual test routines are specified,
  different systems (e.g. materials) each having their own
  subdirectory

Any interatomic potential that can be instantiated as an ASE
 calculator can be tested. The potential models (each in its own
 directory) are kept in a directory structure separate from the
 testing-framework structure above. Examples are given under the
 `example_run_dir` directory, the first level subdirectories specify
 the system, and under each are `models` and `run_dir` subdirectories,
 the former contains the different potentials, the latter is where
 the tests are actually run and the test results appear.

Running tests
-------------

The following is an example of how tests are run:

```
cd testing-framework/example_run_dir/CSiGe/run_dir
../../../scripts/run-model-test.py -s CSiGe Tersoff bulk_Si_diamond
```

This will produce three things:

- The result of the test, in
  `CSiGe-model-Tersoff-test-bulk_Si_diamond-properties.json`
- The log of the standard output during the test in
  `CSiGe-model-Tersoff-test-bulk_Si_diamond.txt`
- A directory called `run_CSiGe-model-Tersoff-test-bulk_Si_diamond` in
  which the intermediate files generated during the test are kept, to
  aid debugging

Running the same test again will first check to see if the `.json`
file exists, and if it does, the corresponding test is skipped.

The `run-all.py` script will find all tests and all models under a
given system and run all those without `.json` files present in the
`run_dir`.


Specification of a model
------------------------

An interatomic potential model is specified by a `model.py` file,
which has to instantiate an ASE Calculator with the variable `calculator`
and give a label to it in the variable `name`.

Specification of a test
-----------------------

Every test has to be called `test.py` and has to create a dictionary in a variable `properties`, and can be as simple as one line, e.g.
```
properties = lattice.do_lattice(os.path.abspath(os.path.dirname(__file__)), 'cubic')
```

This one above is using the `do_lattice` function from the `share`
directory. This routine loads up a file called `bulk.xyz` from the same test directory 
(this is why the full path name of the `test.py` script has to be
passed as the first argument), and calculates the Energy-Volume curve
of that structure, and also computes the elastic constants. 
