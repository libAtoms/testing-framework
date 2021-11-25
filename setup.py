from setuptools import setup, find_packages

PACKAGENAME = "testingframework"
DESCRIPTION = "Module for testing various interatomic potentials"
AUTHOR = "Libatoms"
AUTHOR_EMAIL = ""

version = {}
with open("version.py") as fp:
    exec(fp.read(), version)

setup(
    name=PACKAGENAME,
    packages=find_packages(),
    version=version["__version__"],
    description=DESCRIPTION,
    long_description=open("README.md").read(),
    install_requires=["scipy", "numpy", "matplotlib", "pandas>=0.21.0", "scipy",],
    setup_requires=["pytest-runner",],
    tests_require=["pytest",],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    package_data={"": ["data/*", "calib/data/*"],},
)
