#!/usr/bin/python2

from quippy import *
import numpy as np

sc_size = 25.0  # size of supercell
angle = np.radians(90.0)

aa = AtomsList()

for h in np.linspace(-1.0, 4.0, 101):
    a = Atoms()
    a.set_lattice(make_lattice(sc_size, sc_size, sc_size, angle, angle, angle))
    f = 2.2 / 2.0
    a.add_atoms([f, 0.0, -f / np.sqrt(2.0)], 15)
    a.add_atoms([-f, 0.0, -f / np.sqrt(2.0)], 15)
    a.add_atoms([0.0, f, f / np.sqrt(2.0) + h], 15)
    a.add_atoms([0.0, -f, f / np.sqrt(2.0) + h], 15)
    a.center()
    aa.append(a)

aa.write("P2_dimer_scan_C2v.xyz")

aa = AtomsList()

for h in np.linspace(0.0, 5.0, 101):
    a = Atoms()
    a.set_lattice(make_lattice(sc_size, sc_size, sc_size, angle, angle, angle))
    f = 2.2 / np.sqrt(8.0 / 3)
    d0 = f * 4.0 / 3 + 1.0
    a.add_atoms([f * np.sqrt(8.0 / 9), 0.0, -f / 3], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), f * np.sqrt(2.0 / 3), -f / 3], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), -f * np.sqrt(2.0 / 3), -f / 3], 15)
    a.add_atoms([0.0, 0.0, f], 15)

    a.add_atoms([f * np.sqrt(8.0 / 9), 0.0, -f / 3 + d0 + h], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), f * np.sqrt(2.0 / 3), -f / 3 + d0 + h], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), -f * np.sqrt(2.0 / 3), -f / 3 + d0 + h], 15)
    a.add_atoms([0.0, 0.0, f + d0 + h], 15)

    a.center()
    aa.append(a)

aa.write("P4_dimer_scan_C3v.xyz")

aa = AtomsList()

for h in np.linspace(0.0, 5.0, 101):
    a = Atoms()
    a.set_lattice(make_lattice(sc_size, sc_size, sc_size, angle, angle, angle))
    f = 2.2 / np.sqrt(8.0 / 3)
    d0 = f * 2.0 / 3 + 1.0
    a.add_atoms([f * np.sqrt(8.0 / 9), 0.0, f / 3], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), f * np.sqrt(2.0 / 3), f / 3], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), -f * np.sqrt(2.0 / 3), f / 3], 15)
    a.add_atoms([0.0, 0.0, -f], 15)

    a.add_atoms([f * np.sqrt(8.0 / 9), 0.0, -f / 3 + d0 + h], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), f * np.sqrt(2.0 / 3), -f / 3 + d0 + h], 15)
    a.add_atoms([-f * np.sqrt(2.0 / 9), -f * np.sqrt(2.0 / 3), -f / 3 + d0 + h], 15)
    a.add_atoms([0.0, 0.0, f + d0 + h], 15)

    a.center()
    aa.append(a)

aa.write("P4_dimer_scan_D3h.xyz")
