#!/usr/bin/env python

import numpy as np
from analyze_utils import *
import matplotlib

matplotlib.use("PDF")
from matplotlib import pyplot

import phonopy
import ase.units

(args, models, tests, default_analysis_settings) = analyze_start(["phonon_bulks"])

data = read_properties(models, tests, args.test_set)

ref_model_name = default_analysis_settings["ref_model"]


def analyze_phonons(m_data):
    THz_per_invcm = ase.units._c * 1.0e-10

    at0 = ase.atoms.Atoms(
        symbols=m_data["symb"],
        scaled_positions=np.asarray(m_data["scaled_pos"]),
        masses=np.asarray(m_data["m"]),
        cell=np.asarray(m_data["c"]),
        pbc=[True] * 3,
    )
    phonopy_atoms = phonopy.structure.atoms.PhonopyAtoms(
        symbols=m_data["symb"],
        scaled_positions=np.asarray(m_data["scaled_pos"]),
        masses=np.asarray(m_data["m"]),
        cell=np.asarray(m_data["c"]),
    )
    phonons = phonopy.Phonopy(
        phonopy_atoms, m_data["n_cell"], factor=m_data["unit_factor"]
    )
    phonons.generate_displacements(distance=m_data["dx"])

    phonons.produce_force_constants(
        forces=np.asarray(m_data["all_forces"]), calculate_full_force_constants=True
    )
    phonons.symmetrize_force_constants()

    # DOS
    n_dos_mesh = 16
    phonons.run_mesh([n_dos_mesh] * 3)

    phonons.run_total_dos()
    frequencies = phonons.get_total_dos_dict()["frequency_points"]
    PHdos = phonons.get_total_dos_dict()["total_dos"]

    # contains by columns the frequencies in cm^{-1} and the vDOS
    # the vDOS ins in units of "number of states/(unit cell x frequency[cm^{-1}])"
    # i.e. if you integrate the vDOS throughout frequency, it will be 3N, where N is the number of atoms in the unit cell
    m_data["DOS"] = {"freq": frequencies / THz_per_invcm, "val": PHdos * THz_per_invcm}

    # band path
    if "band_path" in m_data:
        band_path = m_data["band_path"]
        band_n = [20]

        lat = at0.get_cell().get_bravais_lattice()
        special_points = lat.get_special_points()
        path_pts = []
        path_labels = []
        pt = []
        for s in " ".join(band_path).split():
            try:
                qi = float(s)
                pt.append(qi)
            except:
                if len(pt) > 0:
                    raise RuntimeError("got non-float after 1 or more floats")
                for p in s:
                    try:
                        path_pts.append(np.asarray(special_points[p]))
                        path_labels.append(p)
                    except KeyError:
                        raise RuntimeError(
                            "Failed to find special point {}, known {}".format(
                                p, special_points.keys()
                            )
                        )
            if len(pt) == 3:
                path_pts.append(np.asarray(pt))
                path_labels.append("{:.2f}_{:.2f}_{:.2f}".format(*pt))
                pt = []

        if len(band_n) == 1:
            band_n *= len(path_pts) - 1
        assert len(band_n) == len(path_pts) - 1

        q_pts = []
        q_pt_labels = []
        for (seg_i, (label, n, p0, p1)) in enumerate(
            zip(path_labels[0:-1], band_n, path_pts[0:-1], path_pts[1:])
        ):
            for p_i in range(n):
                x = float(p_i) / float(n)

                q_pts.append((1.0 - x) * p0 + x * p1)

                if p_i == 0:
                    q_pt_labels.append(label)
                else:
                    q_pt_labels.append(".")

        q_pts.append(path_pts[-1])
        q_pt_labels.append(path_labels[-1])

        phonons.run_band_structure([q_pts])
        bs = phonons.get_band_structure_dict()
        m_data["BAND_PATH"] = {
            "positions": bs["distances"][0],
            "frequencies": bs["frequencies"][0] / THz_per_invcm,
            "labels": q_pt_labels,
        }


if ref_model_name in data and "phonon_bulks" in data[ref_model_name]:
    for (bulk_i, bulk_struct_test) in enumerate(data[ref_model_name]["phonon_bulks"]):
        print("analyze ref model", ref_model_name, bulk_struct_test)
        analyze_phonons(data[ref_model_name]["phonon_bulks"][bulk_struct_test])

bulk_struct_tests = []
for model_name in models:
    if model_name in data and "phonon_bulks" in data[model_name]:
        for bulk_struct_test in data[model_name]["phonon_bulks"]:
            if bulk_struct_test not in bulk_struct_tests:
                bulk_struct_tests.append(bulk_struct_test)

fig_DOS = pyplot.figure()
ax_DOS = fig_DOS.add_subplot(1, 1, 1)
fig_BP = pyplot.figure()
ax_BP = fig_BP.add_subplot(1, 1, 1)

# label_dict = Dict(
#
# )
colors = ["red", "green", "blue", "cyan"]

k = 0

for (j, model_name) in enumerate(models):
    if model_name == ref_model_name and len(models) > 1:
        continue
    print("analyze model", model_name)
    for (bulk_i, bulk_struct_test) in enumerate(bulk_struct_tests):
        if "phonon_bulks" not in data[model_name]:
            continue

        try:
            ref_model_data = data[ref_model_name]["phonon_bulks"][bulk_struct_test]
        except:
            ref_model_data = None
        try:
            model_data = data[model_name]["phonon_bulks"][bulk_struct_test]
        except:
            continue

        print("analyze model-bulk", model_name, bulk_struct_test)
        analyze_phonons(model_data)

        if ref_model_data is not None and "BAND_PATH" in ref_model_data:
            for i in range(ref_model_data["BAND_PATH"]["frequencies"].shape[1]):
                ax_BP.plot(
                    ref_model_data["BAND_PATH"]["positions"],
                    ref_model_data["BAND_PATH"]["frequencies"][:, i],
                    "-",
                    color="black",
                    label=ref_model_name if k == 0 else None,
                )
                k += 1

        if "BAND_PATH" in model_data:
            for i in range(model_data["BAND_PATH"]["frequencies"].shape[1]):
                ax_BP.plot(
                    model_data["BAND_PATH"]["positions"],
                    model_data["BAND_PATH"]["frequencies"][:, i],
                    "-",
                    color=colors[j],
                    label=model_name if i == 0 else None,
                )

            ax_BP.set_xticks(
                [
                    p
                    for p, l in zip(
                        model_data["BAND_PATH"]["positions"],
                        model_data["BAND_PATH"]["labels"],
                    )
                    if l != "."
                ]
            )
            ax_BP.set_xticklabels(
                [l for l in model_data["BAND_PATH"]["labels"] if l != "."]
            )
            ax_BP.set_xlim(
                model_data["BAND_PATH"]["positions"][0],
                model_data["BAND_PATH"]["positions"][-1],
            )

            # ax_BP.set_ylim(ylim)
            ylim = ax_DOS.get_ylim()
            if bulk_i == len(bulk_struct_tests) - 1:
                ylim = ax_BP.get_ylim()
                for p, l in zip(
                    ref_model_data["BAND_PATH"]["positions"][1:-1],
                    ref_model_data["BAND_PATH"]["labels"][1:-1],
                ):
                    if l != ".":
                        ax_BP.plot(
                            [p, p], ylim, "-", color="black", label=None, linewidth=0.5
                        )

ax_BP.set_xlabel("BZ point")
ax_BP.set_ylim(ylim)
ax_BP.legend(loc="upper right")
ax_BP.axhline(color="black", linewidth=0.5)
ax_BP.set_ylabel("freq. (cm$^{-1}$)")
ax_DOS.set_xlabel("freq (cm$^{-1}$)")
ax_DOS.set_ylabel("DOS (arb. units)")
ylim = ax_DOS.get_ylim()
ax_DOS.set_ylim(0.0, ylim[1])
ax_DOS.legend()

fig_DOS.savefig("phonon_bulks_DOS-PAPER.pdf")
fig_BP.savefig("phonon_bulks_BAND_PATH-PAPER.pdf")
pyplot.clf()
