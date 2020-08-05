#!/usr/bin/env python

from analyze_utils import *
import numpy as np
from ase.units import _k
import matplotlib.pyplot as plt

(args, models, tests, default_analysis_settings) = analyze_start('bulk_Si_diamond')
bulk_data = read_properties(models, tests, args.test_set)

(args, models, tests, default_analysis_settings) = analyze_start('surface_Si_diamond_*')
surface_data = read_properties(models, tests, args.test_set)

(args, models, tests, default_analysis_settings) = analyze_start('point_defect_Si_diamond_*')
point_defect_data = read_properties(models, tests, args.test_set)

ref_linestyles=[ "-", "--" ]
other_linestyles=[ ":", "-." ]
struct_colors = [ "black", "red", "blue", "orange", "green", "brown", "grey", "magenta","cyan" ]
ref_model_name = default_analysis_settings["ref_model"]

all_data = {}

all_data = {}

for model in models:
    all_data[model] = {}

    all_data[model]["B"] = bulk_data[model]["bulk_Si_diamond"]["B"]
    all_data[model]["c11"] = bulk_data[model]["bulk_Si_diamond"]["c11"]
    all_data[model]["c12"] = bulk_data[model]["bulk_Si_diamond"]["c12"]
    all_data[model]["c44"] = bulk_data[model]["bulk_Si_diamond"]["c44"]

    all_data[model]["surf_E_100"] = surface_data[model]["surface_Si_diamond_100"]["Ef"]
    all_data[model]["surf_E_110"] = surface_data[model]["surface_Si_diamond_110"]["Ef"]
    all_data[model]["surf_E_111"] = surface_data[model]["surface_Si_diamond_111"]["Ef"]

    all_data[model]["point_defect_Si_diamond_interstitial_tetr"] = point_defect_data[model]["point_defect_Si_diamond_interstitial_tetr"]["defects"]["Z_14"]["Ef"]
    all_data[model]["point_defect_Si_diamond_vacancy"] = point_defect_data[model]["point_defect_Si_diamond_vacancy"]["defects"]["ind_0_Z_14"]["Ef"]
    all_data[model]["point_defect_Si_diamond_interstitial_hex"] = point_defect_data[model]["point_defect_Si_diamond_interstitial_hex"]["defects"]["Z_14"]["Ef"]

latex_dict = {
    'B' : '$B$',
  'c11' :'$C_{11}$',
  'c12' :'$C_{12}$',
  'c44' :'$C_{44}$',
 'point_defect_Si_diamond_vacancy' : '$\\mathrm{vac}$',
  'point_defect_Si_diamond_interstitial_hex' : '$\\mathrm{hex.\\ int.}$',
  'point_defect_Si_diamond_interstitial_tetr' : '$\\mathrm{tetr.\\ int.}$',
   'surf_E_111' : '$(111)$',
    'surf_E_110': '$(110)$',
    'surf_E_100': '$(100)$'
}

test_names = []

for model in models:
    if model != ref_model_name:
        for obs in all_data[model].keys():
            if obs != "surf_E_110": ####
                test_names.append(latex_dict[obs])
                all_data[model][obs] = ((all_data[model][obs] - all_data[ref_model_name][obs])/all_data[ref_model_name][obs]) * 100

print(all_data)

fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(111, projection='3d')

_x = np.arange(1 )
_y = np.arange(len(all_data[ref_model_name].keys())-1) ###

_xx, _yy = np.meshgrid(_x, _y)
y, x = _xx.ravel(), _yy.ravel()

width = 0.5
depth = 0.3

model_vals = []
model_ticks = []

models.remove(ref_model_name)

for (i,model) in enumerate(models):
    #if model != ref_model_name:
    plot_d = [(key, value) for (key,value) in all_data[model].items() if key != "surf_E_110"] ####
    top = [abs(d[1]) for d in plot_d]
    bottom = np.zeros_like(top)

    ax1.bar3d(x, y, bottom, width, depth, top, shade=True)

    model_vals.append(i*width)
    model_ticks.append(model)

    y = [_y + width for _y in y]

ax1.set_ylim(0,((width+depth)/2.0)*(len(models)))

ax1.set_zlabel("Percentage Error [%]")

print(test_names)

plt.xticks([i for i in range(len(all_data[ref_model_name].items())-1)], test_names, rotation=30, ha='right')
plt.yticks(model_vals, model_ticks, ha='left', va='center')

plt.savefig("bar_plot.pdf")

# fig = plt.figure(figsize=(8, 3))
# for model in models:
#     if model != ref_model_name:
#         x.append()


# print(bulk_data)
# print(surface_data)
# print(point_defect_data)
