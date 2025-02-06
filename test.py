import sys
import os
sys.path.append(os.path.join(os.path.abspath(os.getcwd()), "/flow_cytometry_class.py"))
import flow_cytometry_class as flowc
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import PolygonSelector, SpanSelector
import matplotlib.widgets as mwidgets
import span_hstogram_selection as sh

fc_data_ceph = '/Volumes/Data_01/Alex Papagiannakis/Microscopy/DRAQ5_staining/CJW7323/Flow_cytometry_CJW7323/02042025_CJW7323_M9glyCAAT_Ceph50.fcs'
fc_data_nont = '/Volumes/Data_01/Alex Papagiannakis/Microscopy/DRAQ5_staining/CJW7323/Flow_cytometry_CJW7323/02042025_CJW7323_M9glyCAAT_exp.fcs'

fn = flowc.flow_cytometry_class(fc_data_nont, 'none')
fc = flowc.flow_cytometry_class(fc_data_ceph, 'none')

# n_df = fn.get_flow_cytometry_dataframe()
# c_df = fc.get_flow_cytometry_dataframe()


# plt.figure(figsize=(7,5))
# plot_df = n_df[n_df['RL1-H']>0]
# plt.hist(plot_df['RL1-H'], bins=np.arange(0,400,5), histtype='step')
# print(plot_df['RL1-H'].median())
# plot_df = c_df[c_df['RL1-H']>0]
# plt.hist(plot_df['RL1-H'], bins=np.arange(0,400,5), histtype='step')
# print(plot_df['RL1-H'].median())
# plt.show()

# plt.figure(figsize=(7,5))
# plot_df = n_df[n_df['RL1-A']>0]
# plt.hist(plot_df['RL1-A'], bins=np.arange(0,450,7), histtype='step')
# print(plot_df['RL1-A'].median())
# plot_df = c_df[c_df['RL1-A']>0]
# plt.hist(plot_df['RL1-A'], bins=np.arange(0,450,7), histtype='step')
# print(plot_df['RL1-A'].median())
# plt.show()

# plt.figure(figsize=(7,5))
# plot_df = n_df[n_df['RL1-W']>0]
# plt.hist(plot_df['RL1-W'], bins=np.arange(0,50,1), histtype='step')
# print(plot_df['RL1-W'].median())
# plot_df = c_df[c_df['RL1-W']>0]
# hist_data = plt.hist(plot_df['RL1-W'], bins=np.arange(0,50,1), histtype='step')
# print(plot_df['RL1-W'].median())
# plt.show()

# sh.return_selected_ranges(plot_df['RL1-A'], np.arange(0,450,2), 'RL1-A')
fcs_df = fn.get_flow_cytometry_dataframe()
print(fcs_df.describe())

fn.histogram_gate('RL2-A', log=False, bin_array=np.arange(0,6000,1))
fcs_df = fn.get_flow_cytometry_dataframe()
print(fcs_df.describe())

fn.histogram_gate('RL2-A', log=False, bin_array=np.arange(0, 6000,1))
fcs_df = fn.get_flow_cytometry_dataframe()
print(fcs_df.describe())

fn.reset_gates()

fcs_df = fn.get_flow_cytometry_dataframe()
print(fcs_df.describe())

fn.histogram_gate('RL2-A', log=False, bin_array=np.arange(0, 6000,1))
fcs_df = fn.get_flow_cytometry_dataframe()
print(fcs_df.describe())