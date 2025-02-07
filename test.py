import sys
import os
sys.path.append(os.path.join(os.path.abspath(os.getcwd()), "/flow_cytometry_class.py"))
import flow_cytometry_class as flowc
import pickle_read_save as prs
import numpy as np

fc_data_ceph = '/Users/alexandros/Downloads/02042025_CJW7323_M9glyCAAT_Ceph50.fcs'
# fc_data_nont = '/Users/alexandros/Downloads/02042025_CJW7323_M9glyCAAT_exp.fcs'

# fn = flowc.flow_cytometry_class(fc_data_nont, 'none')
fc = flowc.flow_cytometry_class(fc_data_ceph, 'none', '02042025_CJW7323_M9glyCAAT_Ceph50')

# sh.return_selected_ranges(plot_df['RL1-A'], np.arange(0,450,2), 'RL1-A')
# fcs_df = fn.get_flow_cytometry_dataframe()
# print(fcs_df.describe())

# fn.histogram_gate('RL2-A', log=False, bin_array=np.arange(0,6000,1))
fcs_df = fc.get_flow_cytometry_dataframe()
print(fcs_df.describe())

fc.histogram_gate('RL2-A', log=True, bin_array=np.arange(0, 5, 0.001))
fcs_df = fc.get_flow_cytometry_dataframe()
print(fcs_df.describe())

fc.reset_gates()

x_variable = 'SSC-A'
y_variable = 'RL2-A'
fc.scatter_gate(x_variable, y_variable, True, True, 10000)
stored_gates = prs.load_data('stored_gates_02042025_CJW7323_M9glyCAAT_Ceph50')

fcs_df = fc.get_flow_cytometry_dataframe()
print(fcs_df.describe())