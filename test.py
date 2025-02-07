import sys
import os
sys.path.append(os.path.join(os.path.abspath(os.getcwd()), "/flow_cytometry_class.py"))
import flow_cytometry_class as flowc
import pickle_read_save as prs
import numpy as np

fc_data_ceph = '.fcs'
fc = flowc.flow_cytometry_class(fc_data_ceph, 'none', 'experiment_id')
# get the Pandas dataframe
fcs_df = fc.get_flow_cytometry_dataframe()
print(fcs_df.describe())
# apply histogram gate with log-transformed data
fc.histogram_gate('RL2-A', log=True, bin_array=np.arange(0, 5, 0.001))
fcs_df = fc.get_flow_cytometry_dataframe()
print(fcs_df.describe())
# remove all gates
fc.reset_gates()
# apply polygon gate to log-transgormed scatter data
x_variable = 'SSC-A'
y_variable = 'RL2-A'
fc.scatter_gate(x_variable, y_variable, True, True, 10000)
stored_gates = prs.load_data('stored_gates_02042025_CJW7323_M9glyCAAT_Ceph50')
# reload Pandas dataframe and describe the stored data and gates
fcs_df = fc.get_flow_cytometry_dataframe()
print(fcs_df.describe())
