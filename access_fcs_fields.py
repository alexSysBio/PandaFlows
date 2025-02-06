
import numpy as np

def access_single_column_data(dataframe, variable, logbool):
    
    if logbool == True:
        return np.log10(dataframe[variable]), '_log'
    else:
        return dataframe[variable], '_lin'
    

def access_double_column_data(dataframe, variables, logbools):
    
    x = dataframe[variables[0]]
    y = dataframe[variables[1]]
    
    xlog = '_lin'
    ylog = 'lin'
    
    if logbools[0] == True:
        x = np.log10(x)
        xlog = '_log'
        
    if logbools[1] == True:
        y = np.log10(y)
        ylog = 'log'
    
    return x, y, xlog+ylog
