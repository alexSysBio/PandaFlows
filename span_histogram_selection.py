
from matplotlib.widgets import PolygonSelector, SpanSelector
import numpy as np
import matplotlib.pyplot as plt
import pickle


def get_middle_of_bins(bin_edges):
    return (bin_edges[:-1]+bin_edges[1:])/2

def get_histogram_data(x, bin_array):
    return np.histogram(x, bins=bin_array)

def return_selected_ranges(x, bin_array, variable, gate_name):
    
    hist_data = get_histogram_data(x, bin_array)
    bin_edges = hist_data[1]
    x = get_middle_of_bins(bin_edges)
    y = hist_data[0]

    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6))

    ax1.plot(x, y)
    ax1.set_ylim(y.min()-0.01*y.min(), y.max()+0.01*y.max())
    ax1.set_title('Press left mouse button and drag '
                'to select a region in the top graph')
    ax1.set_xlabel(variable)

    line2, = ax2.plot([], [])


    def onselect(xmin, xmax):
        indmin, indmax = np.searchsorted(x, (xmin, xmax))
        indmax = min(len(x) - 1, indmax)

        region_x = x[indmin:indmax]
        region_y = y[indmin:indmax]

        if len(region_x) >= 2:
            line2.set_data(region_x, region_y)
            ax2.set_xlim(region_x[0], region_x[-1])
            ax2.set_ylim(region_y.min(), region_y.max())
            fig.canvas.draw_idle()
        
        print(region_x)
        with open(gate_name, 'wb') as handle:
            pickle.dump([region_x[0], region_x[-1], 'histogram_gate'], handle)

    span = SpanSelector(
        ax1,
        onselect,
        "horizontal",
        useblit=True,
        props=dict(alpha=0.5, facecolor="tab:blue"),
        interactive=True,
        drag_from_anywhere=True
    )
    # Set useblit=True on most backends for enhanced performance.
    plt.show()