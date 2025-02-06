import numpy as np
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector
import matplotlib.pyplot as plt
import pickle_read_save as prs
from scipy.stats import gaussian_kde
import shapely.geometry as shape
import pandas as pd

# Source: https://matplotlib.org/stable/gallery/widgets/polygon_selector_demo.html
class SelectFromCollection:
    """
    Select indices from a matplotlib collection using `PolygonSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.poly = PolygonSelector(ax, self.onselect, draw_bounding_box=True)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.poly.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def density_scatter(x,y):
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    return z


def remove_nans_and_infs(x,y, sample_size):
    
    gate_df = pd.DataFrame()
    gate_df['x'] = x
    gate_df['y'] = y
    gate_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    gate_df.dropna(inplace=True)
    
    gate_df = gate_df.sample(sample_size)
    
    x = gate_df.x
    y = gate_df.y
    
    return x, y 


def return_markers_in_polygon(x, y, x_variable, y_variable, gate_name, sample_size, colormap):
    
    fig, ax = plt.subplots(figsize=(6,6))
    
    x, y = remove_nans_and_infs(x,y, sample_size)
    
    z = density_scatter(x,y)
    pts = ax.scatter(x, y, c=z, cmap=colormap, s=2)
    ax.set_xlabel(x_variable)
    ax.set_ylabel(y_variable)
    
    selector = SelectFromCollection(ax, pts)

    print("Select points in the figure by enclosing them within a polygon.")
    print("Press the 'esc' key to start a new polygon.")
    print("Try holding the 'shift' key to move all of the vertices.")
    print("Try holding the 'ctrl' key to move a single vertex.")

    plt.show()

    selector.disconnect()

    # After figure is closed print the coordinates of the selected points
    print('\nSelected points:')
    print(selector.xys[selector.ind])
    print(selector.ind)
    
    polygon_verts = selector.poly.verts
    # gated_data = selector.xys[selector.ind].tolist()
    # gated_data.extend(polygon_verts)  # Add the polygon vertices as a list
    gated_data = [polygon_verts, 'polygon_gate']
    prs.save_data(gated_data, gate_name)
    

def get_markers_inside_gate(x,y, polygon):
    if np.isnan(x) or np.isnan(y):
        return np.nan
    else:
        if polygon.contains(shape.Point(x,y)):
            return 1
        else:
            return 0

def get_polygon_from_coordinates(coord_list):
    return shape.Polygon(coord_list)