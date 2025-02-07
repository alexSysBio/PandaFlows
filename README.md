# flowio_to_pandas
This repository includes functions that can be used to transfer flow cytometry data into Pandas DataFrames and build arrays from the cell images taken by the flow cytometer

The class uses the FlowIO fcs parser:
https://github.com/whitews/FlowIO

## Initialization
The initializatio nof the class will store a flow cytometry database which can be accessed using the 
self.get_flow_cytometry_dataframe() function.

## Gating data
*Histogram gates*
<br>To apply a gate on histogram data use the self.histogram_gate() function, specifying the selected variable, if it is linear or log-transformed, and the array of the histogram bins

*Polygon gates*
<br>To apply a polygon gate on scatter-plot data, use the self.scatter_gate() function, specifying the x and y variables, if they are linear or log-transformed and the number of sampled data that will be used to draw the polygon.

*Gate operations*
<br>Reseting the gates, to remove all gate data per experiment, or importing and applying gates from a different experiment is also possible.

## Cell image operations
*Cell image segmentation*
<br>The class also includes functions to read and segment cell images, linked to the ImageFlag column in the dataframe. Useful statistics such as the cell are are also returned.
