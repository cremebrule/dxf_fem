"""Library that contains help messages for the user to read from femgui"""

help = "\n'help' Selected:\n\n\
This is a simple FEM program for the specific plane stress \
problem shown to the left. Associated parameters can be edited above as shown. \
\n\nIf a parameter study would like to be run, please select the corresponding \
button and select the parameter(s) to be varied for the study. The files are \
automatically exported to the .vtk format which can then be opened and read in \
ParaView.\n\n \
Enter 'help' plus the following to get more detailed information \
about the specific options:\
\n\n -buttons:  Provides information on the button functions within the GUI. \
\n\n -inputs: Provides information on the user-specified inputs within the GUI. \
\n\n -param: Provides information on the parameter study functionality. \
\n\n -dxf: Provides information on formatting, importing, and using dxf files. \
\n\n -vtk: Provides information on viewing vtk files from a parameter study. \
\n\n -toolbar: Provides information on the toolbar options available at the top. \
\n\n -overview: Provides information on the scope and functionality of this program. \
\n\n"

buttons = "\n'help -buttons' Selected:\n\n \
The buttons provided in the GUI allow the user to specify the scope of the FEM study being conducted. \
The buttons' individual functions are specified as follows:\
\n\n 'Param Study is XX': Specifies whether a parameter study or normal (single) execution will be run. \
\n\n 'Geometry is XX': Option to toggle display of FEM geometry figure after calculation. \
\n\n 'Mesh is XX': Option to toggle wireframe display on outputted figures after calculation. \
\n\n 'Nodal Values is XX': Option to toggle display of model displacement after calculation. \
\n\n 'Element Values is XX': Option to toggle display of Von Mises Stress distribution after calculation. \
\n\n 'Display Model Image': Displays the current loaded DXF file geometry. \
\n\n 'Clear Screen': Clears the interface textbox screen of all text. \
\n\n"

inputs = "\n'help -inputs' Selected:\n\n \
The inputs provided in the GUI allow the user to specify various parameters related to the geometry, \
boundary conditions, and meshing properties. Enter 'help' plus the following to get more detailed information \
about the specific options:\
\n\n -geometry, -mesh, -bc\
\n\n\
"

geo = "\n'help -geometry' Selected:\n\n \
Instead of using the native DXF dimensions, these inputs can specify aspects of the model geometry prior to \
running the FEM calculation. The geometry inputs' individual functions are specified as follows:\
\n\n 'Thickness': Specifies the thickness of the plate.\
\n\n 'Dim X Name': Specifies the name of the dimension to be edited. This is the name of the corresponding layer that \
the geometry edge(s) belong to.\
\n\n '[Dim] Value': Specifies the value of the dimension to be edited. \
\n\n 'Param X Name': Functions in the same way as 'Dim X Name'.\
\n\n '[Param] start, end': If Param Study is set to ON and Vary is toggled, then these values represent the range \
between which the \
corresponding parameter will be toggled during the study. If Param Study is set to OFF, then the 'start' value \
functions the same as [Dim] Value, and 'end' value is ignored.\
\n\n 'Vary': Determines whether the specified parameter is to be toggled between the corresponding range of \
values during the parameter study. If Param Study is set to OFF, this radio button does nothing.\
\n\n\
"

mesh = "\n'help -mesh' Selected:\n\n \
These inputs alter the meshing quality as follows: \
\n\n 'Max El Size': Determines the maximum element size (in m^2) for the meshing elements. \
\n\n 'Refine': If set, increases the meshing fidelity surrounding interior geometry and holes. \
\n\n\
"

bc = "\n'help -bc' Selected:\n\n \
These inputs specify the material properties and boundary conditions of the FEM calculation, and are described \
as follows: \
\n\n 'Elastic Modulus': Determines the elastic modulus of the object. \
\n\n 'Poisson's Ratio': Determines the poisson ratio of the object. \
\n\n 'Force Name': Specifies the name of the force to be applied. This is the name of the coresponding layer that \
the geometry edge(s) belong to.*\
\n\n 'Force': Specifies the magnitude of the force applied. \
\n\n 'Angle': Specifies the direction of the force applied (in degrees). \
\n\n 'Boundary Name': Specifies the name of the boundary condition applied. This is the name of the corresponding \
layer that the geometry edge(s) belong to.*\
\n\n 'Boundary Value': Specifies the magnitude of the boundary condition (in displacement units). \
\n\n *NOTE* : The Force and Boundary Names cannot share the same name, but can be the same as \
dimension / parameter names.\
\n\n\
"

param = "\n'help -param' Selected:\n\n \
A Parameter study can be run by toggling the 'Param Study is XX' button. Up to two parameters can be varied for the \
parameter study; they are activated by toggling their respective 'Vary' radio buttons. The number of cases run is \
specified by the 'Param Step' option; this determines the number of linearly-stepped values each active parameter \
will be varied by (including the start/end points). Running the parameter study will result in an output of \
sequentially-numbered .vtk files that each represent a specific case run. The results can then be viewed in \
ParaView. For more information on viewing the .vtk files, please enter 'help -vtk'.\
\n\n\
"

dxf = "\n'help -dxf' Selected:\n\n \
The core functionality of this FEM program is centered around reading DXF files to be used for the FEM study. \
It is important that the DXF file is formatted correctly for the program to successfully read. The aspects, \
limitations, and features of the DXF-reading functionality is described below:\
\n\n Edges are specifically distinguished by the layer they belong to. Defining a new layer for a specific edge(s) \
will allow the user to specify those edge(s) as a dimension, force surface, etc.\
\n\n For readibility, it is highly recommended to include text within the same layer that displays the layer name and \
place the labels next to the corresponding edge(s) that belong to the layer. This will allow the user to easily see \
which edge(s) belong to which layer when they want to 'Display Model Image' within the GUI. \
\n\n To specify the bounding box of the object, draw a rectangle (using the rectangle tool) defining the frame \
of the object. At the top-left \
vertex of the rectangle, attach a text (note: TEXT, not MTEXT) that reads 'drawing'. Make sure both the rectangle \
frame and text belong to the layer 'svgframe'. No other geometry should belong to this layer. This will tell the FEM \
program the bounding range for the object.\
\n\n To specify the units of the object, enter 'units' and then specify the unit under the 'Insertion Scale Units' \
option.\
\n\n Currently, other than the rectangle for the frame, only the CIRCLE and LINE tools are supported for defining \
the object geometry. The CIRCLE tool can only be used for defining interior holes. \
\n\n A couple of example .dxf files have been included for reference. They are 'default.dxf' and 'frame1.dxf'\
\n\n\
"

vtk = "\n'help -vtk' Selected:\n\n \
Files outputted from a parameter study are stored in the .vtk format, and to be read by ParaView. To view the files \
in ParaView, follow the steps below: \
\n\n Under the toolbar, click 'Open' and select the .vtk cluster of the outputted files. \
\n Click the 'Apply' button on the left hand side of the GUI. \
\n Toggling through the various cases can be done by pressing the 'Rewind', 'Play', or 'Fast-Forward' buttons at \
the top of the GUI.\
\n\n\
"

toolbar = "\n'help -toolbar' Selected:\n\n \
The toolbar provides options to load and save files to the program. Below, the individual toolbar command \
functionalities are specified:\
\n\n 'New': Opens a new (default) template. This is based on 'default.dxf' file and 'default.json' template.\
\n\n 'Open': Opens a saved template (.json file). Has an internal reference to the .dxf file to be used as well.\
\n\n 'Save': Saves the current GUI inputs as a template (.json file). Has an internal reference to the .dxf file to \
be used as well. Will only ask the user to enter a name the first time this is option is used.\
\n\n 'Save As': Works the same as 'Save', but will always ask the user to enter a name for the saved template.\
\n\n 'Load DXF File': Loads a new DXF file that will be referenced as the object upon which the FEM calculation will \
be run. Note: The GUI inputs are NOT reset after this command.\
\n\n 'Execute': Executes the FEM calculation / parameter study.\
\n\n\
"

overview = "\n'help -overview' Selected:\n\n \
This specific model problem assumes the following:\
\n\n A plate object is assumed to be suspended in zero gravity.\
\n The object geometry is assumed to be constant with respect to its depth.\
\n Plane stress is assumed to be applied. \
\n Objects are assumed to remain in their elastic range (i.e.: constant proportional stress-strain relationship).\
\n\n\
"

menu = {
    'help\n' : help,
    'help -buttons\n' : buttons,
    'help -inputs\n' : inputs,
    'help -param\n' : param,
    'help -dxf\n' : dxf,
    'help -vtk\n' : vtk,
    'help -toolbar\n' : toolbar,
    'help -overview\n' : overview
}