# DXF FEM Calculator 1.0
A Python-based FEM Program to calculate plane-stress on DXF-based geometry.

## Installation
To install, please copy the git repo into a folder of your choosing. Then, source a Conda environment within which the program is to be run. The following packages must be included within your conda environment and can be sourced as follows:

```
pip install dxf2svg
pip install dxfgrabber
pip install ezdxf
conda install -c anaconda pyqt
pip install calfem-python
```

Any additional dependencies that cannot be found can usually be installed with `pip install {pkg_name}`. A quick google search usually clears up most errors during this stage.

Next, a couple small adjustments must be made to fix some bugs in the installed packages.

In package `dxf2svg`, within module `pycore.py`, replace any instance of:

```
get_rstrip_points()
```

with:

```
get_points()
```

Next, to solve a retina screen viewing bug in VisVis, please see the following thread to implement the quick fix:
https://github.com/almarklein/visvis/issues/97

A couple of additional software programs must be installed to run and view the output FEM data -- **Gmsh** and **ParaView**. Gmsh can be installed from http://gmsh.info/, and ParaView from https://www.paraview.org/. On Mac, the default download location is sufficient; on Windows or Linux, the file path must be specified. In module `femModel.py`, insert the following line at line 175 (after `meshGen.returnBoundaryElements = True`), where `{path_to_gmsh application executable}` is the absolute filepath location to the GMSH application executable:

```
meshGen.gmsh_exec_path='{path_to_gmsh application executable}'
```

DXF files are assumed to be created in AutoCAD; this software can be downloaded from https://www.autodesk.com/products/autocad/overview

Now you should be all set! To run the program, in your appropriately sourced conda environment and working directory that includes this git execute the following:
```
python femgui.py
```


## Known Issues
### Program crashes due to 'singular matrix' error
This is most likely due to bad geometry. You should probably check your units and provided dimensions, and make sure that you don't have any overlapping edges (e.g.: a hole that overlaps an external boundary!)

### Program seems to time out during execution stage
This is most likely due to an over-meshed setting. Try reducing (that is, increasing) the 'Max El. Size' setting and/or toggling off the 'Refine' setting before trying to execute the calculation again.

### Program crashes due to Segmentation Fault
This is an ongoing issue with an unclear cause or solution. It seems to be caused after attempting to run a calculation after a period of idleness, or after attempting a computationally-intensive study (e.g: highly-refined mesh setting). Regardless, it never occurs during the first calculation run, nor does it seem to be a deterministic error -- trying the same sequence of FEM calculations sometimes results in a Seg fault, and sometimes doesn't. Anyone with potential ideas towards a solution is more than welcome to request a commit!
