"""Classes that will be used to define the various suborutines necessary to perform FEM calculations

The following classes will be used:

InputData: Stores the input variables required for calculations
Solution: Impelements routine solution for the problem
OutputData: Stores the results from the Solution object
Report: Manages the I/O reports for the program"""

import numpy as np
import math
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import json
import pyvtk as vtk
import femDXF

import visvis as vv

global globalVisVisApp

global visApp
visApp = vv.use('qt5') # use qt4

#### Define global constants ####
# Units to SI units: [m, N/m, Pa], from Imperial [in, lbf/in, psi]
U2SI = {
    "SI" : [1.0, 1.0, 1.0e9],
    "IMPERIAL" : [0.0254, 175.127, 6894.76]
}


class InputData(object):
    """Class for defining all the input data for the model"""
    def __init__(self):

        self.version = 1

        # ------ Element Properties
        self.E = None  # Elastic Modulus, psi
        self.t = None  # thickness, in
        self.d = None  # Dimensions [[dim1 name, dim1 val], [dim2 name, dim2 val], ...]
        self.c = None  # additional constants, in (dictionary)
        self.refineMesh = None  # Determine whether to refine mesh
        self.v = None  # Poisson's Ratio, unitless
        self.ptype = None   # problem type (1 = plane stress)
        self.ep = None      # element properties [ptype, t]
        self.mp = None      # mesh parameters [elType, dofsPerNode, elSizeFactor]
        self.bp = None      # boundary condition parameters [[marker(s)], [values]]
        self.fp = None      # force parameters [[marker(s)], [values], [angles]]
        self.dxf = femDXF.InputDXF()
        self.dxf_filename = None
        self.units = None   # Units for the system (SI, or IMPERIAL)

        # ------ Extra Properties for Parameter Study
        self.paramFilename = None       # Filename for parameter study
        self.paramSteps = None          # Number of steps for parameter study

    def updateparams(self):
        """Updates internal parameters that depend on other internal parameters"""
        self.ep = [self.ptype, self.t]

    def save(self, filename):
        # ------ Used for saving the data to a file

        inputData = {}
        inputData["version"] = self.version
        inputData["E"] = self.E
        inputData["t"] = self.t
        inputData["d"] = self.d
        inputData["v"] = self.v
        inputData["c"] = self.c
        inputData["ptype"] = self.ptype
        inputData["mp"] = self.mp
        inputData["fp"] = self.fp
        inputData["bp"] = self.bp
        inputData["paramSteps"] = self.paramSteps
        inputData["dxf_filename"] = self.dxf_filename
        inputData["units"] = self.units

        # Open and write to the file
        ofile = open(filename, "w")
        json.dump(inputData, ofile, sort_keys=True, indent=4)
        ofile.close()

    def load(self, filename):
        # ------ Used for loading the data from a file

        ifile = open(filename, "r")
        inputData = json.load(ifile)
        ifile.close()

        # Convert the data to a usable format
        self.version = inputData["version"]
        self.E = inputData["E"]
        self.t = inputData["t"]
        self.d = inputData["d"]
        self.v = inputData["v"]
        self.c = inputData["c"]
        self.ptype = inputData["ptype"]   # problem type (1 = plane stress)
        self.mp = inputData["mp"]
        self.fp = inputData["fp"]
        self.bp = inputData["bp"]
        self.paramSteps = inputData["paramSteps"]
        self.dxf_filename = inputData["dxf_filename"]
        self.dxf.readDXF(self.dxf_filename)
        self.dxf.convertDXFtoSVG()
        self.units = inputData["units"]
        self.updateparams()


class OutputData(object):
    """Class for storing the results from the calculation"""
    def __init__(self):
        self.disp = None
        self.stress = None
        self.geometry = None
        self.a = None
        self.coords = None
        self.edof = None
        self.mp = None
        self.meshGen = None
        self.paramnum = 0
        self.statistics = None      # [max VMstress, max disp, curveID(s) of layers, anchor location, w/h of drawing]


class Solver(object):
    """Class to handle the solver algorithm of our solution model"""
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData

    def execute(self):
        # ------ Transfer model variables to local variables
        self.inputData.updateparams()
        version = self.inputData.version
        units = self.inputData.units
        v = self.inputData.v
        ep = self.inputData.ep
        E = self.inputData.E
        mp = self.inputData.mp
        fp = self.inputData.fp
        bp = self.inputData.bp
        ep[1] = ep[1] * U2SI[units][0]
        E = E * U2SI[units][2]
        for i in range(len(fp[0])):
            fp[1][i] = fp[1][i] * U2SI[units][1]
        for i in range(len(bp[0])):
            bp[1][i] = bp[1][i] * U2SI[units][0]

        # Get most updated dxf dimensions and import model geometry to calfem format
        self.inputData.dxf.readDXF(self.inputData.dxf_filename)
        for dim in self.inputData.d:
            ("Adjusting Dimension {0} with val {1}".format(dim[0], dim[1] * U2SI[units][0]))
            self.inputData.dxf.adjustDimension(dim[0], dim[1] * U2SI[units][0])
        self.inputData.dxf.adjustDimension(self.inputData.c['aName'], self.inputData.c['a'] * U2SI[units][0])
        self.inputData.dxf.adjustDimension(self.inputData.c['bName'], self.inputData.c['b'] * U2SI[units][0])
        dxf = self.inputData.dxf

        if self.inputData.refineMesh:
            geometry, curve_dict = dxf.convertToGeometry(max_el_size=mp[2])
        else:
            geometry, curve_dict = dxf.convertToGeometry()

        # Generate the mesh
        meshGen = cfm.GmshMeshGenerator(geometry)
        meshGen.elSizeFactor = mp[2]                  # Max Area for elements
        meshGen.elType = mp[0]
        meshGen.dofsPerNode = mp[1]
        meshGen.returnBoundaryElements = True

        coords, edof, dofs, bdofs, elementmarkers, boundaryElements = meshGen.create()

        # Add the force loads and boundary conditions
        bc = np.array([], int)
        bcVal = np.array([], int)
        nDofs = np.size(dofs)
        f = np.zeros([nDofs, 1])

        for i in range(len(bp[0])):
            bc, bcVal = cfu.applybc(bdofs, bc, bcVal, dxf.markers[bp[0][i]], bp[1][i])
        for i in range(len(fp[0])):
            xforce = fp[1][i] * np.cos(np.radians(fp[2][i]))
            yforce = fp[1][i] * np.sin(np.radians(fp[2][i]))
            cfu.applyforce(bdofs, f, dxf.markers[fp[0][i]], xforce, dimension=1)
            cfu.applyforce(bdofs, f, dxf.markers[fp[0][i]], yforce, dimension=2)

        # ------ Calculate the solution

        print("")
        print("Solving the equation system...")

        # Define the elements coordinates
        ex, ey = cfc.coordxtr(edof, coords, dofs)

        # Define the D and K matrices
        D = (E / (1 - v**2))*np.matrix([
            [1, v, 0],
            [v, 1, 0],
            [0, 0, (1-v)/2]
             ])
        K = np.zeros([nDofs, nDofs])

        # Extract element coordinates and topology for each element
        for eltopo, elx, ely in zip(edof, ex, ey):
            Ke = cfc.plante(elx, ely, ep, D)
            cfc.assem(eltopo, K, Ke)

        # Solve the system
        a, r = cfc.solveq(K, f, bc, bcVal)

        # ------ Determine stresses and displacements

        print("Computing the element forces")

        # Extract element displacements
        ed = cfc.extractEldisp(edof, a)

        # Determine max displacement
        max_disp = [[0, 0], 0]       # [node idx, value]
        for idx, node in zip(range(len(ed)), ed):
            for i in range(3):
                disp = math.sqrt(node[2*i]**2 + node[2*i+1]**2)
                if disp > max_disp[1]:
                    max_disp = [[idx, 2*i], disp]

        # Determine Von Mises stresses
        vonMises = []
        max_vm = [0, 0]  # [node idx, value]
        for i in range(edof.shape[0]):
            es, et = cfc.plants(ex[i, :], ey[i, :], ep, D, ed[i, :])
            try:
                vonMises.append(math.sqrt(pow(es[0, 0], 2) - es[0, 0] * es[0, 1] + pow(es[0, 1], 2) + 3 * es[0, 2]))
                if vonMises[-1] > max_vm[1]:
                    max_vm = [i, vonMises[-1]]
            except ValueError:
                vonMises.append(0)
                print("CAUGHT MATH EXCEPTION with es = {0}".format(es))
            # Note: es = [sigx sigy tauxy]

        # ------ Store the solution in the output model variables
        self.outputData.disp = ed
        self.outputData.stress = vonMises
        self.outputData.geometry = geometry
        self.outputData.a = a
        self.outputData.coords = coords
        self.outputData.edof = edof
        self.outputData.mp = mp
        self.outputData.meshGen = meshGen
        self.outputData.statistics = [max_vm, max_disp, curve_dict, self.inputData.dxf.anchor, self.inputData.dxf.wh]
        if self.inputData.paramFilename is None:
            print("Solution completed.")

    def executeParamStudy(self):
        """Method that runs a parameter study; still utilizes execute() as main calculation method"""

        # --- Create vectors of varying values for each respective parameter that was selected for the study
        if self.inputData.c['paramA']:
            aRange = np.linspace(self.inputData.c['aStart'], self.inputData.c['aEnd'], self.inputData.paramSteps)
        else:
            aRange = [self.inputData.c['a']]
        if self.inputData.c['paramB']:
            bRange = np.linspace(self.inputData.c['bStart'], self.inputData.c['bEnd'], self.inputData.paramSteps)
        else:
            bRange = [self.inputData.c['b']]

        # --- Run the parameter study for each combination of parameters
        i = 0
        for a in aRange:
            for b in bRange:
                print("\nExecuting combination a: {0} and b: {1}" .format(a, b))
                i = i + 1
                self.inputData.c['a'] = float(a)
                self.inputData.c['b'] = float(b)
                if i < 10:
                    num_str = "0" + str(i)
                else:
                    num_str = str(i)
                vtk_filename = self.inputData.paramFilename + "_" + num_str
                self.execute()
                self.exportVtk(vtk_filename)
        print("Successfully completed %i studies." % i)
        self.outputData.paramnum = i

    def exportVtk(self, filename):
        """Method for exporting fem calculation output to VTK-compatible format"""
        print("Exporting results to '%s'..." % filename)

        # --- Create points and polygon definitions from our node network
        points = self.outputData.coords.tolist()

        # --- Make sure topology is VTK-compatible; i.e.: 0-based
        #polygons = (self.outputData.edof-1).tolist()
        topo = np.zeros([self.outputData.edof.shape[0], 3], dtype=int)
        for i in range(self.outputData.edof.shape[0]):
            topo[i, 0] = self.outputData.edof[i,1]/2 - 1
            topo[i, 1] = self.outputData.edof[i, 3] / 2 - 1
            topo[i, 2] = self.outputData.edof[i, 5] / 2 - 1

        polygons = (topo).tolist()

        # --- Specify both vector and scalar data for each element
        #pointData = vtk.PointData(vtk.Scalars(self.outputData.a.tolist(), name="Displacement"))
        #cellData = vtk.CellData(vtk.Scalars(max(self.outputData.stress), name="maxvmstress"),\
        #                        vtk.Vectors(self.outputData.stress, "stress"))
        cellData = vtk.CellData(vtk.Scalars(self.outputData.stress, name="Von Mises"))

        # --- Create the structure of the element network
        structure = vtk.PolyData(points=points, polygons=polygons)

        # --- Store everything in a vtk instance
        #vtkData = vtk.VtkData(structure, pointData, cellData)
        vtkData = vtk.VtkData(structure, cellData)

        # --- Save the data to the specified file
        vtkData.tofile(filename, "ascii")

class Visualization(object):
    """Class for visualizing the results from the Solver"""
    def __init__(self, inputData, outputData, figsOn):
        self.inputData = inputData
        self.outputData = outputData

        # --- Variables for references to the various gmsh figures
        self.geomFig = figsOn[0]
        self.meshFig = figsOn[1]
        self.nodeValueFig = figsOn[2]
        self.elValueFig = figsOn[3]

    def show(self):
        # ------ Shows the geometry
        geometry = self.outputData.geometry
        a = self.outputData.a
        vonMises = self.outputData.stress
        coords = self.outputData.coords
        edof = self.outputData.edof
        dofsPerNode = self.outputData.mp[1]
        elType = self.outputData.mp[0]
        meshGen = self.outputData.meshGen
        stats = self.outputData.statistics
        fp = self.inputData.fp
        bp = self.inputData.bp
        w = stats[4][0]
        h = stats[4][1]
        units = self.inputData.units

        # Create the figure
        print("Visualizing...")
        cfv.close_all()
        if(self.geomFig):
            cfv.figure()
            cfv.drawGeometry(geometry, title="Geometry", drawPoints=False, labelCurves=True)

        if(self.elValueFig):
            cfv.figure()
            cfv.drawElementValues(vonMises, coords, edof, dofsPerNode, elType, a, doDrawMesh=self.meshFig,
                              doDrawUndisplacedMesh=False, title="Effective (Von Mises) Stress (Pa)")

            # ------ Add extra text
            node_x = [coords[int((edof[stats[0][0]][0]-1)/2), 0] + a[edof[stats[0][0]][0]-1, 0],
                      coords[int((edof[stats[0][0]][2]-1)/2), 0] + a[edof[stats[0][0]][2]-1, 0],
                      coords[int((edof[stats[0][0]][4]-1)/2), 0] + a[edof[stats[0][0]][4]-1, 0]]
            node_y = [coords[int((edof[stats[0][0]][0]-1)/2), 1] + a[edof[stats[0][0]][1]-1, 0],
                      coords[int((edof[stats[0][0]][2]-1)/2), 1] + a[edof[stats[0][0]][3]-1, 0],
                      coords[int((edof[stats[0][0]][4]-1)/2), 1] + a[edof[stats[0][0]][5]-1, 0]]
            node_centroid = [sum(node_x) / 3, sum(node_y) / 3]
            cfv.addText("Max Stress: {0:.2f} MPa".format(
                stats[0][1]/1e6), (stats[3][0] + 0.725*w, stats[3][1] + 0.925*h), fontSize=15, color='w')
            cfv.vv.plot([stats[3][0] + 0.715*w, node_centroid[0]], [stats[3][1] + 0.925*h, node_centroid[1]], lc='w')

        if(self.nodeValueFig):
            cfv.figure()
            cfv.drawDisplacements(a, coords, edof, dofsPerNode, elType, doDrawUndisplacedMesh=True, title="Displacements (m)")

            # Add markers
            symbols = {
                "rightarrow"    : "\u2192",
                "leftarrow"     : "\u2190",
                "uparrow"       : "\u2191",
                "downarrow"     : "\u2193",
                "nearrow"       : "\u2197",
                "nwarrow"       : "\u2196",
                "swarrow"       : "\u2199",
                "searrow"       : "\u2198",
                "fixed"         : "\u2215"
            }

            forceNodes = set()
            bcNodes = set()
            for j in range(len(fp[0])):
                for curveID in stats[2][fp[0][j]]:
                    forceNodes = forceNodes.union(set(meshGen.nodesOnCurve[curveID]))
                for i in forceNodes:
                    x = coords[i, 0] + a[i * 2, 0]    # Position of node with displacements
                    y = coords[i, 1] + a[i * 2 + 1, 0]
                    cfv.addText(symbols["rightarrow"], (x, y), angle=fp[2][j], fontSize=20, color='g')
                for curveID in stats[2][bp[0][j]]:
                    bcNodes = bcNodes.union(set(meshGen.nodesOnCurve[curveID]))
                for i in bcNodes:
                    x = coords[i, 0] + a[i * 2, 0]  # Position of node with displacements
                    y = coords[i, 1] + a[i * 2 + 1, 0]
                    cfv.addText(symbols["fixed"], (x, y), fontSize=15, color='r')

            # --- Add additional text
            cfv.addText("Forces Applied: {0:6.2f} kN/m".format(self.inputData.fp[1][0]*U2SI[units][1]/1e3),
                        (stats[3][0] + 0.725*w, stats[3][1] + 0.925*h), fontSize=15, color='g')
            cfv.addText("Boundary Condition: {0:6.2f}m displacement".format(
                self.inputData.bp[1][0]*U2SI[units][0]),
                (stats[3][0] + 0.625*w, stats[3][1] + 0.95*h), fontSize=15, color='r')
            node_x = coords[int((edof[stats[1][0][0]][stats[1][0][1]] - 1) / 2), 0] +\
                     + a[edof[stats[1][0][0]][stats[1][0][1]] - 1, 0]
            node_y = coords[int((edof[stats[1][0][0]][stats[1][0][1]] - 1) / 2), 1] +\
                     + a[edof[stats[1][0][0]][stats[1][0][1]+1] - 1, 0]
            cfv.addText("Max Displacement: {0:6.2f} mm".format(
                stats[1][1]*1e3), (stats[3][0] + 0.725*w, stats[3][1] + 0.9*h), fontSize=15, color='w')
            cfv.vv.plot([stats[3][0] + 0.715*w, node_x], [stats[3][1] + 0.9*h, node_y], lc='w')

    def closeAll(self):
        # ------ Closes all windows and resets the character variables
        cfv.close_all()
        self.geomFig = 0
        self.meshFig = 0
        self.elValueFig = 0
        self.nodeValueFig = 0

    def wait(self):
        # ------ Ensures that the windows are kept updated and will return when the last window is closed
        cfv.showAndWait()


class Report(object):
    """Class for the report of input and output data in report form"""
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData
        self.report = ""

    def clear(self):
        self.report = ""

    def addText(self, text=""):
        self.report += str(text) + "\n"

    def __str__(self):
        # --- String to output to user
        units = self.inputData.units
        self.clear()
        self.addText()
        self.addText("---------------- MODEL INPUT -----------------")
        self.addText()
        self.addText("// Material Properities //")
        self.addText()
        self.addText("Modulus of Elasticity: " + str(self.inputData.E * U2SI[units][2] * 1.0e-9) + " GPa")
        self.addText("Poisson's Ratio: " + str(self.inputData.v))
        self.addText("Plate Thickness: " + str(self.inputData.t * U2SI[units][0]) + " m")
        for d in self.inputData.d:
            if d[1] > 0 and d[0] in self.inputData.dxf.markers.keys():
                self.addText("Plate Dimension '{0}': {1} m".format(d[0], d[1] * U2SI[units][0]))
        if self.inputData.c['a'] > 0 and self.inputData.c['aName'] in self.inputData.dxf.markers.keys():
            self.addText("Plate Dimension '{0}': {1} m".format(self.inputData.c['aName'],
                                                               self.inputData.c['a'] * U2SI[units][0]))
        if self.inputData.c['b'] > 0 and self.inputData.c['bName'] in self.inputData.dxf.markers.keys():
            self.addText("Plate Dimension '{0}': {1} m".format(self.inputData.c['bName'],
                                                               self.inputData.c['b'] * U2SI[units][0]))
        self.addText()
        self.addText("// Mesh Properties //")
        self.addText()
        self.addText("Mesh Parameters [elType, dofsPerNode, elSizeFactor]:")
        self.addText()
        self.addText(self.inputData.mp)
        self.addText()
        self.addText("Boundary Condition Parameters [[marker(s)], [values (m)]]:")
        self.addText()
        self.addText(self.inputData.bp)
        self.addText()
        self.addText("Applied Load Parameters [[marker(s)], [values (N/m)], [angles]]:")
        self.addText()
        self.addText(self.inputData.fp)
        self.addText()
        self.addText()
        self.addText("---------------- MODEL OUTPUT ----------------")
        self.addText()
        self.addText("Max Displacement [[element, vertex], disp] (m): ")
        self.addText()
        self.addText("{0}".format(self.outputData.statistics[1]))
        self.addText()
        self.addText("Max Element Von Mises Stress [Element, stress] (Pa): ")
        self.addText()
        self.addText("{0}".format(self.outputData.statistics[0]))

        # -- More elegant to not include tons of data -- can re-implement if necessary in the future
        """
        self.addText("Displacements [DOF, disp] (m): ")
        self.addText()
        for element in self.outputData.disp:
            self.addText(element)
        #self.addText(self.outputData.disp)
        self.addText()
        self.addText("Element Stresses (Pa): ")
        self.addText()
        for element in self.outputData.stress:
            self.addText(element)
        #self.addText(self.outputData.stress)
        """
        self.addText()
        self.addText()
        return self.report
