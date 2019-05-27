"""Classes that will be used to read and parse DXF files to then be read by the femgui.py and femModel.py

The following classes will be used:

InputDXF: Reads a DXF file and returns/manipulates specific objects from the file

"""


import dxfgrabber
from dxf2svg.pycore import save_svg_from_dxf, extract_all
import calfem.geometry as cfg
import numpy as np

# ------ Define global constants
# ACAD Version Dictionary
ACVERSIONDICT = {"AC1009" : "AutoCAD R12", "AC1015" : "AutoCAD R2000",\
                 "AC1018" : "AutoCAD R2004", "AC1021" : "AutoCAD R2007",\
                 "AC1024" : "AutoCAD R2010", "AC1027" : "AutoCAD R2013",\
                 "AC1032" : "AutoCAD R2018"}

# ACAD Unit Dictionary
ACUNITS = {
    0 : 'UNITLESS',
    1 : 'INCHES',
    2 : 'FEET',
    3 : 'MILES',
    4 : 'MILLIMETERS',
    5 : 'CENTIMETERS',
    6 : 'METERS',
    7 : 'KILOMETERS',
    8 : 'MICROINCHES',
    9 : 'MILS',
    10 : 'YARDS',
    11 : 'ANGSTROMS',
    12 : 'NANOMETERS',
    13 : 'MICRONS',
    14 : 'DECIMETERS',
    15 : 'DECAMETERS',
    16 : 'HECTOMETERS',
    17 : 'GIGAMETERS',
    18 : 'ASTRONOMICAL UNITS',
    19 : 'LIGHT YEARS',
    20 : 'PARSECS'
}

# Unit conversion Dictionary from (Unit -> Meters)
ACU2M = {
    'UNITLESS' : 1.0,
    'INCHES' : 0.0254,
    'FEET' : 0.3048,
    'MILES' : 1609.34,
    'MILLIMETERS' : 0.001,
    'CENTIMETERS' : 0.01,
    'METERS' : 1.0,
    'KILOMETERS' : 1000.0,
    'MICROINCHES' : 0.0000000254,
    'MILS' : 0.0000254,
    'YARDS' : 0.9144,
    'ANGSTROMS' : 1.0e-10,
    'NANOMETERS' : 1.0e-9,
    'MICRONS' : 1.0e-6,
    'DECIMETERS' : 0.1,
    'DECAMETERS' : 10.0,
    'HECTOMETERS' : 100.0,
    'GIGAMETERS' : 1.0e9,
    'ASTRONOMICAL UNITS' : 1.496e11,
    'LIGHT YEARS' : 9.461e15,
    'PARSECS' : 3.086e16
}

class InputDXF(object):
    """Class for interfacing and defining input DXF data"""
    def __init__(self):

        # ----- Class properties
        self.filename = None
        self.drawingname = None
        self.autocadversion = None
        self.exterior = None
        self.interior = None
        self.markers = None
        self.anchor = None
        self.wh = [0, 0]              # [width, height]
        self.units = None

    def readDXF(self, filename):
        """Reads a DXF file and stores it as an object"""
        dxf = dxfgrabber.readfile(filename)
        print("DXF version: {}".format(ACVERSIONDICT[dxf.dxfversion]))
        entities = dxf.entities
        entity_count = len(dxf.entities)
        print("Reading DXF file: '{}' ...".format(dxf.filename))
        self.filename = dxf.filename

        headers = dxf.header
        self.units = ACUNITS[headers['$INSUNITS']]
        units = ACU2M[self.units]

        # ------ Read list of layers
        layers = dxf.layers.names()
        marker_val = 10
        marker_dict = {}
        for layer in layers:
            marker_dict.update({layer : marker_val})
            marker_val = marker_val + 10
        self.markers = marker_dict

        # ------ Define lists to be filled with line entities (raw), and sorted (exterior, holes, etc.)
        dxf_lines = []              # Entries formatted: [start_tuple, end_tuple, Layer Name, Dxftype, extra]
        dxf_exterior = []           # Entries formatted: [start_tuple, end_tuple, Layer Name, Dxftype, idx]
        dxf_holes = []              # Entries formatted: [start_tuple, end_tuple, Layer Name, Dxftype, idx, hole_num, extra]
        num_ext = 0
        num_holes = 0
        num_lines = 0

        # ------ First, read all the line values and save them to dxf_lines
        for ent in entities:
            if(ent.layer == 'svgframe' and ent.dxftype == 'LWPOLYLINE'):
                # This is the frame text point, use this to anchor text for the future

                # Check to see which point is in the lower left (min x, min y)
                anchor = ent.points[0]
                for point in ent.points:
                    if point[0] - anchor[0] < 0.000000001 and point[1] - anchor[1] < 0.000000001:
                        anchor = point
                # Get width and height
                if abs(ent.points[0][0] - ent.points[1][0]) < 0.000000001:      # x is same, this is the height
                    self.wh[0] = abs(ent.points[0][0] - ent.points[-1][0]) * units
                    self.wh[1] = abs(ent.points[0][1] - ent.points[1][1]) * units
                else:       # reverse is true
                    self.wh[0] = abs(ent.points[0][0] - ent.points[1][0]) * units
                    self.wh[1] = abs(ent.points[0][1] - ent.points[-1][1]) * units
                for point in ent.points:
                    # Now, get width and height
                    if abs(point[0] - anchor[0]) > 0.000000001 and abs(point[1] - anchor[1]) < 0.000000001:
                        self.wh[0] = abs(point[0] - anchor[0]) * units
                    elif abs(point[0] - anchor[0]) < 0.000000001 and abs(point[1] - anchor[1]) > 0.000000001:
                        self.wh[1] = abs(point[1] - anchor[1]) * units
                self.anchor = [i * units for i in anchor]

            if(ent.dxftype == 'LINE'):
                num_lines = num_lines + 1
                start_point = []
                end_point = []
                for x, y in zip(ent.start, ent.end):
                    start_point.append(x * units)
                    end_point.append(y * units)
                start_point = tuple(start_point)
                end_point = tuple(end_point)
                # ------ Add to the number counts and seed the exterior and holes arrays
                # if they haven't been done so already
                if ent.line_weight >= 50:
                    num_ext = num_ext + 1
                    if len(dxf_exterior) == 0:
                        dxf_exterior.append([start_point, end_point, ent.layer, ent.dxftype, 0])
                        continue
                else:
                    num_holes = num_holes + 1
                dxf_lines.append([start_point, end_point, ent.layer, ent.dxftype])

            elif(ent.dxftype == 'CIRCLE'):
                center = []
                radius = ent.radius * units
                for z in ent.center:
                    center.append(z * units)
                center = tuple(center)
                right_point = tuple(sum(x) for x in zip(center, (radius, 0, 0)))
                left_point = tuple(sum(x) for x in zip(center, (-radius, 0, 0)))
                top_point = tuple(sum(x) for x in zip(center, (0, radius, 0)))
                bottom_point = tuple(sum(x) for x in zip(center, (0, -radius, 0)))
                center_point = tuple(center)
                dxf_lines.append([right_point, left_point, ent.layer, ent.dxftype,
                                  center_point, top_point, bottom_point])
                num_holes = num_holes + 1

        # ------ Now, sort the dxf_lines
        while(len(dxf_exterior) < num_ext):
            for idx, ext_line in zip(range(len(dxf_lines)), dxf_lines):
                if(ext_line[0] == dxf_exterior[-1][1]):
                    # The two lines are connected, append the line to the list and restart the loop
                    dxf_exterior.append([ext_line[0], ext_line[1], ext_line[2], ext_line[3], len(dxf_exterior)])
                    del dxf_lines[idx]
                    break
                elif(ext_line[1] == dxf_exterior[-1][1]):
                    # The two lines are connected (but start/end reversed),
                    # append the line to the list and restart the loop
                    dxf_exterior.append([ext_line[1], ext_line[0], ext_line[2], ext_line[3], len(dxf_exterior)])
                    del dxf_lines[idx]
                    break

        current_hole = -1
        current_hole_idx = 0
        circle_add_idx = 0
        while len(dxf_holes) < num_holes:
            # Update the current_hole values, and add a new "seed" for the next set of holes
            if (len(dxf_lines) > 0):

                current_hole = current_hole + 1
                # If the first element in dxf_lines is still a circle, restart the loop
                if dxf_lines[0][3] == 'CIRCLE':
                    # If this is the first time going through this loop (current_hole == 0), set initial idx values
                    if current_hole == 0:
                        dxf_holes.append([dxf_lines[0][0], dxf_lines[0][1], dxf_lines[0][2], dxf_lines[0][3],
                                          num_ext, current_hole, dxf_lines[0][4], dxf_lines[0][5], dxf_lines[0][6]])
                    else:
                        dxf_holes.append([dxf_lines[0][0], dxf_lines[0][1], dxf_lines[0][2], dxf_lines[0][3],
                                          dxf_holes[0][4] + len(dxf_holes) + circle_add_idx, current_hole,
                                          dxf_lines[0][4], dxf_lines[0][5], dxf_lines[0][6]])
                    del dxf_lines[0]
                    circle_add_idx = circle_add_idx + 4
                    continue

                # If this is the first time going through this loop (current_hole == 0), set initial idx values
                if current_hole == 0:
                    dxf_holes.append([dxf_lines[0][0], dxf_lines[0][1], dxf_lines[0][2], dxf_lines[0][3],
                                      num_ext, current_hole])
                else:
                    dxf_holes.append([dxf_lines[0][0], dxf_lines[0][1], dxf_lines[0][2], dxf_lines[0][3],
                                      dxf_holes[0][4] + len(dxf_holes) + circle_add_idx, current_hole])
                del dxf_lines[0]
                current_hole_idx = len(dxf_holes) - 1

            while(dxf_holes[current_hole_idx][0] != dxf_holes[-1][1]):
                # --- While the end coordinates are not equal to the start coordinates for the current hole, hole has
                #       not been completed yet. Must continue to iterate

                for idx, int_line in zip(range(len(dxf_lines)), dxf_lines):
                    if (int_line[0] == dxf_holes[-1][1]):
                        # The two lines are connected, append the line to the list and restart the loop
                        dxf_holes.append([int_line[0], int_line[1], int_line[2], int_line[3],
                                        dxf_holes[0][4] + len(dxf_holes) + circle_add_idx, current_hole])
                        del dxf_lines[idx]
                        break
                    elif (int_line[1] == dxf_holes[-1][1]):
                        # The two lines are connected (but start/end reversed),
                        # append the line to the list and restart the loop
                        dxf_holes.append([int_line[1], int_line[0], int_line[2], int_line[3],
                                          dxf_holes[0][4] + len(dxf_holes) + circle_add_idx, current_hole])
                        del dxf_lines[idx]
                        break

        self.exterior = dxf_exterior
        self.interior = dxf_holes
        return dxf_exterior, dxf_holes


    def convertToGeometry(self, max_el_size=1):
        """Takes exterior and interior line lists generated from readDXF to convert to Calfem-readable geometry"""
        # ------ First, store the dxf lines as local variables
        outline = self.exterior
        holes = self.interior

        # ------ Create a geometry instance to store our geometry parameters
        g = cfg.geometry()

        # ------ Create points for the model
        for line in outline:
            g.addPoint([line[0][0], line[0][1]], ID=line[4])
        for hole in holes:
            # Check to see if it's a circle, if so, do things differently
            if hole[3] == 'CIRCLE':
                g.addPoint([hole[0][0], hole[0][1]], ID=hole[4], elSize=max_el_size * 1.5)              # Right side
                g.addPoint([hole[7][0], hole[7][1]], ID=hole[4]+1, elSize=max_el_size * 1.5)            # Top side
                g.addPoint([hole[1][0], hole[1][1]], ID=hole[4]+2, elSize=max_el_size * 1.5)            # Left side
                g.addPoint([hole[8][0], hole[8][1]], ID=hole[4]+3, elSize=max_el_size * 1.5)            # Bottom side
                g.addPoint([hole[6][0], hole[6][1]], ID=hole[4]+4, elSize=max_el_size * 1.5)            # Center
            else:
                g.addPoint([hole[0][0], hole[0][1]], ID=hole[4], elSize=max_el_size * 1.5)

        # ------ Create links between the points of the model

        # --- First, define a dictionary to fill with all curves belonging to each layer
        curve_list = []
        for i in range(len(self.markers)):
            curve_list.append([])
        curve_dict = {}
        exterior_splines = []

        # Define exterior lines
        for line in outline:
            if(line[4] == len(outline) - 1):
                # For the final element, make sure to loop back to zero
                g.spline([line[4], outline[0][4]], ID=line[4], marker=self.markers[line[2]])
            else:
                g.spline([line[4], line[4] + 1], ID=line[4], marker=self.markers[line[2]])
            # Convert marker number to curve list idx num, i.e.: 10 -> 0, 20 -> 1, etc.
            curve_list[int(self.markers[line[2]]/10 - 1)].append(line[4])
            exterior_splines.append(line[4])

        # Define interior lines (holes)
        hole_num = 0
        hole_start_idx = 0
        hole_temp = []
        hole_splines = []
        for idx, hole in zip(range(len(holes)), holes):
            # Check if it's a circle, if so, add the necessary splines
            if hole[3] == 'CIRCLE':
                g.addCircle([hole[4], hole[4] + 4, hole[4] + 1], ID=hole[4],
                            marker=self.markers[hole[2]])   # Right to top
                g.addCircle([hole[4] + 1, hole[4] + 4, hole[4] + 2], ID=hole[4] + 1,
                            marker=self.markers[hole[2]])   # Top to left
                g.addCircle([hole[4] + 2, hole[4] + 4, hole[4] + 3], ID=hole[4] + 2,
                            marker=self.markers[hole[2]])   # Left to bottom
                g.addCircle([hole[4] + 3, hole[4] + 4, hole[4]], ID=hole[4] + 3,
                            marker=self.markers[hole[2]])   # Bottom to left
                curve_list[int(self.markers[hole[2]]/10 - 1)].append(hole[4])
                curve_list[int(self.markers[hole[2]]/10 - 1)].append(hole[4] + 1)
                curve_list[int(self.markers[hole[2]]/10 - 1)].append(hole[4] + 2)
                curve_list[int(self.markers[hole[2]]/10 - 1)].append(hole[4] + 3)
                hole_num = hole_num + 1
                hole_temp.append(hole[4])
                hole_temp.append(hole[4] + 1)
                hole_temp.append(hole[4] + 2)
                hole_temp.append(hole[4] + 3)
                hole_splines.append(hole_temp)
                hole_temp = []
                hole_start_idx = idx + 1
                continue
            elif hole[1] == holes[hole_start_idx][0]:
                # For the final element, make sure to loop back to the start idx
                g.spline([hole[4], holes[hole_start_idx][4]], ID=hole[4], marker=self.markers[hole[2]])
                curve_list[int(self.markers[hole[2]]/10 - 1)].append(hole[4])
                hole_num = hole_num + 1
                hole_temp.append(hole[4])
                hole_splines.append(hole_temp)
                hole_temp = []
                hole_start_idx = idx + 1
                continue
            else:
                g.spline([hole[4], hole[4] + 1], ID=hole[4], marker=self.markers[hole[2]])
                curve_list[int(self.markers[hole[2]]/10 - 1)].append(hole[4])
            hole_temp.append(hole[4])

        # ------ Define the surface on which the mesh is to be generated
        g.addSurface(exterior_splines, hole_splines)

        # ------ Convert the curve list to a dictionary
        for mark, curves in zip(self.markers, curve_list):
            curve_dict.update({mark : curves})

        # ------ Return the created geometry and curve dictionary
        return g, curve_dict


    def adjustDimension(self, marker, value):
        """Takes the numerical 'value' and adjusts all lines of the given marker to the value (assumes equal extensions
        in the same direction that the line currently is in)"""
        if value == -1 or marker not in self.markers.keys():
            return
        for idx, ext_line in zip(range(len(self.exterior)), self.exterior):
            if(ext_line[2] == marker):
                dy = ext_line[1][1] - ext_line[0][1]
                dx = ext_line[1][0] - ext_line[0][0]
                if dx == 0 and dy == 0:
                    continue
                old_length = np.sqrt(dy ** 2 + dx ** 2).item()
                angle = np.arctan2(dy, dx).item()
                delta = value - old_length
                deltax = delta * np.cos(angle).item()
                deltay = delta * np.sin(angle).item()

                # ---- Add half the delta value to each end of the line
                self.exterior[idx][1] = (ext_line[1][0] + deltax / 2, ext_line[1][1] + deltay / 2, 0.0)
                self.exterior[idx][0] = (ext_line[0][0] - deltax / 2, ext_line[0][1] - deltay / 2, 0.0)

                # ---- Do the same manipulations to the adjoining line segments
                self.exterior[idx - 1][1] = (self.exterior[idx - 1][1][0] - deltax / 2,
                                             self.exterior[idx - 1][1][1] - deltay / 2, 0.0)

                if(idx == len(self.exterior)):
                    # ---- Make sure the idx doesn't go out of bounds if at end of the list
                    self.exterior[0][0] = (self.exterior[0][0][0] + deltax / 2,
                                           self.exterior[0][0][1] + deltay / 2, 0.0)
                else:
                    self.exterior[idx + 1][0] = (self.exterior[idx + 1][0][0] + deltax / 2,
                                             self.exterior[idx + 1][0][1] + deltay / 2, 0.0)

        hole_start_idx = 0
        hole_end_idx = 0
        hole_num = -1
        for idx, int_line in zip(range(len(self.interior)), self.interior):
            # ---- If we are at a new hole, update the hole variables accordingly
            if int_line[5] != hole_num:
                hole_num = hole_num + 1
                hole_start_idx = idx
                # ---- First, check to see if we're at a CIRCLE, if so, pass
                if int_line[3] != 'CIRCLE':
                    # ---- Check to see if we're at the end of the hole list
                    if self.interior[-1][5] == hole_num:
                        hole_end_idx = len(self.interior) - 1
                    else:
                        for i, line in zip(range(len(self.interior)), self.interior):
                            if line[5] != hole_num:
                                hole_end_idx = i - 1
                                break

            if int_line[2] == marker:
                dy = int_line[1][1] - int_line[0][1]
                dx = int_line[1][0] - int_line[0][0]
                if dx == 0 and dy == 0:
                    continue
                old_length = np.sqrt(dy ** 2 + dx ** 2).item()
                angle = np.arctan2(dy, dx).item()
                delta = value - old_length
                deltax = delta * np.cos(angle).item()
                deltay = delta * np.sin(angle).item()

                # ---- Add half the delta value to each end of the line
                self.interior[idx][1] = (int_line[1][0] + deltax / 2, int_line[1][1] + deltay / 2, 0.0)
                self.interior[idx][0] = (int_line[0][0] - deltax / 2, int_line[0][1] - deltay / 2, 0.0)

                # ---- If a CIRCLE, do special set of calculations (adjust the top, bottom values as well)
                if int_line[3] == 'CIRCLE':
                    self.interior[idx][7] = (int_line[7][0] - deltay / 2, int_line[7][1] - deltax / 2, 0.0) # top
                    self.interior[idx][8] = (int_line[8][0] + deltay / 2, int_line[8][1] + deltax / 2, 0.0) # bottom

                # ---- Do the same manipulations to the adjoining line segments if NOT a circle
                else:
                    if(idx == hole_start_idx):
                        # ---- Make sure the idx doesn't go out of bounds if at end of the list
                        self.interior[hole_end_idx][1] = (self.interior[hole_end_idx][1][0] - deltax / 2,
                                                         self.interior[hole_end_idx][1][1] - deltay / 2, 0.0)
                        self.interior[idx + 1][0] = (self.interior[idx + 1][0][0] + deltax / 2,
                                                     self.interior[idx + 1][0][1] + deltay / 2, 0.0)
                    elif(idx == hole_end_idx):
                        self.interior[idx - 1][1] = (self.interior[idx - 1][1][0] - deltax / 2,
                                                     self.interior[idx - 1][1][1] - deltay / 2, 0.0)
                        self.interior[hole_start_idx][0] = (self.interior[hole_start_idx][0][0] + deltax / 2,
                                                            self.interior[hole_start_idx][0][1] + deltay / 2, 0.0)
                    else:
                        self.interior[idx - 1][1] = (self.interior[idx - 1][1][0] - deltax / 2,
                                                     self.interior[idx - 1][1][1] - deltay / 2, 0.0)
                        self.interior[idx + 1][0] = (self.interior[idx + 1][0][0] + deltax / 2,
                                                     self.interior[idx + 1][0][1] + deltay / 2, 0.0)

    def convertDXFtoSVG(self):
        self.drawingname = save_svg_from_dxf(self.filename, frame_name="drawing")



if __name__ == '__main__':

    # --- Create InputDXF instance
    inputdxf = InputDXF()
    exterior, interior = inputdxf.readDXF("frame1.dxf")
    inputdxf.convertToGeometry()
    inputdxf.adjustDimension('a', 6.0)
    print(inputdxf.exterior)
    print(inputdxf.interior)
    inputdxf.convertDXFtoSVG()