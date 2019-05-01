# -*- coding: utf-8 -*-
#******************************************************************************
#
# bflat.py
#        
bflat_VERSION = "0.2.0"
# ---------------------------------------------------------
# Python GUI module for flattening a X-Plane mesh at a given airport.
#
# For more details refert to GitHub: https://github.com/nofaceinbook/betterflat
#
# WARNING: This code is still under development and may still have some errors.
#
# Copyright (C) 2019 by schmax (Max Schmidt)
#
# This code is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.  
#
# A copy of the GNU General Public License is available at:
#   <http://www.gnu.org/licenses/>. 
#
#******************************************************************************

import logging
from xplnedsf import *
from bflatKMLexport import *
import os
from tkinter import *
from tkinter.filedialog import askopenfilename
from math import sin, cos, sqrt, atan2, radians


def displayHelp(win):
    newwin = Toplevel(win)
    Label(newwin, anchor=W, justify=LEFT, text=
          "This program flattens the mesh of X-Plane at an given airport.\n"
          "It finds all tringles that either intersect the airport boundary \n"
          "or the runways, both defined in an airport data file apt.dat.\n"
          "and sets them all to a given height.\n\n"
          "So first select the apt.dat including your airport. In case this apt.dat\n"
          "includes several airport you have to specify it with the identifier\n"
          "usually the 4 letter ICAO code.\n"
          "Next you have to specify and read the dsf-file including the mesh for the\n"
          "area of the airport.\n" 
          "Now you can write and later use the updated dsf-file. You can\n"
          "also directly specify a backup file of the original dsf-file.\n\n"
          "With the .kml-button you can generate a kml-file where you can\n"
          "view the mesh and adapted triangles e.g. in Google Earth.\n\n"
          "WARNING: Existing files are directly overwritten.\n"
          "Strongly recommended to BACKUP original files first.\n\n" 
          "This program is published under GNU General Public License.\n"
          "So you have permission for us but NO Warrenty and NO Liability.\n"
          "Info on license of used 3rd party libraries included in separate\n"
          "license.txt or refer to info below.\n\n"
          "MORE INFORMATION, source code and contact info are\n"
          "available at GitHub: https://github.com/nofaceinbook/betterflat/\n\n"
          "Hope the tool helps you.    (c) 2019 by schmax (Max Schmidt)"   
          ).grid(row=0, pady=10, padx=10)


def defineLog(logname, logLevelStream='INFO', logLevelFile='INFO', mode='w+'):
    #
    # defines log-Handler, if Level is set to None no File/Stream handler will be created
    # The file for the stream is written to logname.log
    # Global variable LogName is defined to be used e.g. to create sub-loggers
    # Returns the created logger log to be used
    #
    global LogName
    LogName = logname
    directory = os.path.dirname(os.path.abspath(__file__))
    logfile = os.path.join(directory, logname + '.log')
    log = logging.getLogger(logname)
    log.handlers = []  # Spyder/IPython currently does not remove existing loggers, this does the job
    if logLevelStream == 'DEBUG' or logLevelFile == 'DEBUG':
        log.setLevel('DEBUG') #so maximum level is DEBUG; needed to see with getEffectiveLevel if DEBUG is enabled or not  
    else:
        log.setLevel('INFO') #set maximum level if not DEBUG before restriction by handlers below
    formatter = logging.Formatter('%(name)s-%(levelname)s: %(message)s')
    if logLevelStream:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logLevelStream)
        stream_handler.setFormatter(formatter)
        log.addHandler(stream_handler)
    if logLevelFile:
        file_handler = logging.FileHandler(logfile, mode)
        file_handler.setLevel(logLevelFile)
        file_handler.setFormatter(formatter)
        log.addHandler(file_handler)
    return log


def readAPT(filename, icao_id=""):
    #
    # Reads boundaries, runways and airport height from airport in given apt.dat file and returns them
    # In case apt.dat contains several airports the right airport is selected by giving the correct ICAO Code
    # In case airport has no icoa code the identifier, use the identiefer on 5th position in airport defintion
    # Each boundary is list of [lon, lat] values of vertices of the boundary
    # Hole definition in boundary are treated the same as a boundary!!!
    # For Bezier nodes of the boundary only the node without Bezier definitions is considered!!!
    # Only land runways  (type 100) are returned in a list with sub-lists of their endpoints, width in meter
    # Returned height is an integer in meter, also retruns flatten flag if set for the airport
    #
    log.info("Reading airport data from: {}".format(filename))
    Airport = False  # else first the correct icoa id has to be found before Airport becomes true
    BoundarySection = False  # is set true when we are in boundary section in .apt file
    bounds = []  # list of list with all boundaries found
    runways = [] # list of runway endpoints
    apt_elev = None # elevation of airport in meters
    apt_flatten = None #includes entry if this airport has a flatten flag set
    apt_name = None
    if not os.path.isfile(filename):
        log.error("Airport File {} does not exist!".format(filename))
        return [], "Error: Airport file does not exist!"
    with open(filename, encoding="utf8", errors="ignore") as f:
        for line in f:
            v = line.split()
            if len(v) == 0:  # don't consider empty lines
                continue
            if len(v) > 4:  # check if correct airport section in file is reached
                if v[0] == '1':
                    if v[4] == icao_id or icao_id =='': #if no icao id is given just first airport is selected
                        Airport = True
                        icao_id = v[4] #set now icao id in case it was '' before
                        apt_elev = round(int(v[1]) * 0.3048)
                        #v.append(" ")
                        #log.info(v)
                        apt_name = v[5] 
                        if len(v) > 6: apt_name = apt_name + " " + v[6]
                        log.info("Airport {} found with elevation {} m.".format(apt_name, apt_elev))
                    else:
                        Airport = False  # change to false in case of new different airport
            if Airport: 
                if v[0] == '130':
                    BoundarySection = True
                    bounds.append([])  # add new list of boundary vertices
                elif v[0] == '100':
                    log.info("Runway from {}, {} to {}, {} with width {} found".format(v[9], v[10], v[18], v[19], v[1]))
                    runways.append( [(float(v[9]), float(v[10])), (float(v[18]), float(v[19])), float(v[1]) ])
                elif v[0] == '1302' and v[1] == 'flatten':
                    apt_flatten = int(v[2])
                    log.warning("Airport includes flatten flag set to: {}".format(apt_flatten))
                elif BoundarySection:
                    if v[0] == '111' or v[0] == '112':
                        bounds[-1].append(
                            [float(v[2]), float(v[1])])  # Note: Bezier definitins are not considered, just the base point
                    elif v[0] == '113' or v[0] == '114':
                        bounds[-1].append(
                            [float(v[2]), float(v[1])])  # Note: Bezier definitins are not considered, just the base point
                        bounds[-1].append(bounds[-1][0])  # #form closed loop by adding again first vertex
                        BoundarySection = False
                        log.info("Boundary no. {} with {} vertices read.".format(len(bounds), len(bounds[-1])))
    if len(bounds) == 0 and len(runways) == 0:
        log.warning("No valid boundary or runway found in file!")
        return [], [], None, None, None, "Warning: No valid boundary or runway found in file!"
    else:
        log.info("Finished reading boundaries.")
        return bounds, runways, apt_elev, apt_flatten, apt_name, None

def distance (p1, p2): #calculates distance between p1 and p2 in meteres where p is pair of longitude, latitude values
    R = 6371009 #mean radius earth in m
    lon1 = radians(p1[1]) ##acutally are lon and lat mixed here, but works like that
    lat1 = radians(p1[0])
    lon2 = radians(p2[1])
    lat2 = radians(p2[0]) 
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))   
    return R * c

def getRunwayBounds (p1, p2, w):
    degree_dist_at_equator = 111120 #for longitude (or 111300?)
    lat_degree_dist_average = 111000
    degree_dist_at_lat = cos (radians(p1[0])) * degree_dist_at_equator
    if round (p1[1], 6) == round (p2[1], 6): #runway exactly east-west direction
        dx = 0   #difference for longitute in meters to reach corner from center end
        dy = w/2 #difference for latitude in meters to reach corner from center end
    elif round (p1[0], 6) == round (p2[0], 6): #runway is exactly north-south direction
        dx = w/2
        dy = 0
    else: 
        m = -1 / ((p2[0] - p1[0]) / (p2[1] - p1[1])) #gradient of perpendicular line Steigung der Senkrechten
        dx = sqrt( ( (w/2)**2) / (1 + m**2)) 
        dy = dx * m 
    dx /= degree_dist_at_lat #convert meters in longitute coordinate difference at geographical latitude
    dy /= lat_degree_dist_average #convert meters in latitude coordinate difference 
    l = []
    l.append([round(p1[1] - dx, 8), round(p1[0] - dy, 8)])
    l.append([round(p1[1] + dx, 8), round(p1[0] + dy, 8)])
    l.append([round(p2[1] + dx, 8), round(p2[0] + dy, 8)])
    l.append([round(p2[1] - dx, 8), round(p2[0] - dy, 8)])
    l.append(l[0]) #add firs corner to form closed loop
    return l

def _linsolve_(a1, b1, c1, a2, b2, c2):
    divisor = (a1 * b2) - (a2 * b1)
    if divisor == 0:
        return (-99999, -99999)  # actually None but with negative values calling intersection function returns None
    return (round(((c1 * b2) - (c2 * b1)) / divisor, 8),
            round(((a1 * c2) - (a2 * c1)) / divisor, 8))  # ROUNDING TO ALWAYS GET SAME CONCLUSION


def intersection(p1, p2, p3, p4):  # checks if segment from p1 to p2 intersects segement from p3 to p4
    s0, t0 = _linsolve_(p2[0] - p1[0], p3[0] - p4[0], p3[0] - p1[0], p2[1] - p1[1], p3[1] - p4[1], p3[1] - p1[1])
    if s0 >= 0 and s0 <= 1 and t0 >= 0 and t0 <= 1:
        return (round((p1[0] + s0 * (p2[0] - p1[0])), 8), round(p1[1] + s0 * (p1[1] - p1[1]),
                                                                8))  ### returns the cutting point as tuple; ROUNDING TO ALWAYS GET SAME POINT
    else:
        return (None)


def PointInPoly(p, poly):  # test wether a point p with [lat, lon] coordinates lies in polygon (list of [lat, lon] pairs
    # counts number of intersections from point outside poly to p on same y-coordinate, if it is odd the point lies in poly
    # to avoid intersection at vertex of poly on same y-coordinate, such points are shifte about 1mm above for testing intersection
    count = 0
    for i in range(len(poly) - 1):  # for all segments in poly
        epsilon0, epsilon1 = (
        0, 0)  # added to p's y coordinate in case p is on same y-coordinate than according vertex of segment
        if poly[i][1] < p[1] and poly[i + 1][1] < p[
            1]:  # if vertices of segment below y-coordinate of p, no intersection
            continue
        if poly[i][1] > p[1] and poly[i + 1][1] > p[
            1]:  # if vertices of segment above y-coordinate of p, no intersection
            continue
        if poly[i][1] == p[1]:
            epsilon0 = 0.00000001
        if poly[i + 1][1] == p[1]:
            epsilon1 = 0.00000001
        x = intersection([poly[i][0], poly[i][1] + epsilon0], [poly[i + 1][0], poly[i + 1][1] + epsilon1], [181, p[1]],
                         p)
        if x:
            count += 1
    if count % 2:
        return True  # odd number of intersections, so p in poly
    else:
        return False  # even number of intersection, so p outside poly


def vertices_of_boundary_intersecting_trias(dsf, poly):
    #
    # Identifies all mesh triangles of dsf file having an intersection with an edge of the boundary (poly with [lon, lat] coordinates)
    # or lie inside boundary (or single triangle in which boundary completely lies)
    # Returns a set with all vertices (tuple of pool_id and vertex number in pool) that belong to such intersecting triangles
    #
    log.info("Reading of dsf-file completed. Checking now for mesh triangles intersecting airport boundary.")
    miny, maxy, minx, maxx = dsf.BoundingRectangle(poly)
    s = set([])
    for p in dsf.PatchesInArea(miny, maxy, minx, maxx):
        for t in p.triangles():
            TriaV = dsf.TriaVertices(t)
            TriaV.append(TriaV[0])  # append first vertex to get a closed connection
            for i in range(3):  # 3 sides of triangle
                for j in range(len(poly) - 1):  # sides of poly
                    if intersection(TriaV[i], TriaV[i + 1], poly[j],
                                    poly[j + 1]):  # current triangle intersects with an poly line
                        s.add((t[0][0], t[0][1]))  # add all 3 vertices of tria to set as tuples
                        s.add((t[1][0], t[1][1]))
                        s.add((t[2][0], t[2][1]))
            if PointInPoly(TriaV[0],
                           poly):  # Check also that not complete Tria lies in poly by checking for first vertex of Tria
                s.add((t[0][0], t[0][1]))  # add all 3 vertices of tria to set as tuples
                s.add((t[1][0], t[1][1]))
                s.add((t[2][0], t[2][1]))
            if PointInPoly(poly[0],
                           TriaV):  # Check also that not complete poly lies in current tria by checking for first vertex of poly
                s.add((t[0][0], t[0][1]))  # add all 3 vertices of tria to set as tuples
                s.add((t[1][0], t[1][1]))
                s.add((t[2][0], t[2][1]))
    log.info("{} vertices of mesh found that belong to mesh triangles intersecting or within boundary.".format(len(s)))
    changedVertices = {}  # set up dictonary with all coords of changed vertices to find additional vertices at such coords for trias out of bound
    for v in s:
        changedVertices[
            (round(dsf.V[v[0]][v[1]][0], 7), round(dsf.V[v[0]][v[1]][1], 7))] = True  # coords round to range of cm
    miny, maxy, minx, maxx = dsf.BoundingRectangle(changedVertices)
    delta = 0.000001  # check also for vertices 0.1m outside boundary
    for p in dsf.PatchesInArea(miny - delta, maxy + delta, minx - delta, maxx + delta):
        for t in p.triangles():
            for v in t:
                if (round(dsf.V[v[0]][v[1]][0], 7), round(dsf.V[v[0]][v[1]][1], 7)) in changedVertices:
                    s.add((v[0], v[1]))  # add tuple of vertex at coords to make sure to get all vertices at the coords
    log.info("References to vertices to be changed: {}".format(s))
    log.info("Unique coords to be changed: {}".format(changedVertices))
    return s, changedVertices


def averageheight(dsf, vertices):
    count = 0
    sum = 0
    for v in vertices:
        count += 1
        sum += dsf.getElevation(dsf.V[v[0]][v[1]][0], dsf.V[v[0]][v[1]][1], dsf.V[v[0]][v[1]][2])
    if count == 0:
        return None
    else:
        return round(sum / count)


class bflatGUI:
    def __init__(self):

        self.dsf = XPLNEDSF(LogName, self.showProgress)  # includes all information about the read dsf-file
        self.boundaries = []  # set of bouundaries where each boundary includes vertices of boundary
        self.runways = [] #list of runways for processed airport (endpoints, width in m)
        self.apt_elev = None #airport elevation as defined in apt.dat in m
        self.vertices = set([])  # set of vertices that should get new height
        self.dsfreadfile = ""  # name of the file that is stored in self.dsf
        self.current_action = None #is set to 'read' or 'write' when reading/writing dsf

        self.window = Tk()
        self.window.title("X-Plane bflat (version: {})".format(bflat_VERSION))

        self.header = Label(self.window, text="Make sure you always have copies of your original files!")
        self.header.grid(row=0, column=0, columnspan=2)
        self.help_button = Button(self.window, text=' Help ', fg='red', command=lambda: displayHelp(self.window))
        self.help_button.grid(row=0, column=3, sticky=W, pady=4, padx=10)

        self.aptfile_label = Label(self.window, text="Airport File (apt.dat):")
        self.aptfile_label.grid(row=1, column=0, sticky=W)
        self.apt_entry = Entry(self.window, width=60)
        self.apt_entry.grid(row=1, column=1, columnspan=2, sticky=W)
        self.apt_select = Button(self.window, text='Select', command=lambda: self.select_file(self.apt_entry))
        self.apt_select.grid(row=1, column=3, sticky=W, pady=4, padx=10)

        self.aptid_label = Label(self.window, text="Airport ICAO Id:")
        self.aptid_label.grid(row=2, column=0, sticky=W)
        self.aptid_entry = Entry(self.window, width=6)
        self.aptid_entry.grid(row=2, column=1, sticky=W)
        self.aptid_info = Label(self.window, text="")
        self.aptid_info.grid(row=2, column=2, sticky=W)
        
        self.boundtype = StringVar()
        self.boundtype_label = Label(self.window, text="Boundary Type:")
        self.boundtype_label.grid(row=3, column=0, sticky=W)
        self.boundtype_radioA = Radiobutton(self.window, text="Airport", variable=self.boundtype, value="airport")
        self.boundtype_radioA.grid(row=3, column=1)
        self.boundtype_radioR = Radiobutton(self.window, text="Runways", variable=self.boundtype, value="runways")
        self.boundtype_radioR.grid(row=3, column=2)
        self.boundtype_radioA.select()
        
        self.apt_read = Button(self.window, text='Read Boundary',
                               command=lambda: self.read_apt(self.apt_entry.get(), self.aptid_entry.get()))
        self.apt_read.grid(row=4, column=0, sticky=E, pady=4)
        self.apt_status_label = Label(self.window, text="ICOA Id only required if file contains several airports!")
        self.apt_status_label.grid(row=4, column=1, columnspan=2, sticky=W)

        self.dsffile_label = Label(self.window, text="DSF File:")
        self.dsffile_label.grid(row=5, column=0, sticky=W)
        self.dsf_entry = Entry(self.window, width=60)
        self.dsf_entry.grid(row=5, column=1, columnspan=2, sticky=W)
        self.dsf_select = Button(self.window, text='Select', command=lambda: self.select_file(self.dsf_entry)).grid(
            row=5, column=3, sticky=W, pady=4, padx=10)

        self.dsf_read = Button(self.window, text='  Read DSF   ', state=DISABLED,
                               command=lambda: self.read_dsf(self.dsf_entry.get()))
        self.dsf_read.grid(row=6, column=0, sticky=E, pady=4)
        self.dsf_status_label = Label(self.window, text="   this can take some minutes....")
        self.dsf_status_label.grid(row=6, column=1, sticky=W)

        self.result_text = Label(self.window, text="Result:")
        self.result_text.grid(row=7, column=0, sticky=E, pady=8)
        self.result_label = Label(self.window, text=" not yet")
        self.result_label.grid(row=7, column=1, sticky=W, pady=8)
        self.kml_create = Button(self.window, text='.kml', state=DISABLED, command=lambda: kmlExport(self.dsf, self.boundaries, self.vertices, self.dsf_entry.get()))
        self.kml_create.grid(row=7, column=3, sticky=W, pady=4, padx=10)
        self.height_label = Label(self.window, text="New height (in m):")
        self.height_label.grid(row=8, column=0, pady=4, sticky=E)
        self.height_entry = Entry(self.window, width=7)
        self.height_entry.grid(row=8, column=1, sticky=W)

        self.newdsf_label = Label(self.window, text="Updated DSF file:")
        self.newdsf_label.grid(row=9, column=0, sticky=E, pady=4)
        self.newdsf_entry = Entry(self.window, width=60)
        self.newdsf_entry.grid(row=9, column=1, columnspan=2, sticky=W)
        self.newdsf_change = Button(self.window, text='Change', command=lambda: self.select_file(self.newdsf_entry))
        self.newdsf_change.grid(row=9, column=3, sticky=W, pady=4, padx=10)

        self.bakdsf_label = Label(self.window, text="Backup orig. DSF:")
        self.bakdsf_label.grid(row=10, column=0, sticky=E, pady=4)
        self.bakdsf_entry = Entry(self.window, width=60)
        self.bakdsf_entry.grid(row=10, column=1, columnspan=2, sticky=W)
        self.bakdsf_change = Button(self.window, text='Change', command=lambda: self.select_file(self.bakdsf_entry))
        self.bakdsf_change.grid(row=10, column=3, sticky=W, pady=4, padx=10)

        self.write_button = Button(self.window, text="Write DSF", state=DISABLED,
                                   command=lambda: self.write_dsf(self.height_entry.get(), self.newdsf_entry.get(),
                                                                  self.bakdsf_entry.get()))
        self.write_button.grid(row=11, column=0, sticky=E, pady=4)
        self.write_status_label = Label(self.window, text="Warning: existing files are overwritten!! Takes time...")
        self.write_status_label.grid(row=11, column=1, sticky=W, pady=4)

        log.info("GUI is set up.")
        mainloop()

    def showProgress(self, percentage):
        if self.current_action == 'read':
            self.dsf_status_label.config(text = "read {} percent".format(percentage))
        elif self.current_action == 'write':
            self.write_status_label.config(text = "written {} percent".format(percentage))
        self.window.update()

    def read_apt(self, filename, icao):
        self.boundaries, self.runways, self.apt_elev, flatten_flag, apt_name, err = readAPT(filename, icao)
        if err:
            self.apt_status_label.config(text=err)
        else:
            self.apt_status_label.config(
                text="{}: {} boundary(ies) and {} runway(s) read.".format(apt_name, len(self.boundaries), len(self.runways)))
            self.dsf_read.config(state="normal")

    def read_dsf(self, filename):
        self.current_action = 'read'
        self.dsfreadfile = filename
        err = self.dsf.read(filename)  # returns value > 0 in case of errors
        if err == 1:
            self.dsf_status_label.config(text="Error: File not found!")
        elif err == 2:
            self.dsf_status_label.config(text="Error: File is 7zipped. Unzip first!")
        elif err == 3:
            self.dsf_status_label.config(text="Error: File is not correct dsf format!")            
        else:
            self.dsf_status_label.config(
                text="DSF with {} pools of vertices and {} patches read.".format(len(self.dsf.V),
                                                                                 len(self.dsf.Patches)))
            self.vertices = set([])  # set to empty in case of previous computations
            coords = {}
            if self.boundtype.get() == 'airport':
                bounds = self.boundaries
                log.info("Using airport boundaries for calculating area to be flattened.")
            else:
                bounds = []
                for r in self.runways:
                    bounds.append(getRunwayBounds(r[0], r[1], r[2] + 2)) #add 2m width for shoulders
                log.info("Using runway boundaries for calculating area to be flattened.")
            for boundary in bounds:  # get intersections for all boundaries
                self.result_label.config(text="Calculate intersections for boundary. Wait....")
                self.window.update()
                vs, cs = vertices_of_boundary_intersecting_trias(self.dsf, boundary)
                self.vertices = self.vertices.union(vs)
                coords.update(cs)
            ah = averageheight(self.dsf, self.vertices)
            if ah == None:
                self.result_label.config(text="No mesh intersection found. Check that files are correct.")
            else:
                self.write_button.config(state="normal")
                self.kml_create.config(state="normal")
                self.height_entry.delete(0, END)
                if self.apt_elev:
                    self.height_entry.insert(0, self.apt_elev)
                self.result_label.config(
                    text="{} vertices at {} coords with average height {} to be adapted.".format(len(self.vertices),
                                                                                                 len(coords), ah))
                self.newdsf_entry.delete(0, END)
                self.newdsf_entry.insert(0, filename)
                self.bakdsf_entry.delete(0, END)
                self.bakdsf_entry.insert(0, filename + ".bak")
        self.current_action = None

    def write_dsf(self, newheight, dsffile, bakfile):
        self.current_action = 'write'
        newheight = int(self.height_entry.get())
        if bakfile != "":
            log.info("Moving original dsf file: {} to backup-file: {}".format(dsffile, bakfile))
            try:
                os.replace(self.dsfreadfile, bakfile)
            except:
                log.error('{} can not be replaced!'.format(self.dsfreadfile))
                self.write_status_label.config(text="Error: Original file {} can not be replaced!".format(self.dsfreadfile))
                self.current_action = None
                return
        log.info("Writing dsf file: {} for height: {}".format(dsffile, newheight))
        for v in self.vertices:
            log.debug(" {} with height {} changed to: ".format(self.dsf.V[v[0]][v[1]], self.dsf.getElevation(self.dsf.V[v[0]][v[1]][0], self.dsf.V[v[0]][v[1]][1], self.dsf.V[v[0]][v[1]][2]), newheight))
            self.dsf.V[v[0]][v[1]][2] = newheight
        self.write_status_label.config(text="Wrtiting changes ...") 
        self.window.update()
        self.dsf.write(dsffile)
        self.write_status_label.config(text="Done.")
        log.info("Done.")
        self.current_action = None

    def select_file(self, entry):
        filename = askopenfilename()
        entry.delete(0, END)
        entry.insert(0, filename)

log = defineLog('bflat', 'INFO', 'INFO')
log.info("Started bflat Version: {}".format(bflat_VERSION))
main = bflatGUI()
