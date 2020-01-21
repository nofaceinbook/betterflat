# -*- coding: utf-8 -*-
#******************************************************************************
#
# bflat.py
#        
bflat_VERSION = "0.4.1 exp"
# ---------------------------------------------------------
# Python GUI module for flattening a X-Plane mesh at a given airport.
#
# For more details refert to GitHub: https://github.com/nofaceinbook/betterflat
#
# WARNING: This code is still under development and may still have some errors.
#
# Copyright (C) 2020 by schmax (Max Schmidt)
#
# This code is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.  
#
# A copy of the GNU General Public License is available at:
#   <http://www.gnu.org/licenses/>. 
#
#******************************************************************************


from logging import StreamHandler, FileHandler, getLogger, Formatter
from xplnedsf import *
from bflatKMLexport import *
from os import path, replace
from shutil import copy2
from tkinter import *
from tkinter.filedialog import askopenfilename
from math import cos, sqrt, radians #for getRunwayBounds calculations

def displayHelp(win):
    newwin = Toplevel(win)
    Label(newwin, anchor=W, justify=LEFT, text=
          "This program flattens the mesh of X-Plane at an given airport.\n"
          "It finds all tringles that either intersect the airport boundary \n"
          "or the runways, both defined in an airport data file apt.dat.\n"
          "and sets them all to a given height.\n"
          "In case you tick the option cut, the triangles will also be cut in new\n"
          "ones that build up the boundary or runway.\n\n"
          "So first select the apt.dat including your airport. In case this apt.dat\n"
          "includes several airport you have to specify it with the identifier\n"
          "usually the 4 letter ICAO code.\n"
          "Next you have to specify and read the dsf-file including the mesh for the\n"
          "area of the airport.\n"
          "IMPORTANT: The dsf-file you need to update is either the default one\n"
          "in X-Plane Global Scenery folder or if you installed mesh, like HD mesh\n"
          "it will be located in the according folder for them mesh.\n"
          "Now you can write and later use the updated dsf-file. You can\n"
          "also directly specify a backup file of the original dsf-file.\n\n"
          "With the .kml-button you can generate a kml-file where you can\n"
          "view the mesh and adapted triangles e.g. in Google Earth.\n\n"
          "In the CONFIGURATION window you can change settings for the\n"
          "cut option. Details for that you have to check in docu.\n\n"
          "WARNING: Existing files are directly overwritten.\n"
          "Strongly recommended to BACKUP original files first.\n\n" 
          "This program is published under GNU General Public License.\n"
          "So you have permission for us but NO Warrenty and NO Liability.\n"
          "Info on license of used 3rd party libraries included in separate\n"
          "license.txt or refer to info below.\n\n"
          "MORE INFORMATION, source code and contact info are\n"
          "available at GitHub: https://github.com/nofaceinbook/betterflat/\n\n"
          "Hope the tool helps you.    (c) 2020 by schmax (Max Schmidt)"   
          ).grid(row=0, pady=10, padx=10)

def ConfigMenu(win, config):
    def update_config():
        config.accuracy = float(accuracy_entry.get())
        config.isprofile = profileType.get()
        config.text = text_entry.get("1.0",END)
        config.terrain = terrain_entry.get()
        config.isTerrainReplaced = replaceTerrainType.get()
        config.interval = float(interval_entry.get())
    configwin = Toplevel(win)
    toplabel = Label(configwin, anchor=W, justify=LEFT, text="Change configurations for cut option:").grid(row=0, column=0, columnspan=2, pady=10, padx=10)
    accuracy_label = Label(configwin, text="Accuracy (in m):")
    accuracy_label.grid(row=1, column=0, pady=4, sticky=E)
    accuracy_entry = Entry(configwin, width=7)
    accuracy_entry.grid(row=1, column=1, sticky=W)
    accuracy_entry.insert(0, config.accuracy)
    terrain_label = Label(configwin, text="TerrainFile:")
    terrain_label.grid(row=2, column=0, pady=4, sticky=E)
    terrain_entry = Entry(configwin, width=40)
    terrain_entry.grid(row=2, column=1, sticky=W)
    terrain_entry.insert(0, config.terrain)
    replaceTerrainType = IntVar() # 1 if terrain should be used to replace existing ones, 0 if not
    replaceTerrainType.set(config.isTerrainReplaced)
    replaceTerrainCB = Checkbutton(configwin, text="Replace terrain under boundary/rwy", variable=replaceTerrainType)
    replaceTerrainCB.grid(row=3, column=0, sticky=E, pady=4)
    profileType = IntVar() # 1 if profile in rwys should be cut, 0 if not
    profileType.set(config.isprofile)
    profileCB = Checkbutton(configwin, text="Cut Runway Profile (always uses terrain above)", variable=profileType)
    profileCB.grid(row=4, column=0, sticky=E, pady=4)
    text_label = Label(configwin, anchor=W, justify=LEFT, text="Profile Definition (if empty raster or tria elevation of dsf is used):").grid(row=5, column=0, columnspan=2, pady=4, padx=10)
    text_entry = Text(configwin, height=2, width=60)
    text_entry.grid(row=6, column=0, columnspan=2, padx=10, pady=4)
    text_entry.insert(END, config.text)
    interval_label = Label(configwin, text="Profile intervals (in m):")
    interval_label.grid(row=9, column=0, pady=4, sticky=E)
    interval_entry = Entry(configwin, width=7)
    interval_entry.grid(row=9, column=1, sticky=W)
    interval_entry.insert(0, config.interval)
    update_button = Button(configwin, text='Update Configuration', command=update_config)
    update_button.grid(row=10, column=0, pady=4)
    buttonlabel = Label(configwin, anchor=W, justify=LEFT, text="Changes are only applied when you press update button.\nClose window when done.").grid(row=11, column=0, columnspan=2, pady=10, padx=10)

def defineLog(logname, logLevelStream='INFO', logLevelFile='INFO', mode='w+'):
    #
    # defines log-Handler, if Level is set to None no File/Stream handler will be created
    # The file for the stream is written to logname.log
    # Global variable LogName is defined to be used e.g. to create sub-loggers
    # Returns the created logger log to be used
    #
    global LogName
    LogName = logname
    directory = path.dirname(path.abspath(__file__))
    logfile = path.join(directory, logname + '.log')
    log = getLogger(logname)
    log.handlers = []  # Spyder/IPython currently does not remove existing loggers, this does the job
    if logLevelStream == 'DEBUG' or logLevelFile == 'DEBUG':
        log.setLevel('DEBUG') #so maximum level is DEBUG; needed to see with getEffectiveLevel if DEBUG is enabled or not  
    else:
        log.setLevel('INFO') #set maximum level if not DEBUG before restriction by handlers below
    formatter = Formatter('%(name)s-%(levelname)s: %(message)s')
    if logLevelStream:
        stream_handler = StreamHandler()
        stream_handler.setLevel(logLevelStream)
        stream_handler.setFormatter(formatter)
        log.addHandler(stream_handler)
    if logLevelFile:
        file_handler = FileHandler(logfile, mode)
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
    # Only land runways (type 100) are returned in a list with sub-lists of their endpoints, width in meter
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
    if not path.isfile(filename):
        log.error("Airport File {} does not exist!".format(filename))
        return [], "Error: Airport file does not exist!"
    with open(filename, encoding="utf8", errors="ignore") as f:
        for line in f:
            v = line.split()
            if len(v) == 0:  # don't consider empty lines
                continue
            if len(v) > 4:  # check if correct airport section in file is reached
                if v[0] == '1' or v[0] == '16' or v[0] == '17':
                    if v[4] == icao_id or icao_id =='': #if no icao id is given just first airport is selected
                        Airport = True
                        icao_id = v[4] #set now icao id in case it was '' before
                        apt_elev = round(int(v[1]) * 0.3048)
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
                        bounds[-1].append([float(v[2]), float(v[1])])  # Note: Bezier definitins are not considered, just the base point
                    elif v[0] == '113' or v[0] == '114':
                        bounds[-1].append([float(v[2]), float(v[1])])  # Note: Bezier definitins are not considered, just the base point
                        bounds[-1].append(bounds[-1][0])  # #form closed loop by adding again first vertex
                        BoundarySection = False
                        log.info("Boundary no. {} with {} vertices read.".format(len(bounds), len(bounds[-1])))
    if len(bounds) == 0 and len(runways) == 0:
        log.warning("No valid boundary or runway found in file!")
        return [], [], None, None, None, "Warning: No valid boundary or runway found in file!"
    else:
        log.info("Finished reading boundaries.")
        return bounds, runways, apt_elev, apt_flatten, apt_name, None


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
        m = -1 / ((p2[0] - p1[0]) / (p2[1] - p1[1])) #gradient of perpendicular line
        dx = sqrt( ( (w/2)**2) / (1 + m**2)) 
        dy = dx * m 
    dx /= degree_dist_at_lat #convert meters in longitute coordinate difference at geographical latitude
    dy /= lat_degree_dist_average #convert meters in latitude coordinate difference 
    l = []
    if (p1[1] <= p2[1] and dy >= 0) or (p1[1] > p2[1] and dy < 0): #make sure to always insert in clockwise order
        l.append([round(p1[1] - dx, 8), round(p1[0] - dy, 8)]) #buttom corner1
        l.append([round(p1[1] + dx, 8), round(p1[0] + dy, 8)]) #buttom corner2
        l.append([round(p2[1] + dx, 8), round(p2[0] + dy, 8)]) #top corner1
        l.append([round(p2[1] - dx, 8), round(p2[0] - dy, 8)]) #top corner2
    else: #insert vertices in different order to assure clockwise orientation
        l.append([round(p1[1] + dx, 8), round(p1[0] + dy, 8)])
        l.append([round(p1[1] - dx, 8), round(p1[0] - dy, 8)])
        l.append([round(p2[1] - dx, 8), round(p2[0] - dy, 8)])
        l.append([round(p2[1] + dx, 8), round(p2[0] + dy, 8)])
    l.append(l[0]) #add first corner to form closed loop
    return l

def gauss_jordan(m, eps = 1.0/(10**10)):
    """Puts given matrix (2D array) into the Reduced Row Echelon Form.
       Returns True if successful, False if 'm' is singular.
       NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
       Written by Jarno Elonen in April 2005, released into Public Domain"""
    (h, w) = (len(m), len(m[0]))
    for y in range(0,h):
      maxrow = y
      for y2 in range(y+1, h):    # Find max pivot
        if abs(m[y2][y]) > abs(m[maxrow][y]):
          maxrow = y2
      (m[y], m[maxrow]) = (m[maxrow], m[y])
      if abs(m[y][y]) <= eps:     # Singular?
        return False
      for y2 in range(y+1, h):    # Eliminate column y
        c = m[y2][y] / m[y][y]
        for x in range(y, w):
          m[y2][x] -= m[y][x] * c
    for y in range(h-1, 0-1, -1): # Backsubstitute
      c  = m[y][y]
      for y2 in range(0,y):
        for x in range(w-1, y-1, -1):
          m[y2][x] -=  m[y][x] * m[y2][y] / c
      m[y][y] /= c
      for x in range(h, w):       # Normalize row y
        m[y][x] /= c
    return True

def lin_equation_solve(M, b):
    """
    solves M*x = b
    return vector x so that M*x = b
    :param M: a matrix in the form of a list of list
    :param b: a vector in the form of a simple list of scalars
    """
    m2 = [row[:]+[right] for row,right in zip(M,b) ]
    result = gauss_jordan(m2)
    return [row[-1] for row in m2] if result else None

def getspline(xp, yp):
    """
    for x values xp with according y values yp a natural cubic spline is defined
    note: x values should be sorted from lowest to highest value!
    note: generates only float values to be compatible with gaus_jordan function used
    returns list with deepcopy list of x-values and list of cubic spline paramters (one sgement after the other)
    """
    points = len(xp)
    segments = points - 1
    if (points != len(yp)) or (points < 3):
        return None
    A = []
    b = []
    xp_returned = [] #deepcopy of returned x-values
    for i in range(points):
        xp_returned.append(float(xp[i]))
    for i in range(4 * segments):
        A.append([])
        b.append(0.0)
        for j in range(4 * segments):
            A[i].append(0.0)
    for i in range(segments):
        #condition for left end of segment
        A[i][4*i+0] = float(xp[i]**3)
        A[i][4*i+1] = float(xp[i]**2)
        A[i][4*i+2] = float(xp[i])
        A[i][4*i+3] = 1.0
        b[i] = float(yp[i])
        #condition for right end of segment
        A[segments+i][4*i+0] = float(xp[i+1]**3)
        A[segments+i][4*i+1] = float(xp[i+1]**2)
        A[segments+i][4*i+2] = float(xp[i+1])
        A[segments+i][4*i+3] = 1.0        
        b[segments+i] = float(yp[i+1])
        if i == 0:
            continue #do outer points later, so one row missing now therefore -1 below
        #condition for first derivation of inner points, setting b value to 0 omitted
        A[2*segments+i-1][4*(i-1)+0] = float(3*xp[i]**2)
        A[2*segments+i-1][4*(i-1)+1] = float(2*xp[i])
        A[2*segments+i-1][4*(i-1)+2] = 1.0
        A[2*segments+i-1][4*(i-1)+4] = float(-3*xp[i]**2)
        A[2*segments+i-1][4*(i-1)+5] = float(-2*xp[i])
        A[2*segments+i-1][4*(i-1)+6] = -1.0
        #condition for second derivation of inner points, setting b value to 0 omitted
        A[3*segments+i-1][4*(i-1)+0] = float(6*xp[i])
        A[3*segments+i-1][4*(i-1)+1] = 2.0
        A[3*segments+i-1][4*(i-1)+4] = float(-6*xp[i])
        A[3*segments+i-1][4*(i-1)+5] = -2.0
    # Now consider derivation for endpoints, here for NATURAL SPLINE so second derivation 0 at ends, setting b to 0 omitted
    A[3*segments-1][0] = float(6*xp[0])
    A[3*segments-1][1] = 2.0
    A[4*segments-1][4*(segments-1)] = float(6*xp[segments])
    A[4*segments-1][4*(segments-1)+1] = 2.0
    #Now solve according linear equations
    x = lin_equation_solve(A, b)
    return [xp_returned, x]

def evalspline(x, spline): #evaluates spline at position x
    #assumes spline in format with two lists, first x-values giving intervals/segments, second all paramters for cubic spline one after the other
    #important: it is assumed that x-values for intervals/segments are ordered from low to high
    for i in range(len(spline[0]) - 1):#splines are one less then x-points
        if x <= spline[0][i+1]: #for all points before second point use first spline segment, for all after the last-1 use last spline segment
            break
    return spline[1][4*i+0] * x**3 + spline[1][4*i+1] * x**2 + spline[1][4*i+2] * x + spline[1][4*i+3] 
         
def interpolateRWYprofile(rwys, dsf, rwyNum, defintion=""): #calculates the interpolated rwy profiles for runwy number rwyNum       
    #definition could include profile defined in "0@98.2 180@97.4 ...." where value before @ is distance and after elevation, use ";" to seperate multiple rwy 
    log.info("Interpolating runway profile for runway number {}".format(rwyNum))
    logprofile = "" #profile as text-definiton in format as in input field   
    rwy = rwys[rwyNum] # runway is selected from severl ones 
    start = (rwy[0][1], rwy[0][0]) #start coordinates rwy
    end = (rwy[1][1], rwy[1][0]) #end coordinates rwy
    l = distance(start, end) #length of runway
    x_points = []
    y_points = []
    if not "@" in defintion: #seems that defintion includes no elevation values in required format so use elevation given by raster in dsf
        log.info("Based on raster elevation of dsf the following profile was generated (distance-from-rwy-start@elevation in m):")
        for d in range(-100, int(l)+101, 100): #starting before and ending after rwy #### TBD: variable width for definition points
            p = (start[0] + d/l * (end[0] - start[0]), start[1] + d/l * (end[1] - start[1]))
            x_points.append(d)
            y_points.append(dsf.getElevation(p[0], p[1]))
            logprofile += "{}@{} ".format(d, round(y_points[-1], 2))
    else: #use defintion for setting profile
        runwayValues = defintion.split(";")
        if len(runwayValues) - 1 >  rwyNum:
            log.error("Definition of profile is not for all runways!! Please enter profiles for all runway seperated with semicolon!!")
            return None ###### TBD: User-friendly output in GUI
        log.info("Based on input text-defintion following profile is created (distance-from-rwy-start@elevation in m):")
        logprofile = runwayValues
        vpairs = runwayValues[rwyNum].split()
        for v in vpairs:
            vx, vy = v.split("@")
            x_points.append(float(vx))
            y_points.append(float(vy))
    #rwySpline = interpolate.splrep(x_points, y_points, s=12) ######### TBD: smooth factor just fixed test value, has to be variable !!!! ###########
    rwySpline = getspline(x_points, y_points)
    log.info(logprofile)
    return rwySpline

def interpolatedRWYelevation(rwy, p, rwySpline): #based on runway it's spline profile, the elevation of a point on the runway is calculated       
    start = (rwy[0][1], rwy[0][0]) #start coordinates rwy 
    end = (rwy[1][1], rwy[1][0]) #end coordinates rwy
    startD = (start[0] - 0.1 * (end[0] - start[0]), start[1] - 0.1 * (end[1] - start[1])) #use starting point 10% of rwy length before to really get value for all points around runway
    endD = (start[0] + 1.1 * (end[0] - start[0]), start[1] + 1.1 * (end[1] - start[1]))   #use end point 10% of rwy length behind to really get value for all points around runway
    inclination_of_ortho = (end[1] - start[1], start[0] - end[0])
    orthoStartD = (p[0] - inclination_of_ortho[0], p[1] - inclination_of_ortho[1]) # Start of orthogonal line of RWY through point p with length double of RWY (to guarentee intersection on center line)
    orthoEndD = (p[0] + inclination_of_ortho[0], p[1] + inclination_of_ortho[1]) # End of orthogonal line of RWY through point p with length double of RWY (to guarentee intersection on center line)
    p_centered = intersection(startD, endD, orthoStartD, orthoEndD) #location of p on center line
    d = distance(start, p_centered)
    #elev = interpolate.splev(d, rwySpline)
    elev = evalspline(d, rwySpline)
    return elev
    #return elev.item()        


def moveSubmPools(dsf, V):
    log.info("Moving vertices to pools allowing submeter definiton and remove vertices in old pool as fas that other indeces are not affected.")
    newV = {} #dictionary mapping old vertices in V to the ones in the new submPool
    Vcoords = [] #list of lon,lat coordinates of all vertices in V
    oldPoolIDs = set() #set of pool ids from pools where vertices moved from; needed to delete not required vertices from these pools
    for v in V:
        oldPoolIDs.add(v[0])
        newSubmV = []
        ScaleForNewSubmV = []
        for i in range(len(dsf.V[v[0]][v[1]])): #This is doing deepcopy for vertex and scaling     
            newSubmV.append(dsf.V[v[0]][v[1]][i])
            ScaleForNewSubmV.append([dsf.Scalings[v[0]][i][0], dsf.Scalings[v[0]][i][1]])
        new_pid, new_vindex = dsf._addSubmVertex_(newSubmV, ScaleForNewSubmV)    
        newV[(v[0], v[1])] = (new_pid, new_vindex)
        log.debug("Moving vertex with index {} from pool {} to pool {} with new index {}.".format(v[1], v[0], new_pid, new_vindex))
        Vcoords.append((dsf.V[v[0]][v[1]][0], dsf.V[v[0]][v[1]][1]))
    for p in dsf.PatchesInArea(*dsf.BoundingRectangle(Vcoords)): #now update all trias in area of V with reference to new Pools
        for t in p.trias: #go through all trias that are in patches that have coordinates in area of Vcoords
            for i in range(3):
                if (t[i][0], t[i][1]) in V:  #if one vertex of tria is in the list of vertices moved to submPools
                    t[i] = newV[(t[i][0], t[i][1])] #change to updated to vertex
        p.trias2cmds()
    for i in oldPoolIDs:
        while (i, len(dsf.V[i]) - 1) in V:
            x = dsf.V[i].pop()
            log.debug("Removed vertex {} with index {} from pool {}.".format(x, len(dsf.V[i]), i))
    return newV.values() #list of vertices with updated poolIDs is returned 


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
                    if intersection(TriaV[i], TriaV[i + 1], poly[j], poly[j + 1]):  # current triangle intersects with an poly line
                        s.add((t[0][0], t[0][1]))  # add all 3 vertices of tria to set as tuples
                        s.add((t[1][0], t[1][1]))
                        s.add((t[2][0], t[2][1]))
            if PointInPoly(TriaV[0],poly):  # Check also that not complete Tria lies in poly by checking for first vertex of Tria
                s.add((t[0][0], t[0][1]))  # add all 3 vertices of tria to set as tuples
                s.add((t[1][0], t[1][1]))
                s.add((t[2][0], t[2][1]))
            if PointInPoly(poly[0], TriaV):  # Check also that not complete poly lies in current tria by checking for first vertex of poly
                s.add((t[0][0], t[0][1]))  # add all 3 vertices of tria to set as tuples
                s.add((t[1][0], t[1][1]))
                s.add((t[2][0], t[2][1]))
    log.info("{} vertices of mesh found that belong to mesh triangles intersecting or within boundary.".format(len(s)))
    return getAllVerticesForCoords(dsf, s)


def getAllVerticesForCoords(dsf, s):  ########### TBD: Call this function from function above
    ### For all vertices in set s that need to be changed (e.g. new elevation), the lon/lat coordinates are taken and all other vertices on these coordinates are determined and returned
    changedVertices = {}  # set up dictonary with all coords of changed vertices to find additional vertices at such coords for trias out of bound
    for v in s:
        changedVertices[(round(dsf.V[v[0]][v[1]][0], 7), round(dsf.V[v[0]][v[1]][1], 7))] = True  # coords round to range of cm
    miny, maxy, minx, maxx = dsf.BoundingRectangle(changedVertices)
    delta = 0.000001  # check also for vertices 0.1m outside boundary
    for p in dsf.PatchesInArea(miny - delta, maxy + delta, minx - delta, maxx + delta):
        for t in p.triangles():
            for v in t:
                if (round(dsf.V[v[0]][v[1]][0], 7), round(dsf.V[v[0]][v[1]][1], 7)) in changedVertices:
                    s.add((v[0], v[1]))  # add tuple of vertex at coords to make sure to get all vertices at the coords
    log.info("References to vertices to be changed: {}".format(s))
    log.info("Unique coords to be changed: {}".format(changedVertices))
    return s, changedVertices    #returns set of all vertices and coordinates and the coordinates



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

class ConfigDetails:
    def __init__(self):
        self.accuracy = 0.1  #value in m used for cutting, when existing vertices are within that distance they will be reused
        self.isprofile = 0 #is set to 1 in case the runways will be cut in segements in order to allow profile, 0 for no profile
        self.text = "" #input text to be used e.g. for definition of runway profile
        self.terrain = "lib/g10/terrain10/grass_tmp_sdry_fl.ter" #value used when terrain of mesh trias is replaced, e.g. under runway profile
        self.isTerrainReplaced = 0 #is set to 1 in case terrain under boundary or runway should be replaced with terrain definition (in case of isprofile = 1 this is always done regardless of this defintion)
        self.interval = 25 #distance of intervals in m in which runway is cut to generate profile

class bflatGUI:
    def __init__(self):

        self.dsf = XPLNEDSF(LogName, self.showProgress)  # includes all information about the read dsf-file
        self.boundaries = []  # set of bouundaries where each boundary includes vertices of boundary
        self.runways = [] #list of runways for processed airport (endpoints, width in m)
        self.apt_elev = None #airport elevation as defined in apt.dat in m
        self.vertices = set([])  # set of vertices that should get new height
        self.dsfreadfile = ""  # name of the file that is stored in self.dsf
        self.current_action = None #is set to 'read' or 'write' when reading/writing dsf
        self.config = ConfigDetails() #additional values will be set, especially for configure cutting of mesh

        self.window = Tk()
        self.window.title("X-Plane bflat (version: {})".format(bflat_VERSION))

        self.header = Label(self.window, text="Make sure you always have copies of your original files!")
        self.header.grid(row=0, column=0, columnspan=2)
        self.config_button = Button(self.window, text=' Config ', fg='black', command=lambda: ConfigMenu(self.window, self.config))
        self.config_button.grid(row=0, column=2, sticky=W, pady=4, padx=10)        
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
        self.cuttype = IntVar() # 1 if mesh should be cut, 0 if just elevation of underlaying trias should be adapted
        self.cuttype_checkB = Checkbutton(self.window, text="cut", variable=self.cuttype)
        self.cuttype_checkB.grid(row=3, column=3)
        self.cuttype_checkB.select()
        
        self.apt_read = Button(self.window, text='Read Boundary', command=lambda: self.read_apt(self.apt_entry.get(), self.aptid_entry.get()))
        self.apt_read.grid(row=4, column=0, sticky=E, pady=4)
        self.apt_status_label = Label(self.window, text="ICOA Id only required if file contains several airports!")
        self.apt_status_label.grid(row=4, column=1, columnspan=2, sticky=W)

        self.dsffile_label = Label(self.window, text="DSF File:")
        self.dsffile_label.grid(row=5, column=0, sticky=W)
        self.dsf_entry = Entry(self.window, width=60)
        self.dsf_entry.grid(row=5, column=1, columnspan=2, sticky=W)
        self.dsf_select = Button(self.window, text='Select', command=lambda: self.select_file(self.dsf_entry)).grid(row=5, column=3, sticky=W, pady=4, padx=10)

        self.dsf_read = Button(self.window, text='  Read DSF   ', state=DISABLED, command=lambda: self.read_dsf(self.dsf_entry.get()))
        self.dsf_read.grid(row=6, column=0, sticky=E, pady=4)
        self.dsf_status_label = Label(self.window, text="   this can take some minutes....")
        self.dsf_status_label.grid(row=6, column=1, sticky=W)

        self.result_text = Label(self.window, text="Result:")
        self.result_text.grid(row=7, column=0, sticky=E, pady=8)
        self.result_label = Label(self.window, text=" not yet")
        self.result_label.grid(row=7, column=1, sticky=W, pady=8)
        self.kml_create = Button(self.window, text='.kml', state=DISABLED,  command=lambda: self.write_kml())
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

        self.write_button = Button(self.window, text="Write DSF", state=DISABLED, command=lambda: self.write_dsf(self.height_entry.get(), self.newdsf_entry.get(), self.bakdsf_entry.get()))
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
            if flatten_flag != None:
                self.removeFlattenFromAPT(filename, icao, flatten_flag)
    
    def removeFlattenFromAPT(self, filename, icao, flatten_flag):
        flattenwin = Toplevel(self.window)
        Label(flattenwin, anchor=W, justify=LEFT, text=
          "WARNING\n\n"
          "The airport definition {} in the apt-file {}\n"
          "includes a definition for flattening '1302 flatten {}'.\n"
          "Better flattening results of this tool will not be shown with this\n"
          "flag set in the airport defintion.\n"
          "You can now directly remove this line from the file.\n"
          "Before doing so you should make a backup of apt-file if not already done.\n"
          "ATTENTION: Existing files will be overwritten\n\n"
          "In case you finished or you don't want to apply changes just close this window.".format(icao, filename, flatten_flag)
          ).grid(row=0, column=0, columnspan=3, pady=10, padx=10)
        Label(flattenwin, text="Name for Backupfile:").grid(row=1, column=0, sticky=W)
        apt_bak_entry = Entry(flattenwin, width=60)
        apt_bak_entry.grid(row=1, column=1, columnspan=2, sticky=W)
        apt_bak_entry.delete(0, END)
        apt_bak_entry.insert(0, filename + '.bak')
        Button(flattenwin, text='Select', command=lambda: self.select_file(apt_bak_entry)).grid(row=1, column=3, sticky=W, pady=4, padx=10)
        bak_status_label = Label(flattenwin, text="")
        bak_status_label.grid(row=2, column=2, sticky=W, pady=4, padx=10)
        Button(flattenwin, text='Create Backup', command=lambda: self.write_apt_backup(filename, apt_bak_entry.get(), bak_status_label)).grid(row=2, column=1, sticky=W, pady=4, padx=10)
        Label(flattenwin, text="Name of updated apt-file:").grid(row=3, column=0, sticky=W)
        apt_update_entry = Entry(flattenwin, width=60)
        apt_update_entry.grid(row=3, column=1, columnspan=2, sticky=W)
        apt_update_entry.delete(0, END)
        apt_update_entry.insert(0, filename)
        Button(flattenwin, text='Select', command=lambda: self.select_file(apt_update_entry)).grid(row=3, column=3, sticky=W, pady=4, padx=10)
        update_status_label = Label(flattenwin, text="(will remove the 1302 flatten definiton)")
        update_status_label.grid(row=4, column=2, sticky=W)
        Button(flattenwin, text='Update apt-file', command=lambda: self.update_apt(apt_update_entry.get(), icao, update_status_label)).grid(row=4, column=1, sticky=W, pady=4, padx=10)

        
    def write_apt_backup(self, orgfile, bakfile, statuslabel):
        log.info("Copying original apt-file: {} to backup-file: {}".format(orgfile, bakfile))
        try:
            copy2(orgfile, bakfile)
            statuslabel.config(text="Backup created successfully.")
        except:
            log.error('{} can not be replaced!'.format(orgfile))
            statuslabel.config(text="Error: Original file {} can not be replaced!".format(orgfile))
    
    def update_apt(self, filename, icao_id, statuslabel):
        log.info("Updating airport data from: {}".format(filename))
        Airport = False  # else first the correct icoa id has to be found before Airport becomes true
        apt_name = None
        try: 
            with open(filename, "r", encoding="utf8") as f:
                lines = f.readlines()           
        except:
            log.error("apt-file {} not readable!".format(filename))
            statuslabel.config(text="Error: apt-file {} not readable!".format(filename))
            return
        with open(filename, "w", encoding="utf8") as f:  
            for line in lines:
                v = line.split()
                #if len(v) == 0:  # don't consider empty lines
                #    continue
                if len(v) > 4:  # check if correct airport section in file is reached
                    if v[0] == '1' or v[0] == '16' or v[0] == '17':
                        if v[4] == icao_id or icao_id =='': #if no icao id is given just first airport is selected
                            Airport = True
                            icao_id = v[4] #set now icao id in case it was '' before
                            apt_name = v[5] 
                            if len(v) > 6: apt_name = apt_name + " " + v[6]
                            log.info("Airport {} found where flattened shall be removed.".format(apt_name))
                        else:
                            Airport = False  # change to false in case of new different airport
                if len(v) > 2 and Airport and v[0] == '1302' and v[1] == 'flatten':
                        log.info("Line with flatten flag skipped for writing and thus removed.")
                else:
                    f.write(line)
        statuslabel.config(text="apt-file updated successfully for {}".format(apt_name))       

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
            self.vertices = set([])  #set to empty in case of previous computations
            coords = {}
            if self.boundtype.get() == 'airport':
                bounds = self.boundaries
                log.info("Using airport boundaries for calculating area to be flattened.")
            else:
                bounds = []
                for r in self.runways:
                    bounds.append(getRunwayBounds(r[0], r[1], r[2] + 2)) #add 2m width for shoulders
                log.info("Using runway boundaries for calculating area to be flattened.")
            boundNumber = -1 #index number for courent boundary processed below    ##### NEW 5 ####
            for boundary in bounds:  # get intersections for all boundaries
                boundNumber += 1
                self.result_label.config(text="Calculate intersections for boundary number {}. Wait....".format(boundNumber))
                self.window.update()
                if self.cuttype.get():
                    if not self.config.isprofile: #just cut without inserting segments for runway profile
                        if self.config.isTerrainReplaced:
                            verticesFromCut = self.dsf.cutTerrainInMesh(boundary, self.config.terrain, self.config.accuracy) #
                        else:
                            verticesFromCut = self.dsf.cutPolyInMesh(boundary, self.config.accuracy) 
                        vs, cs = getAllVerticesForCoords(self.dsf, verticesFromCut)  
                        self.vertices = self.vertices.union(vs)
                        coords.update(cs)                        
                    else: #cut also runways in segments to allow sloped profile    
                        RWYvertices = {} #set of vertices for the current runway  
                        newRWYvertics = [] #list of the new vertices after been inserted in submeter Pool   
                        RWYvertices = self.dsf.cutRwyProfileInMesh(boundary, self.config.terrain, self.config.interval, self.config.accuracy)                     
                        vs, cs = getAllVerticesForCoords(self.dsf, RWYvertices)
                        coords.update(cs)
                        rwySpline = interpolateRWYprofile(self.runways, self.dsf, boundNumber, self.config.text) #tbd: update for several runways ######
                        for v in RWYvertices: #updates heigt of vertices based on spline profile
                            newelev = interpolatedRWYelevation(self.runways[boundNumber], (self.dsf.V[v[0]][v[1]][0], self.dsf.V[v[0]][v[1]][1]), rwySpline)  
                            newelev = round(newelev, 3)
                            dist = round(distance((self.dsf.V[v[0]][v[1]][0], self.dsf.V[v[0]][v[1]][1]), (self.runways[boundNumber][0][1], self.runways[boundNumber][0][0]))) #### CHANGED from runway 0 to boundNumber !! ########
                            log.debug("Vertex {} with coordinates {} at approx distance {}m set to elevation {}m.".format(v, (self.dsf.V[v[0]][v[1]][0], self.dsf.V[v[0]][v[1]][1]), dist, newelev))
                            self.dsf.V[v[0]][v[1]][2] = newelev
                        for v in RWYvertices:
                            log.debug("Old vertex {} is: {}".format(v, self.dsf.V[v[0]][v[1]]))
                        newRWYvertices = moveSubmPools(self.dsf, RWYvertices) #moves vertices to new pool allowing submeter elevation #### WARNING: Returns list of dictionary values, no set !!!!!
                        for v in newRWYvertices:
                            self.vertices.add(v)  #add new vertices list element to self.vertices set for further processing  
                            log.debug("New vertex {} is: {}".format(v, self.dsf.V[v[0]][v[1]]))
                else:
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

    def write_kml(self):
        if self.boundtype.get() == 'airport':
            log.info("Writing kml data for airport boundary to file {}.".format(self.dsf_entry.get() + '.kml'))
            kmlExport(self.dsf, self.boundaries, self.vertices, self.dsf_entry.get())
        else:
            bounds = []
            for r in self.runways:
                bounds.append(getRunwayBounds(r[0], r[1], r[2] + 2)) #add 2m width for shoulders
            log.info("Writing kml data for runway boundary to file {}.".format(self.dsf_entry.get() + '.kml'))
            kmlExport(self.dsf, bounds, self.vertices, self.dsf_entry.get())

    def write_dsf(self, newheight, dsffile, bakfile):
        self.current_action = 'write'
        newheight = int(self.height_entry.get())
        if bakfile != "":
            log.info("Moving original dsf file: {} to backup-file: {}".format(dsffile, bakfile))
            try:
                replace(self.dsfreadfile, bakfile)
            except:
                log.error('{} can not be replaced!'.format(self.dsfreadfile))
                self.write_status_label.config(text="Error: Original file {} can not be replaced!".format(self.dsfreadfile))
                self.current_action = None
                return
        if not (self.cuttype.get() and self.config.isprofile): #no profile is cut, so now the height of vertices needs to be adapted.
            log.info("Writing dsf file: {} for height: {}".format(dsffile, newheight))
            for v in self.vertices:
                log.debug(" {} with height {} changed to: {}".format(self.dsf.V[v[0]][v[1]], self.dsf.getElevation(self.dsf.V[v[0]][v[1]][0], self.dsf.V[v[0]][v[1]][1], self.dsf.V[v[0]][v[1]][2]), newheight))
                self.dsf.V[v[0]][v[1]][2] = newheight
        else: #in case of runway profile height is already set for the new vertices
            log.info("Writing dsf file: {} with runway profile.".format(dsffile))
        self.write_status_label.config(text="Writing changes ...") 
        self.window.update()
        self.dsf.write(dsffile)
        self.write_status_label.config(text="Done.")
        log.info("Done.")
        self.current_action = None

    def select_file(self, entry):
        filename = askopenfilename()
        entry.delete(0, END)
        entry.insert(0, filename)
        entry.focus()



log = defineLog('bflat', None, 'INFO') #no log on console for EXE version
log.info("Started bflat Version: {}".format(bflat_VERSION))
main = bflatGUI()
