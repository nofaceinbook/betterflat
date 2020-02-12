# -*- coding: utf-8 -*-
#******************************************************************************
#
# mexp.py
#        
mexp_VERSION = "0.0.1 exp"
# ---------------------------------------------------------
# Python Tool for X-Plane Mesh Editing
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
from xplnedsf2 import *
from muxpKMLexport import *
from muxpKMLexport2 import *
from os import path, replace
from shutil import copy2
from math import sin, cos, sqrt, radians #for getRunwayBounds calculations

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


def extractMeshArea(dsf, latS, latN, lonW, lonE):
    """
    Extracts an area [latS, latN, lonW, lonE] from mesh in dsf and returns it as list of tria lists
    in the following form where for each tria the following lists are included:
    [0 - 2] list of all coordinates of tria vertices including s/t ##### ONLY REFERENCE NO DEEPCOPY / BETTER JUST DEEPCOPY of x,y corrds ?????
    [3 - 5] list of pool and vertex id in pool in dsf
    [6] index to patch tria was in dsf
    Important: This function really removes trias from dsf
    """
    triaCount = 0 #counts all trias in dsf
    extracts = [] #list of extracted trias
    for p in dsf.Patches:
        trias = p.triangles()
        tremoved = [] #trias that should be removed
        for t in trias:
            # get for each tria the bounding rectangle in miny, maxy, minx, maxxx
            minx = min(dsf.V[t[0][0]][t[0][1]][0], dsf.V[t[1][0]][t[1][1]][0], dsf.V[t[2][0]][t[2][1]][0])
            maxx = max(dsf.V[t[0][0]][t[0][1]][0], dsf.V[t[1][0]][t[1][1]][0], dsf.V[t[2][0]][t[2][1]][0])
            miny = min(dsf.V[t[0][0]][t[0][1]][1], dsf.V[t[1][0]][t[1][1]][1], dsf.V[t[2][0]][t[2][1]][1])
            maxy = max(dsf.V[t[0][0]][t[0][1]][1], dsf.V[t[1][0]][t[1][1]][1], dsf.V[t[2][0]][t[2][1]][1])
            #now check if bounding rectangle of tria intersects with area
            if not (minx < lonW and maxx < lonW): #x-range of box is not completeley West of area     
                if not (minx > lonE and maxx > lonE): #x-range of box is not completele East of area
                    if not (miny < latS and maxy < latS): #y-range is not completele South of area
                        if not (miny > latN and maxy > latN): #y-range is not conmpletele North of ares
                            extracts.append([ dsf.V[t[0][0]][t[0][1]], dsf.V[t[1][0]][t[1][1]], dsf.V[t[2][0]][t[2][1]] ])
                            extracts[-1].extend(t)  #so we have an intersection of tria box with area and append the tria
                            extracts[-1].append(dsf.Patches.index(p))
                            tremoved.append(t)
                            log.info("Tria {} with latS: {}  latN: {}  lonW: {} and lonE: {} in area.".format(t, miny, maxy, minx, maxx))
            triaCount += 1
        for t in tremoved:
            trias.remove(t)
        p.trias2cmds(trias) #updates dsf patch with trias not including removed ones
    log.info("dsf has {} trias and {} have been extracted!".format(triaCount, len(extracts)))
    return extracts


def _linsolve_(a1, b1, c1, a2, b2, c2):  
    divisor = (a1 * b2) - (a2 * b1)
    if divisor == 0:  ## WARNING: for colineear setting special returns special value no error!!!!!
        return -99999, -99999  ### points are colinear, might intersect or not BUT here with negative values calling intersection function returns None
    return round(((c1 * b2) - (c2 * b1)) / divisor, 8), round(((a1 * c2) - (a2 * c1)) / divisor, 8)  # ROUNDING TO ALWAYS GET SAME CONCLUSION

def intersection(p1, p2, p3, p4):  # checks if segment from p1 to p2 intersects segement from p3 to p4   ### NEW - taken from bflat ###
    s0, t0 = _linsolve_(p2[0] - p1[0], p3[0] - p4[0], p3[0] - p1[0], p2[1] - p1[1], p3[1] - p4[1], p3[1] - p1[1])
    if s0 >= 0 and s0 <= 1 and t0 >= 0 and t0 <= 1:
        return (round((p1[0] + s0 * (p2[0] - p1[0])), 8), round(p1[1] + s0 * (p2[1] - p1[1]), 8))  ### returns the cutting point as tuple; ROUNDING TO ALWAYS GET SAME POINT 
    else:                    
        return (None)

def distance(p1, p2): #calculates distance between p1 and p2 in meteres where p is pair of longitude, latitude values          
    R = 6371009 #mean radius earth in m
    lon1 = radians(p1[0]) # WARNING: changed order as it is in bflat now !!!!!!!!
    lat1 = radians(p1[1])
    lon2 = radians(p2[0])
    lat2 = radians(p2[1]) 
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))   
    return R * c

def PointLocationInTria(p, t): #delivers location of point p in Tria t by vectors spanned by t from last vertex in tria
    denom = ((t[1][1] - t[2][1])*(t[0][0] - t[2][0]) + (t[2][0] - t[1][0])*(t[0][1] - t[2][1]))
    if denom == 0: ### to be checked when this is the case!!!!
        return -0.01, -0.01 ###NEW NEW NEW, should actually never be the case, but for isPointInTria it delviers NO!
    nom_a = ((t[1][1] - t[2][1])*(p[0] - t[2][0]) + (t[2][0] - t[1][0])*(p[1] - t[2][1]))
    nom_b = ((t[2][1] - t[0][1])*(p[0] - t[2][0]) + (t[0][0] - t[2][0])*(p[1] - t[2][1]))
    a = nom_a / denom
    b = nom_b / denom
    return a, b #returns multiplier for vector (t2 - t0) and for vector (t2 - t1) starting from point t2

def isPointInTria(p, t): #delivers True if p lies in t, else False
    a, b = PointLocationInTria(p, t)
    c = 1 - a - b
    return (0 <= a <= 1 and 0 <= b <= 1 and 0 <= c <= 1)

def createFullCoords(x, y, t):
    """
    returns for coordinates (x, y) in tria t (list of 3 vertex list with all coords of vertices of t) the full coordinates
    """
    v = [x, y]
    l0, l1 = PointLocationInTria(v, t) # returns length for vectors from point t3 with l0*(t2-t0) and l1*(t2-t1)
    elevation = t[2][2] + l0 * (t[0][2] - t[2][2])  + l1 * (t[1][2] - t[2][2])
    if elevation < -32765:
        v.append(-32768.0) #append correct elevation (stays -32768.0 in case of raster)
    else:
        v.append(elevation)
    v.extend([0, 0]) ### leave normal vectors to 0 ##### TO BE ADAPTED IN CASE OF NO RASTER !!!!!
    for i in range(5, len(t[0])): #go through optional coordinates s/t values based on first vertex in tria
        v.append(t[2][i] + l0 * (t[0][i] - t[2][i])  + l1 * (t[1][i] - t[2][i]))
    return v
        
def cutEdges(area, v, w, accuracy = 10): #cuts all edges of trias that intersect with segment from vertex v to w in area, if existing vertices/edges are closer than accuracy in m, they will be used  
    ## Note: If edge is completely within a tria then nothing is cut. In that case the segment has to be inserted, then by just inserting the vertices of the edge 
    #### TBD: function adds new vertices for each tria. Already added trias in neighbour vertex should be re-used !!!!!! ###########################
    cPs = [] #functions returns a list of vertices ) that are lying on intersection line v to w (either created or used within accuracy)
    new_trias = [] #trias to be added when cutting in patch
    old_trias = [] #old trias to be removed when cutting in patch
    for t in area:
        iv = [] #list of intersection vertices between line vw and edges of triangle t, could be between 0 and 3
        for edge in range(3): # go through edges of tria by index numbers
            cuttingPoint = intersection(v, w, t[edge][0:2], t[(edge+1)%3][0:2]) # modulo 3 returns to first vertex to close tria
            if cuttingPoint:
                existing_vertex_close = False
                for i in [edge, (edge+1)%3]: #check if v is too close to vertex of tria ONLY FOR ADJECENT vertices on edge  
                    if distance(t[i][0:2], cuttingPoint) < accuracy:
                        log.info("   Cutting Point {} not inserted as too close to vertex of the triangle it is in.".format(cuttingPoint))
                        cPs.append(t[i][0:2]) ### Attention: adds all coordinates of t that are within accuracy, not just one
                        existing_vertex_close = True
                if not existing_vertex_close:
                        cuttingPoint = (cuttingPoint[0], cuttingPoint[1], edge) #store number of edge cut as third element
                        iv.append(cuttingPoint)
                        cPs.append((cuttingPoint[0], cuttingPoint[1]))
        if len(iv) == 2:
            if iv[1][0] > iv[0][0] or (iv[1][0] == iv[0][0] and iv[1][1] > iv[0][1]):
                iv[0], iv[1] = iv[1], iv[0] #make sure to always start with most western (or most southern if both have same west coordinate) cutting Point in order to always get same mesh for overlays
            log.info("   Two intersections found at {}.".format(iv))
            
            edge = iv[0][2] ## start with edge with west/southern cutting point having index v_Index-1 (was inserted first)
            if iv[1][2] == (iv[0][2] + 1) % 3: #second cutting point lies on next edge
                new_v0 = createFullCoords(iv[0][0], iv[0][1], t)
                new_v1 = createFullCoords(iv[1][0], iv[1][1], t)
                new_trias.append([ new_v0, t[(edge+1)%3], new_v1, None, None, None, t[6]]) 
                new_trias.append([ new_v0, new_v1, t[(edge+2)%3], None, None, None, t[6]])
                new_trias.append([ new_v0, t[(edge+2)%3], t[edge], None, None, None, t[6]])
                #new_trias.append([ [iv[0][0], iv[0][1], 0, 0, 0], [t[(edge+1)%3][0], t[(edge+1)%3][1], 0, 0, 0], [iv[1][0], iv[1][1], 0, 0, 0], 0, 0, 0, 0])
                #new_trias.append([ [iv[0][0], iv[0][1], 0, 0, 0], [iv[1][0], iv[1][1], 0, 0, 0], [t[(edge+2)%3][0], t[(edge+2)%3][1], 0, 0, 0], 0, 0, 0, 0 ])
                #new_trias.append([ [iv[0][0], iv[0][1], 0, 0, 0], [t[(edge+2)%3][0], t[(edge+2)%3][1], 0, 0, 0], [t[edge][0], t[edge][1], 0, 0, 0], 0, 0, 0, 0 ])
            else: #second cutting point must lie on previous edge
                new_v0 = createFullCoords(iv[0][0], iv[0][1], t)
                new_v1 = createFullCoords(iv[1][0], iv[1][1], t)
                new_trias.append([ new_v0, t[(edge+1)%3], t[(edge+2)%3], None, None, None, t[6] ])
                new_trias.append([ new_v0, t[(edge+2)%3], new_v1, None, None, None, t[6] ])
                new_trias.append([ new_v0,  new_v1, t[edge], None, None, None, t[6] ])
                #new_trias.append([  [iv[0][0], iv[0][1], 0, 0, 0], [t[(edge+1)%3][0], t[(edge+1)%3][1], 0, 0, 0],  [t[(edge+2)%3][0], t[(edge+2)%3][1], 0, 0, 0], 0, 0, 0, 0 ])
                #new_trias.append([  [iv[0][0], iv[0][1], 0, 0, 0], [t[(edge+2)%3][0], t[(edge+2)%3][1], 0, 0, 0], [iv[1][0], iv[1][1], 0, 0, 0], 0, 0, 0, 0 ])
                #new_trias.append([  [iv[0][0], iv[0][1], 0, 0, 0], [iv[1][0], iv[1][1], 0, 0, 0], [t[edge][0], t[edge][1], 0, 0, 0], 0, 0, 0, 0 ])
            old_trias.append(t)                   
        elif len(iv) == 1:
            edge = iv[0][2]
            log.info("   One intersection found at {}.".format(iv))
            new_v0 = createFullCoords(iv[0][0], iv[0][1], t)
            new_trias.append([ t[edge], new_v0, t[(edge+2)%3], None, None, None, t[6] ])
            new_trias.append([ new_v0, t[(edge+1)%3], t[(edge+2)%3], None, None, None, t[6] ])
            #new_trias.append([ [t[edge][0], t[edge][1], 0, 0, 0], [iv[0][0], iv[0][1], 0, 0, 0],  [t[(edge+2)%3][0], t[(edge+2)%3][1], 0, 0, 0], 0, 0, 0, 0 ])
            #new_trias.append([ [iv[0][0], iv[0][1], 0, 0, 0], [t[(edge+1)%3][0], t[(edge+1)%3][1], 0, 0, 0], [t[(edge+2)%3][0], t[(edge+2)%3][1], 0, 0, 0], 0, 0, 0, 0 ])
            old_trias.append(t)
    for nt in new_trias:
        area.append(nt)
    for ot in old_trias:
        area.remove(ot)
    return cPs 




def showProgress(percentage):
    if percentage%10 == 0:
        log.info("dsf processing at {}%".format(percentage))


dsf_filename = "D:\\Programmierung\\mexp_0_1\\test.dsf"
log = defineLog('mexp', 'INFO', 'INFO') #no log on console for EXE version
log.info("Started mexp Version: {}".format(mexp_VERSION))

dsf = XPLNEDSF(LogName, showProgress)
log.info("Start reading dsf file: {}".format(dsf_filename))
dsf.read(dsf_filename)  # returns value > 0 in case of errors
log.info("dsf file with {} pools of vertices and {} patches read.".format(len(dsf.V), len(dsf.Patches)))

exarea = extractMeshArea(dsf, 50.4085, 50.4106, -125.1370, -125.1262)
print(exarea)
cps = cutEdges(exarea, [-125.147, 50.412, ], [-125.119, 50.407])
print(cps)
print("----")
log.info(exarea)
kmlExport(dsf, 50.4085, 50.4106, -125.1370, -125.1262, "D:\\Programmierung\\mexp_0_1\\test_without.dsf")
kmlExport2(dsf, exarea, "D:\\Programmierung\\mexp_0_1\\test_extract.dsf")
