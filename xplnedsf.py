#******************************************************************************
#
# xplnedsf.py        Version 0.4.3
# ---------------------------------------------------------
# Python module for reading and writing X_Plane DSF files.
#
#
# WARNIG: This code is still under development and may still have some errors.
#         In case you use it be very careful!
#
# Copyright (C) 2020 by schmax (Max Schmidt)
#
# This source is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option)
# any later version.
#
# This code is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.  
#
# A copy of the GNU General Public License is available at:
#   <http://www.gnu.org/licenses/>. 
#
#******************************************************************************

# NEW Version 0.4.2 Function PointLocationInTria returns also tuple if denom is zero
# NEW Version 0.4.3 Avoid errors when encoding pools with negative value or above 65535; values with -1 and 65536 come probably from rounding
# NEW Version 0.4.3 Avoid error when top triangle of runway does not exist, due to existing verctices in mesh not adding enough vertices on runway side A

from os import path, stat #required to retrieve length of dsf-file
from struct import pack, unpack #required for binary pack and unpack
from hashlib import md5 #required for md5 hash in dsf file footer
from logging import StreamHandler, getLogger, Formatter #for output to console and/or file
from io import BytesIO #required to go through bytes of a read 7ZIP-File
from math import sin, cos, sqrt, atan2, radians # for distance calculation etc.

try:
    import py7zlib
except ImportError:
    PY7ZLIBINSTALLED = False
else:
    PY7ZLIBINSTALLED = True


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

def PointInPoly(p, poly):  # test wether a point p with [lat, lon] coordinates lies in polygon (list of [lat, lon] pairs       
    # counts number of intersections from point outside poly to p on same y-coordinate, if it is odd the point lies in poly
    # to avoid intersection at vertex of poly on same y-coordinate, such points are shifte about 1mm above for testing intersection
    count = 0
    for i in range(len(poly) - 1):  # for all segments in poly
        epsilon0, epsilon1 = (0, 0)  # added to p's y coordinate in case p is on same y-coordinate than according vertex of segment
        if poly[i][1] < p[1] and poly[i + 1][1] < p[1]:  # if vertices of segment below y-coordinate of p, no intersection
            continue
        if poly[i][1] > p[1] and poly[i + 1][1] > p[1]:  # if vertices of segment above y-coordinate of p, no intersection
            continue
        if poly[i][1] == p[1]:
            epsilon0 = 0.00000001
        if poly[i + 1][1] == p[1]:
            epsilon1 = 0.00000001
        x = intersection([poly[i][0], poly[i][1] + epsilon0], [poly[i + 1][0], poly[i + 1][1] + epsilon1], [181, p[1]], p)
        if x:
            count += 1
    if count % 2:
        return True  # odd number of intersections, so p in poly
    else:
        return False  # even number of intersection, so p outside poly

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

         

class XPLNEpatch:
    def __init__(self, flag, near, far, poolIndex, defIndex):
        self.flag = flag
        self.near = near
        self.far = far
        self.defIndex = defIndex
        self.cmds = []
        self.trias = []  ###### NEW: store triangles calculated with triangles()  ###### TBD: really needed or use always functions triangles() and directly trias2cmds() to store changes?? ### NEW ## TBD ###
        self.minx = 181  #set out of bound values first  
        self.maxx = -181 #will be set to bounding rectangle values for patch
        self.miny = 91   #to quickly see if the patch applies for a certain ares
        self.maxy = -91  #values set correctly with function _setPatchBoundary()_ after being read  
    def triangles(self): #returns triangles as a list l of [3 x vertexes] that are defined by commands c of thte patch where each vertex of triangle is a pair of index to pool p and vertex        
        l = []
        p = None #current pool needs to be defined with first command
        for c in self.cmds:
            if c[0] == 1: # Pool index changed within patch, so change
                p = c[1]
            elif c[0] == 23: # PATCH TRIANGLE
                for v in range (1, len(c), 3): #skip value (command id)
                    l.append( [ [p, c[v]], [p, c[v + 1]], [p, c[v + 2]] ])
            elif c[0] == 24: # PATCH TRIANGLE CROSS POOL
                for v in range (1, len(c), 6): #skip value (command id)
                    l.append( [ [c[v], c[v+1]], [c[v+2], c[v+3]], [c[v+4], c[v+5]] ] )
            elif c[0] == 25: # PATCH TRIANGLE RANGE 
                for v in range(c[1], c[2] - 1, 3): #last index has one added 
                    l.append( [ [p, v], [p, v + 1], [p, v + 2] ] )
            elif c[0] == 26: # PATCH TRIANGLE STRIP
                for v in range (3, len(c)): #skip first three values - including command id are not needed any more
                    if v % 2: ##### NEW: Strip 1,2,3,4,5 refers to triangles 1,2,3 2,4,3 3,4,5 
                        l.append( [ [p, c[v-2]], [p, c[v-1]], [p, c[v]] ] )  
                    else:                                                    
                        l.append( [ [p, c[v-2]], [p, c[v]], [p, c[v-1]] ] )  
            elif c[0] == 27: # PATCH TRIANGLE STRIP CROSS POOL             
                for v in range (6, len(c), 2): #skip first values - including command id are not needed any more
                    if v % 4: ##### NEW: Strip 1,2,3,4,5 refers to triangles 1,2,3 2,4,3 3,4,5 
                        l.append( [ [c[v-5], c[v-4]], [c[v-3], c[v-2]], [c[v-1], c[v]] ] )  
                    else:                                                                   
                        l.append( [ [c[v-5], c[v-4]], [c[v-1], c[v]], [c[v-3], c[v-2]] ] )  
            elif c[0] == 28:  # PATCH TRIANGLE STRIP RANGE
                for v in range(c[1], c[2] - 2): # last index has one added so -2 as we go 2 up and last in range is not counted
                    if (v - c[1]) % 2:  ##### NEW: Strip 1,2,3,4,5 refers to triangles 1,2,3 2,4,3 3,4,5 
                        l.append( [ [p, v], [p, v + 2], [p, v + 1] ] ) 
                    else:                                              
                        l.append( [ [p, v], [p, v + 1], [p, v + 2] ] ) 
            elif c[0] == 29: # PATCH TRIANGLE FAN                          ##(all FAN commands not yet tested!!!!!!!!)##
                for v in range (3, len(c)): #skip first three values - including command id are not needed any more
                    l.append( [ [p, c[1]], [p, c[v-1]], [p, c[v]] ] ) #at c[1] is the center point id of the fan
            elif c[0] == 30: # PATCH TRIANGLE FAN CROSS POOL
                for v in range (6, len(c), 2): #skip first values - including command id are not needed any more
                    l.append( [ [c[1], c[2]], [c[v-3], c[v-2]], [c[v-1], c[v]] ] ) #at c[1] is the pool index for fan center and at c[2] the point id
            elif c[0] == 31:  # PATCH FAN RANGE
                for v in range(c[1], c[2] - 2): # at c[1] center point; last index has one added so -2 as we go 2 up and last in range is not counted
                    l.append( [ [p, c[1]], [p, v + 1], [p, v + 2] ] )                    
        return l
              
    def trias2cmds(self):  ## For the moment only uses single tirangle CMDS for different pools ## To be updeated 
        self.cmds = [self.cmds[0]]
        i = 0 #counts number of trias
        c = [] #builds single commands
        for t in self.trias:
            if (i % 85) == 0: #start new triangle command at beginning and when max length is reached (3 value pairs per triangle, so 255 / 3 = 85)
                if len(c) > 1: #are there already values for commands
                    self.cmds.append(c) #append them the the commands
                c = [24]
            c.extend(t[0])
            c.extend(t[1])
            c.extend(t[2])
            i += 1
        if len(c) > 1:
            self.cmds.append(c) #add the final command


class XPLNEraster: #Stores data of Raster Atoms (each dsf could have serverl raster layers)
    def __init__(self):
        self.ver = None #Version of Raster
        self.bpp = None #bytes per pixel in raster data
        self.flags = None #flags about Raster centricity and what data type used
        self.width = None #of area in pixel
        self.height = None #of area in pixel
        self.scale = None #scale factor for height values
        self.offset = None #offset for heigt values
        self.data = [] #will store final raster heigt values (after scaling and adding offset) in 2-dimensional list: [pixel x] [pixel y]
        

class XPLNEDSF:   
    def __init__(self, logname='__XPLNEDSF__', statusfunction = "stdout"):
        if logname != "_keep_logger_": self._log_ = self._setLogger_(logname)
        if statusfunction != "_keep_statusfunction_": self._statusfunction_ = statusfunction
        self._DEBUG_ = True if self._log_.getEffectiveLevel() < 20 else False #have DEBUG value in order to call only logger for debug if DEBUG enabled  --> saves time
        self._progress_ = [0, 0, 0] #progress as 3 list items (amount of bytes read/written in percant and shown, read/written but not yet shown as number of bytes) and number of bytes to be processed in total
        self._Atoms_ = {} #dictonary containg for every atom in file the according strings
        self._AtomStructure_ = {'DAEH' : ['PORP'], 'NFED' : ['TRET', 'TJBO', 'YLOP', 'WTEN', 'NMED'], 'DOEG' : ['LOOP', 'LACS', '23OP', '23CS'], 'SMED' : ['IMED', 'DMED'], 'SDMC' : []}
        self._AtomList_ = ['DAEH', 'NFED', 'DOEG', 'SMED', 'SDMC', 'PORP', 'TRET', 'TJBO', 'YLOP', 'WTEN', 'NMED', 'LOOP', 'LACS', '23OP', '23CS', 'IMED', 'DMED']
        self._AtomOfAtoms_ = ['DAEH', 'NFED', 'DOEG', 'SMED']
        self._MultiAtoms_ = ['LOOP', 'LACS', '23OP', '23CS', 'IMED', 'DMED'] #These Atoms can occur severl times and therefore are stored as list in self._Atoms_
        self._CMDStructure_ = {1 : ['H'], 2 : ['L'], 3 : ['B'], 4 : ['H'], 5 : ['L'], 6 : ['B'], 7 : ['H'], 8 : ['HH'], 9 : ['', 'B', 'H'], 10 : ['HH'], 11 : ['', 'B', 'L'], 12 : ['H', 'B', 'H'], 13 : ['HHH'], 15 : ['H', 'B', 'H'], 16 : [''], 17 : ['B'], 18 : ['Bff'], 23 : ['', 'B', 'H'], 24 : ['', 'B', 'HH'], 25 : ['HH'], 26 : ['', 'B', 'H'], 27 : ['', 'B', 'HH'], 28 : ['HH'], 29 : ['', 'B', 'H'], 30 : ['', 'B', 'HH'], 31 : ['HH'], 32 : ['', 'B', 'c'], 33 : ['', 'H', 'c'], 34 : ['', 'L', 'c']}
        self._CMDStructLen_ = {1 : [2], 2 : [4], 3 : [1], 4 : [2], 5 : [4], 6 : [1], 7 : [2], 8 : [4], 9 : [0, 1, 2], 10 : [4], 11 : [0, 1, 4], 12 : [2, 1, 2], 13 : [6], 15 : [2, 1, 2], 16 : [0], 17 : [1], 18 : [9], 23 : [0, 1, 2], 24 : [0, 1, 4], 25 : [4], 26 : [0, 1, 2], 27 : [0, 1, 4], 28 : [4], 29 : [0, 1, 2], 30 :  [0, 1, 4], 31 : [4], 32 : [0, 1, 1], 33 : [0, 2, 1], 34 : [0, 4, 1]}
        self.CMDS = [] #unpacked commands
        self.Patches = [] #mesh patches, list of objects of class XPLNEpatch
        self.V = [] # 3 dimensional list of all vertices whith V[PoolID][Vertex][xyz etc coordinates]
        self.V32 = [] # same as V but for 32-bit coordinates
        self.Scalings = [] # 3 dimensional list of all scale multipliers and offsets for all pools and planes with Scalings[PoolID][coordinate Plane][m / o]
        self.Scal32 = [] # same as Scalings but for vertices with 32-bit coordinates
        self.Raster = [] #raster layers of file
        self.Polygons = [] #for each Polygon Definition (same order as in DefPolygon) a list of PoolId, CommandData values of Polygon
        self.Objects = [] #for each Object Definition (same order as in DefObjects) a list of PoolId, CommandData values of Polygon
        self.Networks = [] #list of network junks where in first positions road subtype, Junction offset, PoolIndex are stored, following by commands for these settings
        self.Properties = {} #dictionary containing the properties of the dsf file stored in HEAD -> PROP
        self.DefTerrains = {} #dictionary containing for each index number (0 to n-1) the name of Terrain definition file
        self.DefObjects = {}  #dictionary containing for each index number (0 to n-1) the name of Object definition file
        self.DefPolygons = {}  #dictionary containing for each index number (0 to n-1) the name of Polygon definition file
        self.DefNetworks = {}  #dictionary containing for each index number (0 to n-1) the name of Network definition file; actually just one file per dsf-file
        self.DefRasters = {}   #dictionary containing for each index number (0 to n-1) the name of Raster definition (for the moment assuing "elevaiton" is first with index 0)
        self.submP = set() #set with indices of vertex pools that allow submeter elevation; not given by dsf file but checked when submeter vertices should be added   ##### NEW ######
        self._log_.info("Class XPLNEDSF initialized.")

    def _setLogger_(self, logname):
        if logname == '__XPLNEDSF__': #define default logger if nothing is set
            logger = getLogger('XPLNEDSF')
            logger.setLevel('INFO')
            logger.handlers = []  # Spyder/IPython currently does not remove existing loggers, this does the job
            stream_handler = StreamHandler()
            stream_handler.setLevel('INFO')
            formatter = Formatter('%(levelname)s: %(message)s')
            stream_handler.setFormatter(formatter)
            logger.addHandler(stream_handler)
        else:
            logger = getLogger(logname + '.' + __name__) #use name of existing Logger in calling module to get logs in according format or change basicConfig for logs
        return logger

    def _updateProgress_(self, bytes): #updates progress with number of read/written bytes and calls show-function
        self._progress_[1] += bytes
        if self._progress_[2] <= 0:
            return 1 #size to be reached not set; no sense to track progress
        currentprogress = round(100 * self._progress_[1] / self._progress_[2])
        if currentprogress:
            self._progress_[0] += currentprogress
            if self._progress_[0] > 100:
                self._progress_[0] = 100 #stopp at 100%
                self._progress_[1] = 0
            else:
                self._progress_[1] -= int (currentprogress * self._progress_[2] / 100)
            if self._statusfunction_ == "stdout":
                print ('[{}%]'.format(self._progress_[0]), end='', flush = True)
            elif self._statusfunction_ != None:
                self._statusfunction_(self._progress_[0])
    
    
    def _TopAtomLength_(self, id): #calculates the lengthe of an AtomOfAtoms with name id including all sub-atoms
        l = 0
        for sub_id in self._AtomStructure_[id]:
            if sub_id in self._MultiAtoms_:
                if sub_id in self._Atoms_: #check that sub-atom is really present 
                    for a in self._Atoms_[sub_id]: #for multiple atoms we have to add length of each instance
                        l += len(a) + 8  # add 8 bytes for each AtomID + Length header
            else: #for single atoms just lengths of atom
                l += len(self._Atoms_[sub_id]) + 8 # add 8 bytes for AtomID + Length header
        return l + 8  # add 8 for header of TopAtom itself

    
    def _GetStrings_(self, atom): #Returning a list of strings packed in a string table atom
        l = []
        i = 0
        while i < len(atom):
            j = i
            while atom[j] != 0 and j < len(atom):
                j += 1
            l.append(atom[i:j].decode("utf-8"))
            i = j + 1
        return l
    
    def _PutStrings_(self, d): #encodes strings in dictionary d to binary to be stored in updated atom   ########## NEW 5 #############
        b = b'' #binary encoding to be returned
        for i in range(len(d)):
            b += d[i].encode("utf-8")
            b += b'\x00'
        return b
            
            
    def _extractProps_(self): #extracts the properties of the dsf file stored in atom PROP within atom HEAD; stores them to dictionary Properties
        l = self._GetStrings_(self._Atoms_['PORP'])
        for i in range(0, len(l), 2):
            self.Properties[l[i]] = l[i+1] #the list contains property(dictionary key) and values one after each other

     
    def _extractDefs_(self): #extracts the definition (DEFN) atoms TERT, OBJT, POLY, NETW, DEMN and stores them in dictionarys
        l = self._GetStrings_(self._Atoms_['TRET'])
        self.DefTerrains = dict(zip(range(len(l)), l))
        l = self._GetStrings_(self._Atoms_['TJBO'])
        self.DefObjects = dict(zip(range(len(l)), l))
        l = self._GetStrings_(self._Atoms_['YLOP'])
        self.DefPolygons = dict(zip(range(len(l)), l))
        l = self._GetStrings_(self._Atoms_['WTEN'])
        self.DefNetworks = dict(zip(range(len(l)), l))        
        l = self._GetStrings_(self._Atoms_['NMED'])
        self.DefRasters = dict(zip(range(len(l)), l))
        self._updateProgress_(self._TopAtomLength_('NFED'))

        
    def _encodeDefs_(self): #encodes the definition dictionaries in dsf object to (DEFN) atoms     
        self._Atoms_['TRET'] = self._PutStrings_(self.DefTerrains)
        self._Atoms_['TJBO'] = self._PutStrings_(self.DefObjects)
        self._Atoms_['YLOP'] = self._PutStrings_(self.DefPolygons)
        self._Atoms_['WTEN'] = self._PutStrings_(self.DefNetworks)
        self._Atoms_['NMED'] = self._PutStrings_(self.DefRasters)
        ###### TBD: other definitions #########################

        
    def _extractPools_(self, bit = 16):  
        if bit == 32:
            self._log_.info("Start to unpack and extract {} pools ({} bit)...".format(len(self._Atoms_['23OP']), bit))
            atomstring = self._Atoms_['23OP']
            V = self.V32
            ctype = "<L"
            size = 4 #bytes read per coordinate in 32 bit pool
        else: #assuming the standard 16bit Pool case
            self._log_.info("Start to unpack and extract {} pools ({} bit)...".format(len(self._Atoms_['LOOP']), bit))
            atomstring = self._Atoms_['LOOP']
            V = self.V
            ctype = "<H"
            size = 2 #bytes read per coordinate in 16 bit pool
        for s in atomstring: #goes through all Pools read; string s has to be unpacked
            nArrays, nPlanes = unpack('<IB', s[0:5])
            if self._DEBUG_: self._log_.debug("Pool number {} has {} Arrays (vertices) with {} Planes (coordinates per vertex)!".format(len(V), nArrays, nPlanes))
            V.append([]) #the current pool starts empty
            for i in range(nArrays): ## span up multi-dimensional array for the new pool of required size (number of vertices in pool)
                V[-1].append([])
            pos = 5 #position in string s
            for n in range(nPlanes):
                encType, = unpack('<B', s[pos : pos + 1])
                pos += 1
                if self._DEBUG_: self._log_.debug("Plane {} is encoded: {}".format(n, encType))
                if encType < 0 or encType > 3: #encoding not defined
                    self._log_.error("Stopp reading pool because not known encoding of plane found!!!")
                    return [] ##This means we return empty pool, which can be used to detect error        
                i = 0  #counts how many arrays = vertices have been read in plane n
                while i < nArrays:
                    if encType >= 2: #this means data is run-length encoded
                        runLength, = unpack('<B', s[pos : pos + 1])
                        pos += 1
                    else: #no run-length encoding (not tested yet!!!!!)
                        runLength = 1 #just read single values until end of this plane  
                    if runLength > 127: #means the following value is repeated
                        v, = unpack(ctype, s[pos : pos + size]) #repeated value
                        pos += size
                        runLength -= 128 #only value without 8th bit gives now the number of repetitions
                        while runLength > 0:
                            V[-1][i].append(v)
                            runLength -= 1
                            i += 1
                    else:
                        while runLength > 0: #now just reading individual values
                            v, = unpack(ctype, s[pos : pos + size])
                            pos += size
                            V[-1][i].append(v)
                            runLength -= 1
                            i += 1
                if encType == 1 or encType == 3: #values are also stored differenced
                    for i in range (1, nArrays):  
                        V[-1][i][n] = (V[-1][i][n] + V[-1][i-1][n]) % 65536  #undo differntiation and modulo for wrapping two byte unsigned integer
                                                         ##### NEW TBD #### FOR 32 Bit Pools adapt Moduls ??!! #################### NEW TBD ######################## NEW TBD ########
            self._updateProgress_(len(s))


    def _encodeRunLength_(self, l):  # yields runlength encoded value pairs of list l
        count = 1 #counting repetitions
        prev = l[0] #the previous value in list starts with first value in the list
        l.append('_END_') #as end of list symbol
        individuals = [] #list of individual values (non repeating ones)
        for value in l[1:]:
            if value != prev or count == 127: #non repeating value or maximum of repeating values reached
                if len(individuals) == 127:  #if maximum length of indivdiual values is reached
                    yield(len(individuals),individuals)
                    individuals = []
                if count > 1:  #we have now repeating values
                    if len(individuals): #individual values before have to be returned if existent
                        yield(len(individuals), individuals) 
                    yield(count + 128, [prev]) #now return the repeating values with highest bit (128) set of byte
                    individuals = []
                else: #non repeating values, so add previous to the individuals
                    individuals.append(prev)
                count = 0
            prev = value
            count += 1
        if len(individuals): #return still existing values in list
            yield(len(individuals), individuals)
        
    
    def _encodePools_(self): #overwrites current Pool atom with actual values of all vertices
        ###### FOR THE MOMENT ONLY WORKING FOR 16 BIT POOLS !!!!!!!!! ################################
        self._log_.info("Start to encode all pools... (read values of vertices are overwritten)")
        self._Atoms_['LOOP'] = [] #start new (future version also think of creating Pool atom in case it new dsf file will be created !!!!!!!!!!)
        ####### This version only stores pools in differentiated run-length encoding !!! #############
        ## Start with differentiation of all values ##
        for p in self.V: # go through all pools
            encpool = pack("<IB",len(p),len(p[0])) ## start string of binary encoded pool number of arrays and number of planes (taken from first vertex)
            for n in range(len(p[0])): # go through all planes; number of planes just taken from first vertex in pool
                plane = [p[0][n]] #set plane to first value (for starts with second value to calculete differences)
                for i in range(1, len(p)): #go through all values of a plane for differntiation and append to current plane
                    plane.append((p[i][n] - p[i-1][n]) % 65536)  #Calculate difference to previous value AND take care of wrapping for two byte unsigned integer
                encpool += pack('<B',3) #plane will be encoded differntiated + runlength
                if p[0][n] < 0:  #### NEW 4.3 #####
                    self._log_.warning("In pool {} negative value {} to be encoded. Set to 0.".format(self.V.index(p), p[0][n]))   
                    p[0][n] = 0
                if p[0][n] > 65535: #### NEW 4.3 #####
                    self._log_.warning("In pool {} exceeding value {}  to be encoded. Set to 65535.".format(self.V.index(p), p[0][n]))
                    p[0][n] = 65535      
                pack('<H',p[0][n])
                ## Now perform run-length encoding ##
                for rlpair in self._encodeRunLength_(plane):
                    encpool += pack('<B', rlpair[0]) #encode runlength value
                    for v in rlpair[1]:
                        if v < 0:   #### NEW 4.3 #####
                            self._log_.warning("In pool {} negative value {} to be encoded. Set to 0.".format(self.V.index(p), v))  
                            v = 0
                        if v > 65535: #### NEW 4.3 #####
                            self._log_.warning("In pool {} exceeding value {}  to be encoded. Set to 65535.".format(self.V.index(p), v))
                            v = 65535                            
                        encpool += pack('<H', v)
            self._updateProgress_(len(encpool)) 
            self._Atoms_['LOOP'].append(encpool)
        self._log_.info("Encoding of pools finished.")
 

    def _extractScalings_(self, bit = 16): #extract scaling values from atoms and writes them to according lists
        self._log_.info("Start to unpack and extract all saclings of {} bit pools.".format(bit))
        if bit == 32:
            atomstring = self._Atoms_['23CS']
            Scalings = self.Scal32
        else: #assuming the standard 16bit Pool case
            atomstring = self._Atoms_['LACS']
            Scalings = self.Scalings
        for s in atomstring: #goes through all Scalings read; string s has to be unpacked
            Scalings.append([]) #new scaling entry for next pool   
            for i in range(int(len(s) / 8) ):  #for each scaling tuple - 4 bytes multiplier and 4 bytes offset
                m, o = unpack('<ff', s[ i*8 : i*8 + 8] )  #retrieve multiplier and offset
                Scalings[-1].append([m, o])

    def _packAllScalings_(self): #packs the all Scalings incl. Scal32 values to binary atom strings to be later written to file 
        self._log_.info("Start to pack and all saclings for pools.")
        self._Atoms_['LACS'] = []
        self._Atoms_['23CS'] = []
        for s in self.Scalings:
            encscal = b''
            for i in s:
                encscal += pack('<ff', i[0], i[1])
            self._Atoms_['LACS'].append(encscal)
        for s in self.Scal32:
            encscal = b''
            for i in s:
                encscal += pack('<ff', i[0], i[1])
            self._Atoms_['23CS'].append(encscal)                    

                
    def _scaleV_(self, bit = 16, reverse = False): #applies scaling to the vertices stored in V and V32
        if reverse:
            self._log_.info("Start to de-scale all {} bit pools.".format(bit))
        else:
            self._log_.info("Start to scale all {} bit pools.".format(bit))
        if bit == 32:
            V = self.V32
            Scalings = self.Scal32
        else: #assuming the standard 16bit Pool case
            V = self.V
            Scalings = self.Scalings
            
        if len(V) != len(Scalings):
            self._log_.error("Amount of Scale atoms does not equal amount of Pools!!")
            return 1
        for p in range(len(V)): #for all Pools
            if V[p] == []: ###There can exist empty pools that have to be skipped for scaling!!!
                self._log_.info("Empty pool number {} not scaled!".format(p))
                break
            if len(V[p][0]) != len(Scalings[p]): #take first vertex as example to determine number of coordinate planes in current pool
                self._log_.error("Amount of scale values for pool {} does not equal the number of coordinate planes!!!".format(p))
                return 2
            for n in range(len(Scalings[p])): #for all scale tuples for that pool = all coordinate planes in pool
                if self._DEBUG_: self._log_.debug("Will now scale pool {} plane {} with multiplier: {} and offset: {}".format(p, n ,Scalings[p][n][0], Scalings[p][n][1]))                
                if float(Scalings[p][n][0]) == 0.0:
                    if self._DEBUG_: self._log_.debug("   Plane will not be scaled because scale is 0!")
                    break
                for v in range(len(V[p])): #for all vertices in current plane
                    if reverse: #de-scale vertices
                        V[p][v][n] = round((V[p][v][n] - Scalings[p][n][1]) * 65535 / Scalings[p][n][0])   #de-scale vertex v in pool p for plane n by subtracting offset and dividing by multiplyer
                    else:  #scale vertices
                        V[p][v][n] = (V[p][v][n] * Scalings[p][n][0] / 65535) + Scalings[p][n][1]  #scale vertex v in pool p for plane n with multiplyer and offset
               
                
    def _extractRaster_(self):  #extracts alll rasters from atoms and stores them in list
        self._log_.info("Extracting {} raster layers...".format(len(self._Atoms_['IMED'])))
        if len(self._Atoms_['IMED']) != len(self._Atoms_['DMED']):
            self._log_.error("Number of raster info atoms not equal to number of raster data atoms!!!")
            return 1           
        for rn in range(len(self._Atoms_['IMED'])):
            R = XPLNEraster()
            R.ver, R.bpp, R.flags, R.width, R.height, R.scale, R.offset = unpack('<BBHLLff', self._Atoms_['IMED'][rn])
            if self._DEBUG_: self._log_.debug("Info of new raster layer: {} {} {} {} {} {} {}".format(R.ver, R.bpp, R.flags, R.width, R.height, R.scale, R.offset))
            if R.flags & 1:  #signed integers to be read
                if R.bpp == 1:
                    ctype = "<b"
                elif R.bpp == 2:  ##### this is the only case tested so far !!!!!!!
                    ctype = "<h"
                elif R.bpp == 4:
                    ctype = "<i"
                else:
                    self._log_.error("Not allowed bytes per pixel in Raster Definition!!!")
                    return 2
            elif R.flags & 2: #unsigned integers to be read
                if R.bpp == 1:
                    ctype = "<B"
                elif R.bpp == 2:
                    ctype = "<H"
                elif R.bpp == 4:
                    ctype = "<I"
                else:
                    self._log_.error("Not allowed bytes per pixel in Raster Definition!!!")
                    return 3
            elif not (R.flags & 1) and not (R.flags & 2): #neither first nor second bit set means that 4 byte float has to be read
                if R.bpp == 4:
                    ctype ="<f"
                else:
                    self._log_.error("Not allowed bytes per pixel in Raster Definition!!!")
                    return 4

            for x in range(0, R.bpp * R.width, R.bpp): #going x-wise from east to west just the bytes per pixes
                line = []
                for y in range(0, R.bpp * R.height * R.width, R.bpp * R.width): #going y-wise from south to north, always jumping over the width of each x-line
                    v, = unpack(ctype, self._Atoms_['DMED'][rn][y + x : y + x + R.bpp]) #unpack R.bpp bytes at position x+y of raster data atom number rn
                    v = v * R.scale + R.offset # APPLYING SCALE + OFFSET 
                    line.append(v) #the pixel appended is to a line from south to north (y-line)
                R.data.append(line) #south to north lines are appended to each other
                self._updateProgress_(R.bpp * R.width) #update progress with number of bytes per raster line
            self.Raster.append(R) #so raster list of list is returned to be indexed by [x][y]
        self._log_.info("Finished extracting Rasters.")
   
    def _unpackCMDS_(self):
        i = 0 #position in CMDS String
        self._log_.info("Start unpacking of Commands.")
        current100kBjunk = 1 #counts processed bytes in 100kB junks
        while i < len(self._Atoms_['SDMC']):
            id, = unpack('<B', self._Atoms_['SDMC'][i : i+1])
            self.CMDS.append([id])
            i += 1
            if id in self._CMDStructure_:
                l = self._CMDStructLen_[id][0] #length of bytes to read
                if l > 0:
                    y = unpack('<' + self._CMDStructure_[id][0], self._Atoms_['SDMC'][i : i + l])
                    self.CMDS[-1].extend(y)
                    i += l
                if len(self._CMDStructLen_[id]) == 3: #read command with variable length
                    l = self._CMDStructLen_[id][1] #length of repeating value n to read
                    n, = unpack('<' + self._CMDStructure_[id][1], self._Atoms_['SDMC'][i : i + l]) #number n of repetitions
                    if id == 15:
                        n += 1 #id = 15 seems a special case that there is one index more than windings  ########??????????
                    i += l
                    l = self._CMDStructLen_[id][2] #length of repeated bytes
                    for j in range(n): 
                       y = unpack('<' + self._CMDStructure_[id][2], self._Atoms_['SDMC'][i : i + l])
                       self.CMDS[-1].extend(y)
                       i += l
            else:
                if id == 14: #special double packed case, which is explicetly treated separate and has special format returning lists inside commands !!!
                             ###### not tested yet !!!!!!!! ########
                    parameter, windings = unpack('<HB', self._Atoms_['SDMC'][i : i + 3])
                    i += 3
                    self.CMDS[-1].extend(parameter)
                    for w in range(windings):
                        indices, = unpack('<B', self._Atoms_['SDMC'][i : i + 1])
                        i += 1
                        windinglist = []
                        for m in range(indices):
                            y, = unpack('<H', self._Atoms_['SDMC'][i : i + 2])
                            i += 2
                            windinglist.extend(y)
                        self.CMDS[-1].append(windinglist)
                else: #command id not tretated here until now
                    self._log_.warning("Unknown command ID {} ignored!".format(id))
                    self.CMDS[-1] #delete already written id
            if self._DEBUG_: self._log_.debug("CMD id {}: {} (string pos next cmd: {})".format(self.CMDS[-1][0], self.CMDS[-1][1:], i))
            if i > current100kBjunk * 100000:
                self._updateProgress_(50000) #count only half of the length, other half by extractCMDS
                current100kBjunk += 1
        self._log_.info("{} commands haven been unpacked.".format(len(self.CMDS)))


    def _extractCMDS_(self): # extract CMDS and stores it as Mesh-Patches, Polygons, ...
        self._log_.info("Start to extract CMDS")
        for i in range(len(self.DefPolygons)):
            self.Polygons.append([]) #span list of empty lists for all defined poygon types
        for i in range(len(self.DefObjects)): 
            self.Objects.append([]) #span list of empty lists for all defined poygon types

        patchPoolIndex = None #poolIndex currently used in current patch; if different from current in CMDS then change command is written to cmds of patch

        flag_physical = None #1 if physical, 2 if overlay
        nearLOD = None  
        farLOD = None
        poolIndex = 0
        defIndex = 0
        subroadtype = 0
        junctionoffset = 0
        counter = 0
        amount_of_two_percent_CMDS = int(len(self.CMDS)/50)  ### NEW ###
        if amount_of_two_percent_CMDS == 0:                  ### NEW ###
            amount_of_two_percent_CMDS = 1 #have 1 as minimal value in order to avoid devision by zero below   ### NEW ###
        for c in self.CMDS:
            if c[0] == 1: # new pool selection
                poolIndex = c[1]
            elif c[0] == 2: # new junction offset
                junctionoffset = c[1]
            elif 3 <= c[0] <= 5: # new definition index
                defIndex = c[1]
            elif c[0] == 6: # new subtype for road
                subroadtype = c[1]
            elif 7 <= c[0] <= 8: #Object Command
                self.Objects[defIndex].append([poolIndex]) #new Polygond added for defIndex type and it starts with poolIndex from which its vertices are
                self.Objects[defIndex][-1].extend(c) #followed by the complete command to build it
            elif 9 <= c[0] <= 11: #Network Commands ### NEW: each Network command put in sublists, addtional [] inclueded !!!! 
                if self.Networks == []: #first network command, so start with first entry
                    self.Networks.append([[subroadtype, junctionoffset, poolIndex]])
                elif self.Networks[-1][0][0] != subroadtype or self.Networks[-1][0][1] != junctionoffset or self.Networks[-1][0][2] != poolIndex: #chang of relevant base settings
                    self.Networks.append([[subroadtype, junctionoffset, poolIndex]]) #sp new entry with new base-settings
                self.Networks[-1].extend([c]) #append complete command to build this network part on current base settings
            elif 12 <= c[0] <= 15: #Polygon Commands
                self.Polygons[defIndex].append([poolIndex]) #new Polygond added for defIndex type and it starts with poolIndex from which its vertices are
                self.Polygons[defIndex][-1].extend(c) #followed by the complete command to build it
            elif 16 <= c[0] <= 18:  # Add new Terrain Patch
                patchPoolIndex = None #poolIndex for new patch needs to be set as first command
                if c[0] == 17: # New Patch with new physical flag
                    flag_physical = c[1]
                elif c[0] == 18: # New Patch with new flag and LOD
                    flag_physical = c[1]
                    nearLOD = c[2]
                    farLOD = c[3]    
                p = XPLNEpatch(flag_physical, nearLOD, farLOD, poolIndex, defIndex)
                self.Patches.append(p)
            elif c[0] >= 23 and c[0] <= 31: # the command is about a patch, so store it to current patch
                if patchPoolIndex != poolIndex: #change of pool index within patch is possible
                    self.Patches[-1].cmds.append( [1, poolIndex] ) #change pool with patch via according command
                    patchPoolIndex = poolIndex #now within patch also current pool will be used
                if self.Patches[-1].defIndex != defIndex:
                    self._log_.error("Definition Index changed within patch. Aborted command extraction!")
                    return 1
                self.Patches[-1].cmds.append(c) 
            counter += 1
            if not counter % amount_of_two_percent_CMDS: #after every 2% processed of commands update progress  --> RAISES ERROR when less than 50 CMDS
                self._updateProgress_(round(len(self._Atoms_['SDMC']) / 100)) #count only half of the length, other half by unpackCMDS
        for p in self.Patches: #add to each patch its trias as list and the rectangle boundary coordinates of patch   
            p.trias = p.triangles()   ##NEW##
            self._setPatchBoundary_(p) ##NEW##
        self._log_.info("{} patches extracted from commands.".format(len(self.Patches)))
        self._log_.info("{} different Polygon types including there definitions extracted from commands.".format(len(self.Polygons)))
        self._log_.info("{} different Objects with placements coordinates extracted from commands.".format(len(self.Objects)))
        self._log_.info("{} different Network subtypes extracted from commands (could include double count).".format(len(self.Networks)))

    def _packCMDS_(self): #packs all CMDS of Object, Polygons, Networks and Patches in binary string to be later written to file #####NEW FUCTION### NEW ####
        ################################# TBD: Check function for polygons and object placements !!!!!! #############################################################
        self._log_.info("Start to pack CMDS")   #### NEW 4 : put definition of variables befor encCMD function, as values like poolIndex will be changed inside ##############
        #SET state variables
        flag_physical = None
        nearLOD = None  
        farLOD = None
        poolIndex = None
        defIndex = None
        subroadtype = None
        junctionoffset = None  #### set to 0 if directly set below as in X-Plane standard dsf files
        enccmds = bytearray() #these will be the encoded CMDS to be returned  
        ### local function for single CMD encoding ###
        def encCMD(c): #encodes single CMD in array c and returns its binary
            ecmd = b'' 
            id = c[0] #id of CMD stored as first value in list
            if id == 3 and c[1] > 255: #avoid errors if definition is too big for command id 3
                id = 4
            if id == 4 and c[1] > 65535: #avoid errors if definition is too big for command id 4
                id = 5
            ecmd += pack('<B', id) #encode CMD id
            i = 1 #index in value list to be packed
            if id in self._CMDStructure_: 
                for valtype in self._CMDStructure_[id][0]: #pack fixed length values of CMD
                    ecmd += pack('<' + valtype, c[i])
                    i += 1
                if len(self._CMDStructLen_[id]) == 3: #pack command with variable length 
                    vlength = int( (len(c) - i) / len(self._CMDStructure_[id][2]) ) #fixed values inclding CMD id are counting not for the variable length, so minus i #### TBD: Test that this value does not exceed 255!!!!!! ####
                    if vlength > 255: #make sure that length is not longer than could be encoded in one byte
                        self._log_.error("Length of CMD with id {} and length {} does exceed 255 and will not be included!!!".format(id, vlength))
                        return b''
                    if id == 15:
                        vlength -= 1 #id = 15 seems a special case that there is one index more than windings 
                    ecmd += pack('<B', vlength) #count packed
                    while i < len(c): #now pack the variable length value list
                        for valtype in self._CMDStructure_[id][2]: #pack fixed length individual value
                            ecmd += pack('<' + valtype, c[i])
                            i += 1
            else:
                ######### TBD: include special code for CMD 14 here!!! ########## TBD ########### TBD ########### TBD ########### TBD ############### TBD ###########
                self._log_.error("CMD id {} is not supported (CMD 14 not implemented yet)! Will be skipped!!".format(id))
                return b''
            return ecmd
        ### end of inner function to pack single CMDS ###   
        for d in self.DefObjects: #for each object definition write according CMDS
            enccmds.extend(encCMD([3, d])) #definition set according to current definition id; function will handle if id > 255
            for c in self.Objects[d]:
                if c[0] != poolIndex: #Pool-Index is written before CMD; adapt index if it changes
                    enccmds.extend(encCMD([1, c[0]]))
                    poolIndex = c[0]
                enccmds.extend(encCMD(c[1:])) #now according command to place objects is encoded  #### TO BE TESTED if it works as part of array???????
        for d in self.DefPolygons: #for each polygon definition write according CMDS
            enccmds.extend(encCMD([3, d])) #definition set according to current definition id; function will handle if id > 255
            for c in self.Polygons[d]:
                if c[0] != poolIndex: #Pool-Index is written before CMD; adapt index if it changes
                    enccmds.extend(encCMD( [1, c[0]] ))
                    poolIndex = c[0]
                enccmds.extend(encCMD(c[1:])) #now according command to place objects is encoded  #### TO BE TESTED if it works as part of array???????        
        for d in self.Networks:
            if d[0][1] != junctionoffset: ## Order in org X-Plane files is with CMD id 2 at first
                enccmds.extend(encCMD([2, d[0][1]]))
                junctionoffset = d[0][1]
            if defIndex != 0:
                enccmds.extend(encCMD([3, 0])) #Currently there is only one Road-Defintion --> DefIndex set to 0
                defIndex = 0
            if d[0][2] != poolIndex:
                enccmds.extend(encCMD([1, d[0][2]]))
                poolIndex = d[0][2]
            if d[0][0] != subroadtype:  ## Order in org X-Plane files is with 6 at last
                enccmds.extend(encCMD([6, d[0][0]]))
                subroadtype = d[0][0]
            for c in d[1:]:
                enccmds.extend(encCMD(c))
        for d in self.Patches:
            if defIndex != d.defIndex:
                enccmds.extend(encCMD([3, d.defIndex])) #definition set according to current definition id; function will handle if id > 255
                defIndex = d.defIndex
            if poolIndex != d.cmds[0][1]: #Pool-Index is defined by first command and required to be defined directly before new Patch is defined!    
                enccmds.extend(encCMD([1, d.cmds[0][1]])) #include command for changing poolIndex
                poolIndex = d.cmds[0][1] #update state variable
            if nearLOD == d.near and farLOD == d.far:
                if flag_physical == d.flag:
                    enccmds.extend(encCMD([16]))
                else:
                    enccmds.extend(encCMD([17, d.flag]))
                    flag_physical = d.flag
            else:
                enccmds.extend(encCMD([18, d.flag, d.near, d.far]))
                farLOD = d.far
                nearLOD = d.near
                flag_physical = d.flag
            for c in d.cmds[1:]:   ##skip first command as this is pool defintion written above
                enccmds.extend(encCMD(c))
                if c[0] == 1: #change of poolIndex can happen within patch commands    
                    poolIndex = c[1] #therefore also state variable has to be adapted        
        self._Atoms_['SDMC'] = enccmds #Commands Atom now set to the now packed CMDS
        self._log_.info("Ended to pack CMDS")
 ########### END of NEW function _packCMDS_() #################               
                

    def _unpackAtoms_(self): #starts all functions to unpack and extract data froms strings in Atoms
        self._log_.info("Extracting properties and definitions.")
        if 'PORP' in self._Atoms_:
            self._extractProps_()
        else:
            self._log_.error("This dsf file has no properties defined!")
        if 'NFED' in self._Atoms_:    
            self._extractDefs_()
        else:
            self._log_.warning("This dsf file has no definitions.") 
        if 'IMED' in self._Atoms_:
            self._extractRaster_()
        else:
            self._log_.info("This dsf file has no raster layers.")
        if 'LOOP' in self._Atoms_:
            self._extractPools_(16)
            self._extractScalings_(16)
            self._scaleV_(16, False) #False that scaling is not reversed
            self._updateProgress_(len(self._Atoms_['LACS']))
        else:
            self._log_.warning("This dsf file has no coordinate pools (16-bit) defined!") 
        if '23OP' in self._Atoms_:
            self._extractPools_(32)
            self._extractScalings_(32)
            self._scaleV_(32, False) #False that scaling is not reversed
            self._updateProgress_(len(self._Atoms_['23CS']))
        else:
            self._log_.info("This dsf file has no 32-bit pools.")
        if 'SDMC' in self._Atoms_:
            self._unpackCMDS_()
            self._extractCMDS_()
        else:
            self._log_.warning("This dsf file has no commands defined.")
        return 0
        

    def _packAtoms_(self): #starts all functions to write all variables to strings (for later been written to file)
        self._log_.info("Preparing data to be written to file.")
        self._log_.info("This version only applies changes to POOL, CMDS and SCAL atoms (all other atoms inc. SC32 are written as read)!")
        self._encodeDefs_() ### TBD: Properties
        self._scaleV_(16, True) #de-scale again       ############### TBD: SC32 scaling !! ##############
        self._encodePools_()
        self._packAllScalings_()
        self._packCMDS_()
        return 0
    
    def _setPatchBoundary_(self, p): #sets min/max values of patch p          
        for t in p.trias:
            for v in t: #all vertexes of each triangle in patch
                if v[1] < 0: ############## NEW 3: ERROR CHECK - TO BE REMOVED LATER ##############
                    self._log_.warning("Index to Pool {} is negative: {}".format(v[0], v[1]))
                    return 1
                if v[1] >= len(self.V[v[0]]):  ############## NEW 3: ERROR CHECK - TO BE REMOVED LATER ##############
                    self._log_.warning("Index to Pool {} is higher {} than number of elements {}.".format(v[0], v[1], len(self.V[v[0]])))
                    return 1
                if self.V[v[0]][v[1]][0] < p.minx:
                    p.minx = self.V[v[0]][v[1]][0]
                if self.V[v[0]][v[1]][0] > p.maxx:
                    p.maxx = self.V[v[0]][v[1]][0]
                if self.V[v[0]][v[1]][1] < p.miny:
                    p.miny = self.V[v[0]][v[1]][1]
                if self.V[v[0]][v[1]][1] > p.maxy:
                    p.maxy = self.V[v[0]][v[1]][1]
        #####################self._log_.info("Patch Boundary x between {} - {} and y between {} - {}.".format(p.minx, p.maxx, p.miny, p.maxy))
        return 0
                                          
    
    def _addVertexToPool_(self, v, t, p): #adds Vertex v = [lon, lat] inside tria t to Pool p ################# NEW ############  
        if self.V[t[0][0]][t[0][1]][2] > -32768 and self.V[t[1][0]][t[1][1]][2] > -32768 and self.V[t[2][0]][t[2][1]][2] > -32768:
            #all vertices have height directly assigned, probably no raster excistent, so also height for new vertex has to be calculated
            l0, l1 = PointLocationInTria(v, self.TriaVertices(t)) # returns length for vectors from point t3 with l0*(t2-t0) and l1*(t2-t1)
            height = self.V[t[2][0]][t[2][1]][2] + l0 * (self.V[t[0][0]][t[0][1]][2] - self.V[t[2][0]][t[2][1]][2])  + l1 * (self.V[t[1][0]][t[1][1]][2] - self.V[t[2][0]][t[2][1]][2])
            #### Note: In case no raster is existent also the vertex normal has to be adapted (to be done in later step for all new vertices)
        else:
            height = -32768.0
        v_new = [ v[0], v[1], height, 0, 0 ] ## actually normal vectors are -1.5259021896696368e-05 instead of 0
        for i in range(5, len(self.V[p][0])): #go through optional coordinates s, t values based on first vertex in according Pool
            l0, l1 = PointLocationInTria(v, self.TriaVertices(t)) # returns length for vectors from point t3 with l0*(t2-t0) and l1*(t2-t1)
            v_new.append(self.V[t[2][0]][t[2][1]][i] + l0 * (self.V[t[0][0]][t[0][1]][i] - self.V[t[2][0]][t[2][1]][i])  + l1 * (self.V[t[1][0]][t[1][1]][i] - self.V[t[2][0]][t[2][1]][i]) )
        self.V[p].append(v_new) #append vertex to pool p
        #### TBD: Make sure Pool Index for pool p does not exceed/reach 2^16 --> Could be handled when pools are converted to binary #### TBD ###
        
 
    def _addSubmVertex_(self, v, scalbase, elevscal=0.05): #adds vertex v to a pool allowing submeter elevation at minimum for required elevation scaling, and scalbase for all other scalings    #### NEW3 #########
        ############## TBD:DEEPCOPY HERE in fucntion below FOR v AMD SCALBASE instead of doint this in calling funciton ############
        v[2] = self.getElevation(v[0], v[1], v[2]) #make sure that v[2] is real elevation and not using raster ##### NEW 6 #########
        if len(self.submP) == 0: #no submeter Pool in the set, then search first all pools to find some   
            self._log_.info("Checking file for existing submeter pools....")
            for i in range(len(self.Scalings)): # go through all pools and check if submeter and if add to submP
                if len(self.Scalings[i]) >= 5: #only vertex pools with elevationa and normal vectors are relevant
                    if self.Scalings[i][2][0] < 65535:  #scaling is x/65535, if x is less then 65535 we have submeter scaling
                        self.submP.add(i)
                        self._log_.info("File already has submeter pool with index {} and offset {}m and scaling {}m.".format(i, self.Scalings[i][2][1], self.Scalings[i][2][0]/65535))
        poolID4v = None #find pool with ID where v could be inserted
        self._log_.info("Checking if vertex {} with elevation scaling {}m and base scale {} fits to existing pool.".format(v, elevscal, scalbase))
        for i in self.submP:
            if len(scalbase) == len(self.Scalings[i]): #does size of scaling vector / number of planes for pool fit
                counter = 0
                for j in range(len(scalbase)): #check for each plane j if v could fit between minimum and maximum reach of scale
                    if v[j] >= self.Scalings[i][j][1] and v[j] <= self.Scalings[i][j][1] + self.Scalings[i][j][0]:
                        counter += 1
                    else:
                        break
                if len(scalbase) == counter:  #all values of v are within range of scale for pool i
                    if elevscal >= self.Scalings[i][2][1] / 65535: #test thate required scale level for elevation could be realized e.g. elevscal of 0.05m (65535*0.05=3276.75 scaling)
                        if len(self.V[i]) < 65534: #Vertex pool has space to add addtional vertex
                            poolID4v = i #existing pool id found fulfilling all requirements
                            self._log_.info("Existing with index {} and scale vector {} found.".format(poolID4v, self.Scalings[poolID4v]))
                            break
        if poolID4v != None: #existing pool was found
            self.V[poolID4v].append(v)
        else: #no existing pool fulfilling requirements was found, so new pool has to be created and added with v
            if len(self.V) >= 65535: #reached maximumb number of pools, no pool could be added any more
                self._log_.error("DSF File already has maximum number of point pools. Addtional pools required for change can not be added!!!")
                return None
            self.V.append([v])
            self.Scalings.append(scalbase)
            self.Scalings[-1][2][0] = 65535 * elevscal #new multiplier based on required scaling for elevation defined for new pool
            self.Scalings[-1][2][1] = int(-500 + int((v[2]+500)/(65535 * elevscal)) * (65535 * elevscal)) #offset for elevation of v   ######## -500m is deepest vaule that can be defined with this routine ####
            poolID4v = len(self.Scalings) - 1 #ID for the pool is the last one added
            self._log_.info("New pool with index {} and scaling {} added to insert vertex.".format(poolID4v, self.Scalings[poolID4v]))
            self.submP.add(poolID4v)
        return poolID4v, len(self.V[poolID4v])-1 #returning ID of Pool and index where v was added

        
    def cutEdges(self, v, w, accuracy = 10): #cuts all edges of trias that intersect with segment from vertex v to w, if existing vertices/edges are closer than accuracy in m, they will be used   ########### NEW #########
        ## Note: If edge is completely within a tria then nothing is cut. In that case the segment will be inserted by just inserting the vertices of the edge
        #### TBD: What happens if segment vw is overlapping an existing edge of the tria?? Function intersection will not return values so only endpoints of vw or other cutting points outside should be inserted and thus fine. But to be checked !!!! ### TBD ###
        #### TBD: function adds new vertices for each tria. Already added trias in neighbour vertex should be re-used !!!!!! ###########################
        cPs = [] #functions returns a list of tuples to vertices (pool ID, index) that are lying on intersection line v to w (either created or used within accuracy)
        l = self.PatchesInArea(*self.BoundingRectangle([v,w]))   ### * converts returned tuple in an argument list
        self._log_.info("Inserting segement from {} to {} into mesh, which is relevant in {} patches.".format(v, w, len(l)))
        if len(l) == 0:
            self._log_.error("No relevant patch for segment from {} to {} found. Nothing to cut!".format(v, w))
            return 0
        for p in l: #now go only through potential patches to check details
            new_trias = [] #trias to be added when cutting in patch
            old_trias = [] #old trias to be removed when cutting in patch
            for t in p.trias:
                tv = self.TriaVertices(t)
                iv = [] #list of intersection vertices between line vw and edges of triangle t, could be between 0 and 3
                for edge in range(3): # go through edges of tria by index numbers
                    cuttingPoint = intersection(v, w, tv[edge], tv[(edge+1)%3]) # modulo 3 returns to first vertex to close tria
                    if cuttingPoint:
                        existing_vertex_close = False
                        for i in [edge, (edge+1)%3]: #check if v is too close to vertex of tria ONLY FOR ADJECENT vertices on edge   ######## NEW 4 ############
                            if distance(tv[i], cuttingPoint) < accuracy:
                                if self._DEBUG_: self._log_.debug("   Cutting Point {} not inserted as too close to vertex of the triangle it is in.".format(cuttingPoint))
                                cPs.append((t[i][0], t[i][1])) ### Attention: adds all vertices of t that are within accuracy, not just one
                                existing_vertex_close = True
                        if not existing_vertex_close:
                                cuttingPoint = (cuttingPoint[0], cuttingPoint[1], edge) #store number of edge cut as third element
                                iv.append(cuttingPoint)
                if len(iv) == 2:
                    if iv[1][0] > iv[0][0] or (iv[1][0] == iv[0][0] and iv[1][1] > iv[0][1]):
                        iv[0], iv[1] = iv[1], iv[0] #make sure to always start with most western (or most southern if both have same west coordinate) cutting Point in order to always get same mesh for overlays
                    if self._DEBUG_: self._log_.debug("   Two intersections found at {}.".format(iv))
                    v_Pool = t[0][0] #select Pool that will be extended by new cutting vertices by first point of tria
                    self._addVertexToPool_(iv[0], t, v_Pool)
                    self._addVertexToPool_(iv[1], t, v_Pool)
                    v_Index = len(self.V[v_Pool]) - 1 #new vertex is at the end, so index of second inserted is length of Pool, for the one added before -1
                    cPs.append((v_Pool, v_Index - 1))
                    cPs.append((v_Pool, v_Index))
                    if self._DEBUG_: self._log_.debug("    Vertex added to Pool {} with index {} and coordinates {}.".format(v_Pool, v_Index, self.V[v_Pool][v_Index]))
                    if self._DEBUG_: self._log_.debug("    Vertex added to Pool {} with index {} and coordinates {}.".format(v_Pool, v_Index-1, self.V[v_Pool][v_Index-1]))
                    edge = iv[0][2] ## start with edge with west/southern cutting point having index v_Index-1 (was inserted first)
                    if iv[1][2] == (iv[0][2] + 1) % 3: #second cutting point lies on next edge
                        new_trias.append([ [v_Pool, v_Index-1], [t[(edge+1)%3][0], t[(edge+1)%3][1]], [v_Pool, v_Index] ])
                        new_trias.append([ [v_Pool, v_Index-1], [v_Pool, v_Index], [t[(edge+2)%3][0], t[(edge+2)%3][1]] ])
                        new_trias.append([ [v_Pool, v_Index-1],  [t[(edge+2)%3][0], t[(edge+2)%3][1]], [t[edge][0], t[edge][1]] ])
                    else: #second cutting point must lie on previous edge
                        new_trias.append([ [v_Pool, v_Index-1], [t[(edge+1)%3][0], t[(edge+1)%3][1]],  [t[(edge+2)%3][0], t[(edge+2)%3][1]] ])
                        new_trias.append([ [v_Pool, v_Index-1], [t[(edge+2)%3][0], t[(edge+2)%3][1]], [v_Pool, v_Index]  ])
                        new_trias.append([ [v_Pool, v_Index-1], [v_Pool, v_Index], [t[edge][0], t[edge][1]] ])
                    old_trias.append(t)                   
                elif len(iv) == 1:
                    edge = iv[0][2]
                    if self._DEBUG_: self._log_.debug("   One intersections found at {}.".format(iv))
                    v_Pool = t[0][0] #select Pool that will be extended by new cutting vertex by first point of tria
                    self._addVertexToPool_(iv[0], t, v_Pool)
                    v_Index = len(self.V[v_Pool]) - 1 #new vertex is at the end, so index is length of Pool
                    cPs.append((v_Pool, v_Index))
                    if self._DEBUG_: self._log_.debug("    Vertex added to Pool {} with index {} and coordinates {}.".format(v_Pool, v_Index, self.V[v_Pool][v_Index]))
                    new_trias.append([ [t[edge][0], t[edge][1]], [v_Pool, v_Index],  [t[(edge+2)%3][0], t[(edge+2)%3][1]] ])
                    new_trias.append([ [v_Pool, v_Index], [t[(edge+1)%3][0], t[(edge+1)%3][1]], [t[(edge+2)%3][0], t[(edge+2)%3][1]] ])
                    old_trias.append(t)
            for nt in new_trias:
                p.trias.append(nt)
                if self._DEBUG_: self._log_.debug("    Added triangle  {} to Patch {}.".format(nt, self.Patches.index(p)))
                if self._DEBUG_: self._log_.debug("         Having coordinates {}.".format(self.TriaVertices(nt)))
            for ot in old_trias:
                p.trias.remove(ot)
                if self._DEBUG_: self._log_.debug("    Removed triangle  {}.".format(ot))
                if self._DEBUG_: self._log_.debug("         Having coordinates {}.".format(self.TriaVertices(ot)))
            p.trias2cmds() ##Update Commands accordingly
        return cPs 
                     
                     
    def insertVertex(self, v, accuracy=10): #adds new vertex to dsf and updates mesh accordingly, returns exsting verticy when it is closer than accuracy in m, returns tuple of Pool, Index as set   ########### NEW ###########
        l = []
        iPs = [] #functions returns a list of tuples to vertices (pool ID, index) that have been inserted  #####NEW3 
        for p in self.Patches:
            if p.minx <= v[0] <= p.maxx and p.miny <= v[1] <= p.maxy:
                l.append(p)
        self._log_.info("Inserting Vertex {} in mesh, which is relevant in {} patches.".format(v, len(l)))
        if len(l) == 0:
            self._log_.error("Vertex {} cannot be inserted as there is no valid Patch found. Check coordinates!".format(v))
            return 0
        for p in l: #now go only through potential patches to check details
            for t in p.trias:
                tv = self.TriaVertices(t)
                if isPointInTria(v, tv):
                    for i in range(3): #check if no vertex of triangle is close to v
                        if distance(tv[i], v) < accuracy: 
                            self._log_.info("Vertex {} not inserted as too close to vertex of the triangle it is in.".format(v))
                            return [ (t[i][0], t[i][1]) ] ########### PPROBLEM: ONLY RETURNS ONE VERTEX and not all others from other terrain patches --> might be issue for insert RWYprofile in mesh
                    for e in range(3): #go through edges of tria and check if vertex v is too close to one of them
                        ortho = [-(tv[e][1] - tv[(e+1)%3][1]), (tv[e][0] - tv[(e+1)%3][0])] #ortogonal vector to edge
                        le, lo = PointLocationInTria(v, [tv[(e+1)%3], ortho, tv[e]]) #getting to now how far to go the edge to come to closest point le*e on the edge to v, lo is how far to go the ortho to reach v
                        closestPoint = [tv[e][0] + le * (tv[(e+1)%3][0] - tv[e][0]), tv[e][1] + le * (tv[(e+1)%3][1] - tv[e][1])]
                        if distance(v, closestPoint) < accuracy:  
                            self._log_.info("Vertex {} is too close to one edge of the triangle it is in --> Edge will be cut at {} and cutting point will be inserted.".format(v, closestPoint))
                            return self.cutEdges(v, [v[0] - lo * 1.001 * (ortho[0] - tv[e][0]),  v[1] - lo * 1.001 * (ortho[1] - v[1])], accuracy) #### NEW2 REDUCED fro,  * 1.000001 to * 1.001 toreally cross edge with ortho vector in - direction ### NEW 2
                            #### PROBLEM: reflects just one Edge it is close to but it could be close to two edges far away from vertex in a very sharp triangle !!!
                            ##### TBD: make multiplier * 1.001 variable dependent how far the point to be inserted is from the edge!!!! ##############
                            ################ BETTER ALTERNATIVE: INSTEAD OF CUTTING create new function that splits edge at given distance !!!!! ##########################
                    #### TBD: Make sure Pool Index does not exceed/reach 2^16 --> Could be handled when pools are converted to binary
                    #self._log_.info("Vertex {} lies in tria {} of patch number {}.".format(v, self.TriaVertices(t), self.Patches.index(p)))
                    #self._log_.info("    This tria is defined by vertices in pool/index {}.".format(t))
                    #self._log_.info("    First vertex of tria is {}.".format(self.V[t[0][0]][t[0][1]]))
                    v_Pool = t[0][0] #selct pool for new vertex based on first vertex of triangle
                    self._addVertexToPool_(v, t, v_Pool) 
                    v_Index = len(self.V[v_Pool]) - 1 #new vertex is at the end, so index is length of Pool
                    iPs.append((v_Pool, v_Index)) #NEW3
                    if self._DEBUG_: self._log_.debug("    Vertex added to Pool {} with index {} and coordinates {}.".format(v_Pool, v_Index, self.V[v_Pool][v_Index]))
                    p.trias.append([ [v_Pool, v_Index], [t[0][0], t[0][1]], [t[1][0], t[1][1]] ])
                    if self._DEBUG_: self._log_.debug("    Added triangle {}.".format(p.trias[-1]))
                    p.trias.append([ [v_Pool, v_Index], [t[1][0], t[1][1]], [t[2][0], t[2][1]] ])
                    if self._DEBUG_: self._log_.debug("    Added triangle {}.".format(p.trias[-1]))
                    p.trias.append([ [v_Pool, v_Index], [t[2][0], t[2][1]], [t[0][0], t[0][1]] ])
                    if self._DEBUG_: self._log_.debug("    Added triangle {}.".format(p.trias[-1]))
                    p.trias.remove(t)
                    if self._DEBUG_: self._log_.debug("    Removed triangle {}.".format(t))
                    break #in same patch there will be no further triangle where v is in
            p.trias2cmds() #Update Commands accordingly
        return iPs
    
    def cutPolyInMesh(self, poly, accuracy=5): ### NEW #### #updates dsf mesh to contain the polygon shape as edges (if existing vertices, edges are closer than accourcy in m they will be used), elevation will be given by mesh  ######## NEW #####
        self._log_.info("Cutting Polygon with {} vertices into the mesh.".format(len(poly) - 1))
        border_vertices = set() #will inclued all vertices of new and existing vertices defining the border of the poly within the mesh
        for i in range(len(poly) - 1):   # assumes that last vertex in poly is same as first, to have a closed polygon   #### NEW 4 insert verices first ####
            for v in self.insertVertex(poly[i], accuracy): 
                border_vertices.add(v)
        for i in range(len(poly) - 1):   # assumes that last vertex in poly is same as first, to have a closed polygon   #### NEW 4 and thenn cut edges ###
            for v in self.cutEdges(poly[i], poly[i+1], accuracy):
                border_vertices.add(v)               
        inner_vertices = set()
        l = self.PatchesInArea(*self.BoundingRectangle(poly))
        for p in l:
            for t in p.trias:
                vt = self.TriaVertices(t)
                for i in range(3):
                    if PointInPoly(vt[i], poly):
                        inner_vertices.add((t[i][0], t[i][1]))
        self._log_.info("Cutting Polygon done. Now there are {} vertices inside and at the border of the Polygon in the mesh.".format(len(border_vertices.union(inner_vertices))))
        return border_vertices.union(inner_vertices)
        # retruned vertices have to be set to required elevation
        ## TBD normal vectors of these vertices need to be updated if no raster is existent / used
        
    def cutTerrainInMesh(self, poly, terrain, accuracy=5):  #cuts Poly but only keeps physical trias with new terrain replaced, puts trias in submeter pools    ##### NEW 8 #####
        self._log_.info("Cutting Polygon with {} vertices into the mesh and only keep physical trias with terrain: {}.".format(len(poly) - 1, terrain))
        V = self.cutPolyInMesh(poly, accuracy) #cuts poly and returns in V all vertices at border and inside poly
        self._log_.info("From pure cutting {} vertices have been returend inside and at border of poly.".format(len(V)))
        #Now identifiy or create terrain id for trias  + create according patch for those trias   #### TBD: Write as funciton to be used alsof from cutRwyProfileInMesh() ####
        terrain_id = 0 #find terrain id for new patch that includes trias for runway
        while terrain_id < len(self.DefTerrains) and self.DefTerrains[terrain_id] != terrain:
            self._log_.info("      Terrain id: {} named:{}.".format(terrain_id, self.DefTerrains[terrain_id]))  ### TESTING #####
            terrain_id += 1
        if terrain_id == len(self.DefTerrains):
            self._log_.info("Terrain {} for runway profile not in current list. Will be added with id {}!".format(terrain, terrain_id))
            self.DefTerrains[terrain_id] = terrain
        else:
            self._log_.info("Terrain {} for runway profile found in current list with id: {}.".format(self.DefTerrains[terrain_id], terrain_id))
        newTerrainPatch = None #will be defined below based on pool of first vertex inserted
        #Now remove trias from poly and insert physical ones in new patch
        Vcoords = [] #get all coordinates of vertices in cutted poly
        for v in V:
            Vcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))
        l = self.PatchesInArea(*self.BoundingRectangle(Vcoords)) 
        keptV = set() #set with all vertices which will remain
        for p in l:
            self._log_.info("Checking patch {} for terrain replacement with flag {}.".format(self.Patches.index(p), p.flag))
            toberemoved = []
            for t in p.trias:
                if (t[0][0], t[0][1]) in V and (t[1][0], t[1][1]) in V and (t[2][0], t[2][1]) in V:
                    toberemoved.append(t)
                    if p.flag == 1: #patch is made from physical trias
                        if newTerrainPatch == None: #Now patch needs to be defined
                            newTerrainPatch = XPLNEpatch(1, 0.0, -1.0, t[0][0], terrain_id) ##use pool index of first vertex as pool defintion
                            newTerrainPatch.cmds.append( [1, t[0][0] ] ) #selection of poolIndex has to be set always as first command ### TBD: to be optimized that not twice set when created and again in cmds
                        
                        newTria = [] #vertices for new tria    
                        for v in [ (t[0][0], t[0][1]), (t[1][0], t[1][1]), (t[2][0], t[2][1]) ]: #go through all vertices in tria   #### avoid deepcopy by excplicit value listing, neccessary ???? ###
                            #just use first five values of vertex and scaling and create a new submeter vertex in dsf   ##### NEW 9 #####
                            newSubmV = []
                            RangeForNewSubmV = []
                            for i in range(5): #This is doing deepcopy for vertex list and scaling list    
                                newSubmV.append(self.V[v[0]][v[1]][i])
                                RangeForNewSubmV.append([self.Scalings[v[0]][i][0], self.Scalings[v[0]][i][1]])
                            pid, n = self._addSubmVertex_(newSubmV, RangeForNewSubmV)
                            keptV.add((pid, n)) #new Vertex with submPool and only 5 cooordinates will be kept
                            newTria.append((pid, n))
                        newTerrainPatch.trias.append(newTria)
                        if self._DEBUG_: self._log_.debug("    New triangle {} for new terrain.".format(newTria))
            for t in toberemoved:
                p.trias.remove(t) ### Remove all triangles in poly
                if self._DEBUG_: self._log_.debug("    Removed triangle {} for new terrain.".format(t))
            p.trias2cmds() #Update Commands accordingly
        newTerrainPatch.trias2cmds() #Update Commands for new terrain
        self._setPatchBoundary_(newTerrainPatch)
        self.Patches.append(newTerrainPatch)
        return keptV
        
    def cutRwyProfileInMesh(self, rwy, terrain, interval=20, accuracy=1):
        #cuts runway in mesh with intersections of runway every interval meter; allows to set granular elevation
        def add2profile(V, shoulderStart, shoulderV, dictV): #adds vertices in V to shoulder and dict    
            for v in V:
                shoulderV.add(v)
                vdist = round(distance(shoulderStart, self.V[v[0]][v[1]]), 3) #dista### NEW 6 ###  ############ IDEA: ROUND on more digits after . ?? ##############
                if len(self.V[v[0]][v[1]]) == 5: #priorize vertices without additional s/t coordinates ### NEW 6 ###
                    dictV[vdist] = v
                elif vdist not in dictV.keys():  #but if nothing was yet inserted in dict for distance insert vertex with s/t instead of nothing   ### NEW 6 ###
                    dictV[vdist] = v
            for k in dictV.keys(): #now we have to adapt vertices in dict that still have s/t coordinates ### NEW 6 ### ######## !!!! NO need to go through complete dict, but just new insertions, at least for loginfo !!!!! #############
                v = dictV[k]
                if len(self.V[v[0]][v[1]]) > 5: #so if additional coordinates
                    #just use first five values of vertex and scaling and create a new submeter vertex in dsf
                    newSubmV = []
                    RangeForNewSubmV = []
                    for i in range(5): #This is doing deepcopy for vertex list and scaling list    ##### NEW 6 #####
                        newSubmV.append(self.V[v[0]][v[1]][i])
                        RangeForNewSubmV.append([self.Scalings[v[0]][i][0], self.Scalings[v[0]][i][1]])
                    pid, n = self._addSubmVertex_(newSubmV, RangeForNewSubmV)
                    dictV[k] = [pid, n] #and replace it in dictionary
                
        self._log_.info("Cutting Runway with boundary {} into the mesh, having intersecitons every {} m.".format(rwy, interval))
        if accuracy > 1.5: #don't allow too high accuracy as this causes errors in profile generation    #### NEW 7 ####
            self._log_.warning("Accuracy {}m is too high for profile generation. Reduced to 1m!".format(accuracy))
            accuracy = 1      
        verticesOnShoulder = [set(), set(), set(), set()] #list of sets, that will inclued all verteices of shoulder per side of RWY
        dictShoulder = [{}, {}, {}, {}] #list of dictonaries with key is distance from point RWY-side-start-Point and pool/vertex id to a vertex as value
        ### IMPORTANT: Order of lists are as order of rwy-corners: rwy button, long side A, rwy top, long side B
        
        #First we cut rectangle for runway into the mesh
        for i,start,end in [(1,1,2), (3,0,3), (0,0,1), (2,2,3)]: #i index for list per RWY, side; start and end are indexes to runway corners in rwy; first insert long sides then button, then top
            self._log_.info("Inserting runway side from corner {} to corner {} in mesh.".format(start, end))
            add2profile(self.cutEdges(rwy[start], rwy[end], accuracy), rwy[start], verticesOnShoulder[i], dictShoulder[i])
            add2profile(self.insertVertex(rwy[start], accuracy), rwy[start], verticesOnShoulder[i], dictShoulder[i]) ##### SHOULD BE INSERTED IN DICT WITH DISTANCE KEY 0
            add2profile(self.insertVertex(rwy[end], accuracy), rwy[start], verticesOnShoulder[i], dictShoulder[i])  ##### SHOULD BE INSERTED IN DICT WITH DISTANCE KEY rwy side length
            ##########Runway long sides may have different length --> always insert last node with same length to avoid differences ???? ##################
        
        #Second we insert vertices on long rwy sides to opposite sides to have symetric points on each side
        self._log_.info("Add to runway shoulder A the vertices of Shoulder B.")
        for k in dictShoulder[3]: #add all distances in Shoulder B (index 3) to Directory of Shoulder A (index 1), so that there are always oposite vertices
            add2profile(self.insertVertex([rwy[1][0] + k / distance(rwy[0], rwy[3]) * (rwy[2][0] - rwy[1][0]), rwy[1][1] + k / distance(rwy[0], rwy[3]) * (rwy[2][1] - rwy[1][1])], accuracy), rwy[1], verticesOnShoulder[1], dictShoulder[1])
        self._log_.info("Now vice versa. add to runway shoulder B the vertices of Shoulder A.")
        ########## IS THIS NO DOUBLE INSERTION NOW AS B WAS ALREADY ADDED TO A????? ------> YES!!!!!! ### TBD: Exclude insertions, where vertex is already close????
        for k in dictShoulder[1]: #add all distances in Shoulder A(index 1) to Directory of Shoulder B (index 3), so that there are always oposite vertices
            add2profile(self.insertVertex([rwy[0][0] + k / distance(rwy[1], rwy[2]) * (rwy[3][0] - rwy[0][0]), rwy[0][1] + k / distance(rwy[1], rwy[2]) * (rwy[3][1] - rwy[0][1])], accuracy), rwy[0], verticesOnShoulder[3], dictShoulder[3])
            ########## Distanc keys for directories are not identical for each side, as they are here calculated each time now from start rwy corner !!!!!! Okay???
             
        #Third we insert on each sides the vertices for the intervalls
        ################### TBD: After insertion on each side check if it was really inserted in dict (dict growing, if not add new entry one mm away)?? Could save verification in step 4?? ######
        ############ TBD: Check that point that point on side A inserted in B is not again inserted in A !!!!! #####################################
        self._log_.info("Now cut runway in intervals.")
        accuracy /= 2 #USE HALF OF ACCURACY, TO HAVE HIGHER CHANCE TO GET ALWAYS A VERTEX ON OPPOSITE INSERTED
        intervalposition = 0
        dictShoulderA = dict() # new dict for side A required to not loop over a manipultated list !!??  ### NEW 6 ###
        dictShoulderB = dict() # new dict for side B required to not loop over a manipultated list !!??  ### NEW 6 ###
        for k in sorted(dictShoulder[1].keys()): ### Shoulder A is index 1
            width = (k-intervalposition) / (int((k - intervalposition) / interval) + 1) #define new width for intervals with maximum length of interval to next point in dict
            intervalposition += width #skip first already included vertex ##TBD: check if could be done better
            while intervalposition < int(k):
                add2profile(self.insertVertex([rwy[1][0] + intervalposition / distance(rwy[1], rwy[2]) * (rwy[2][0] - rwy[1][0]), rwy[1][1] + intervalposition / distance(rwy[1], rwy[2]) * (rwy[2][1] - rwy[1][1])], accuracy), rwy[1], verticesOnShoulder[1], dictShoulderA)
                #rwyApv needs v to be added
                # to be done: rwyAcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))
                intervalposition += width
            dictShoulderA[k] = dictShoulder[1][k] #current vertex is now at position in dictonary
            ##rwyApv.append(v)
            ##rwyAcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))
            intervalposition = int(k)
        intervalposition = 0
        for k in sorted(dictShoulder[3].keys()): ### Shoulder B is index 3
            width = (k-intervalposition) / (int((k - intervalposition) / interval) + 1) #define new width for intervals with maximum length of interval to next point in dict
            intervalposition += width #skip first already included vertex ##TBD: check if could be done better
            while intervalposition < int(k):
                add2profile(self.insertVertex([rwy[0][0] + intervalposition / distance(rwy[0], rwy[3]) * (rwy[3][0] - rwy[0][0]), rwy[0][1] + intervalposition / distance(rwy[0], rwy[3]) * (rwy[3][1] - rwy[0][1])], accuracy), rwy[0], verticesOnShoulder[3], dictShoulderB)                                               
                #rwyBpv needs v to be added
                #to be done rwyBcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))
                intervalposition += width
            dictShoulderB[k] = dictShoulder[3][k] #current vertex is now at position in dictonary
            ##rwyBpv.append(v)
            ##rwyBcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))
            intervalposition = int(k)
           
        #Forth: Check if on each side are the same number of vertices, if not insert the missing ones on opposite site
        distancesA = sorted(dictShoulderA.keys()) ######### DEEPCOPY ????? #####
        distancesB = sorted(dictShoulderB.keys()) ######### DEEPCOPY ????? #####
        accuracy *= 2 ### NOW GO UP AGAIN FOR WANTED accuracy  OKAY???? 
        
        i = 0
        self._log_.info("DISTANCES-ARRAY-LENGTH A:{}   B: {}".format(len(distancesA), len(distancesB))) ################ JUST ERROR CHECKING #################
        while i < len(distancesA) and i < len(distancesB):
            self._log_.info("CurrentDistances A: {}  B: {}".format(distancesA[i], distancesB[i]))
            i += 1
        
        i, j = 0, 0
        while i < len(distancesA) or j < len(distancesB): ######## PROBLEM !!!!!!! Insertions will not be on required position, in case on other site vertex is close --> insert with accuracy 0.001 --> Done below!! ###################
            #### ===> for now use 0.001 instead of accuracy !!!!
            self._log_.info("IJ-STEP  i: {}   j: {}   lenDistA: {}   lenDistB: {}".format(i, j, len(distancesA), len(distancesB)))   ################ JUST ERROR CHECKING #################
            if i >= len(distancesA): #side A finished, but still vertices on side B
                self._log_.info("SHOULDER-WARNING A END: Points on shoulder A finished in distance: {} and B in distance: {} has still vertex to be inserted in A!".format(distancesA[-1], distancesB[j]))
                add2profile(self.insertVertex([rwy[1][0] + distancesB[i] / distance(rwy[1], rwy[2]) * (rwy[2][0] - rwy[1][0]), rwy[1][1] + distancesB[i] / distance(rwy[1], rwy[2]) * (rwy[2][1] - rwy[1][1])], 0.001), rwy[1], verticesOnShoulder[1], dictShoulderA)
                j += 1
            elif j >= len(distancesB):# side B finsished, but still vertices on side A
                self._log_.info("SHOULDER-WARNING B END: Points on shoulder B finished in distance: {} and A in distance: {} has still vertex to be inserted in B!".format(distancesA[i], distancesB[-1]))
                add2profile(self.insertVertex([rwy[0][0] + distancesA[i] / distance(rwy[0], rwy[3]) * (rwy[3][0] - rwy[0][0]), rwy[0][1] + distancesA[i] / distance(rwy[0], rwy[3]) * (rwy[3][1] - rwy[0][1])], 0.001), rwy[0], verticesOnShoulder[3], dictShoulderB)
                i += 1              
            elif distancesA[i] > distancesB[j] + accuracy: #side A has a too big distance jump compared with side B; insert also vertex again on side A with distance of side B
                self._log_.info("SHOULDER-WARNING A: Points on shoulder A in distance: {} and B in distance: {} are too far away (above accuracy: {}). Insert new vertex in shoulder A! ".format(distancesA[i], distancesB[j], accuracy))
                add2profile(self.insertVertex([rwy[1][0] + distancesB[j] / distance(rwy[1], rwy[2]) * (rwy[2][0] - rwy[1][0]), rwy[1][1] + distancesB[j] / distance(rwy[1], rwy[2]) * (rwy[2][1] - rwy[1][1])], 0.001), rwy[1], verticesOnShoulder[1], dictShoulderA)
                j += 1
            elif distancesA[i] + accuracy < distancesB[j]: #side B has a too big distance jump compared with side A; insert also vertex again on side B with distance of side A
                self._log_.info("SHOULDER-WARNING B: Points on shoulder A in distance: {} and B in distance: {} are too far away (above accuracy: {}). Insert new vertex in shoulder B! ".format(distancesA[i], distancesB[j], accuracy))
                add2profile(self.insertVertex([rwy[0][0] + distancesA[i] / distance(rwy[0], rwy[3]) * (rwy[3][0] - rwy[0][0]), rwy[0][1] + distancesA[i] / distance(rwy[0], rwy[3]) * (rwy[3][1] - rwy[0][1])], 0.001), rwy[0], verticesOnShoulder[3], dictShoulderB)
                i += 1
            else: #distances on A and B are okay
                self._log_.info("Points on shoulder A in distance: {} and B in distance: {} are at similar (below accuracy: {}).".format(distancesA[i], distancesB[j], accuracy))
                i += 1
                j += 1            
           
        #Fith: Created sorted list of vertices and coordinates on each long runway side
        rwyAcoords = [] #sorted coordinates of RWY border A from buttom to top of RWY
        rwyBcoords = [] #sorted coordinates of RWY border B from buttom to top of RWY
        rwyApv = [] #sorted list on pool/vertex ids to build new rwy trias for RWY border A
        rwyBpv = []  #sorted list on pool/vertex ids to build new rwy trias for RWY border B
        for k in sorted(dictShoulderA.keys()):
            v = dictShoulderA[k]
            rwyApv.append(v)
            rwyAcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))
        for k in sorted(dictShoulderB.keys()):
            v = dictShoulderB[k]
            rwyBpv.append(v)
            rwyBcoords.append((self.V[v[0]][v[1]][0], self.V[v[0]][v[1]][1]))           
        self._log_.info("CREATED RUNWAY PROFILE:")
        for i in range(len(rwyAcoords)):
            self._log_.info("Interval at {}m {} coords:{}    opposite {}m {} coords:{} ".format(round(distance(rwyAcoords[0], rwyAcoords[i]), 3), rwyApv[i], rwyAcoords[i], round(distance(rwyBcoords[0], rwyBcoords[i]), 3), rwyBpv[i], rwyBcoords[i]))
        
        #Sixth: Identify all vertices inside of runway and add to border verticse to be all removed below for new profile
        inner_vertices = set() # Find vertices inside rwy boundary
        l = self.PatchesInArea(*self.BoundingRectangle(rwy))
        for p in l:
            for t in p.trias:
                vt = self.TriaVertices(t)
                for i in range(3):
                    if PointInPoly(vt[i], rwy):
                        inner_vertices.add((t[i][0], t[i][1]))
        allRWYvertices = verticesOnShoulder[1].union(verticesOnShoulder[3]) #vertices on boundaries from long runway sides 
        allRWYvertices = allRWYvertices.union(verticesOnShoulder[0]) #also vertices inserted from button     ###### NEW 4 #######
        allRWYvertices = allRWYvertices.union(verticesOnShoulder[2]) #and vertices from top of runway        ###### NEW 4 #######
        allRWYvertices = allRWYvertices.union(inner_vertices)            ######## NEW 4: always assign result of union again to allRWYvertices !! ####################
        self._log_.info("ALL RWY VERTICES: {}.".format(allRWYvertices))
        
        #Seventh: Remove all trias inside of runway
        for p in l:
            toberemoved = []
            for t in p.trias:
                if self._DEBUG_: self._log_.debug("Checking if tria {} in RWY.".format(t))
                if (t[0][0], t[0][1]) in allRWYvertices and (t[1][0], t[1][1]) in allRWYvertices and (t[2][0], t[2][1]) in allRWYvertices:
                    toberemoved.append(t)
            for t in toberemoved:
                p.trias.remove(t) ### Remove all triangles of RWY
                if self._DEBUG_: self._log_.debug("    Removed RWY triangle {}.".format(t))
            p.trias2cmds() #Update Commands accordingly
            
        #Eigth: Identifiy or create terrain id for trias under runway building the profile + create according patch for those trias
        terrain_id = 0 #find terrain id for new patch that includes trias for runway
        while terrain_id < len(self.DefTerrains) and self.DefTerrains[terrain_id] != terrain:
            terrain_id += 1
        if terrain_id == len(self.DefTerrains):
            self._log_.info("Terrain {} for runway profile not in current list. Will be added with id {}!".format(terrain, terrain_id))
            self.DefTerrains[terrain_id] = terrain
        else:
            self._log_.info("Terrain {} for runway profile found in current list with id: {}.".format(self.DefTerrains[terrain_id], terrain_id))
        rwyPatch = XPLNEpatch(1, 0.0, -1.0, rwyApv[0][0], terrain_id) ##use pool index of first rwyAborder vertex as
        rwyPatch.cmds.append( [1, rwyApv[0][0] ] ) #selection of poolIndex has to be set always as first command ### TBD: to be optimized that not twice set when created and again in cmds
        
        #Nineth: Create new patch for trias under runway
        for i in range(1, len(rwyApv)): #go through all vertices on runway side A  (function above checked that same numbers of vertices are on both side)
            #IMPORTANT: insert vertices of trias in clockwise order --> was checked when creating box for runway
            rwyPatch.trias.append([ rwyApv[i-1], rwyApv[i], rwyBpv[i-1] ])
            rwyPatch.trias.append([ rwyBpv[i-1], rwyApv[i], rwyBpv[i] ])

        for i in range(len(rwyApv), len(rwyBpv)):#add trias in case on side B are still more vertices than on side A ######### NEW 0.4.3 ##########
            self._log_.info("RWY PROFILE ON SIDE B HAS ADDITIONAL VERTEX: at {}m {} coords: {}. Additional tria for that vertex added.".format(round(distance(rwyBcoords[0], rwyBcoords[i]), 3), rwyBpv[i], rwyBcoords[i]))
            rwyPatch.trias.append([ rwyBpv[i-1], rwyApv[-1], rwyBpv[i] ])

        #Now insert vertices on buttom in first tria of runway 
        self._log_.info("INSERT VERTICES ON BUTTOM OF RUNWAY:")
        rwyPatch.trias.remove([rwyApv[0], rwyApv[1], rwyBpv[0]]) #remove first tria which is replaced below by new inserted trias for button vertices
        previousVonShoulder = rwyBpv[0]
        for k in sorted(dictShoulder[0].keys()):
            if k >= 0.01: #overjumps first corner of rwy on buttom, this is already in previousVonShoulder
                rwyPatch.trias.append([ previousVonShoulder, dictShoulder[0][k], rwyApv[1] ])
                previousVonShoulder = dictShoulder[0][k]
                #assumes that last k is key to second corner of runway button
        #Now insert vertices on top in last tria of runway 
        self._log_.info("INSERT VERTICES ON TOP OF RUNWAY:")
        if [rwyBpv[-2], rwyApv[-1], rwyBpv[-1]] in rwyPatch.trias: #check that top tria really exists an give error if not ## NEW 4.3 
            rwyPatch.trias.remove([rwyBpv[-2], rwyApv[-1], rwyBpv[-1]]) #remove last tria which is replaced below by new inserted trias for top vertices
        else:
            self._log_.error("Top triangle of RUNWAY not existing! Might result in erroroes mesh!")
        previousVonShoulder = rwyApv[-1]
        for k in sorted(dictShoulder[2].keys()):
            if k >= 0.01: #overjumps first corner of rwy on buttom, this is already in previousVonShoulder
                rwyPatch.trias.append([ previousVonShoulder, dictShoulder[2][k], rwyBpv[-2] ])
                previousVonShoulder = dictShoulder[2][k]
                #assumes that last k is key to second corner of runway buton                
        rwyPatch.trias2cmds()
        self._setPatchBoundary_(rwyPatch)
        self.Patches.append(rwyPatch)
        return verticesOnShoulder[1].union(verticesOnShoulder[3], verticesOnShoulder[0], verticesOnShoulder[2]) #also include vertices of buttom/top to rwy vertices    #### NEW 4 ######
    ########### END of cutRunwayProfileInMesh ##################################################################
      
    def getElevation(self, x, y, z = -32768): #gets Elevation at point (x,y) from raster grid Elevation
        if int(z) != -32768: #in case z vertex is different from -32768 than this is the right height and not taken from raster
            return z
        if "sim/west" not in self.Properties:
            self._log_.error("Cannot get elevation as properties like sim/west not defined!!!")
            return None
        if not (int(self.Properties["sim/west"]) <= x <= int(self.Properties["sim/east"])):
            self._log_.error("Cannot get elevation as x coordinate is not within boundaries!!!")
            return None
        if not (int(self.Properties["sim/south"]) <= y <= int(self.Properties["sim/north"])):
            self._log_.error("Cannot get elevation as y coordinate is not within boundaries!!!")
            return None
        if len(self.DefRasters) == 0: #No raster defined, use elevation from trias  #### NEW 7 ####
            self._log_.info("dsf includes no raster, elevation taken from trias!")  #### ERROR CHECKING - TO BE REMOVED ###
            l = self.PatchesInArea(y, y, x, x)  #find all patches relevant for point (x, y)
            for p in l:
                for t in p.trias:
                    tv = self.TriaVertices(t)
                    if isPointInTria([x, y], tv):
                        l0, l1 = PointLocationInTria([x, y], tv)
                        elev = self.V[t[2][0]][t[2][1]][2] + l0 * (self.V[t[0][0]][t[0][1]][2] - self.V[t[2][0]][t[2][1]][2])  + l1 * (self.V[t[1][0]][t[1][1]][2] - self.V[t[2][0]][t[2][1]][2]) #elvation is plane 2
                        self._log_.info("Elevation at x: {}  y: {}  is  {}m.".format(x, y, elev))  #### ERROR CHECKING - TO BE REMOVED ###
                        return elev
            self._log_.error("No triangle in mesh for position x: {}  y: {}  found in order to retrieve elevation!".format(x, y))  
        else: # use raster to get elevation; THIS VERSION IS ASSUMING THAT ELEVATION RASTER IS THE FIRST RASTER-LAYER (index 0), if it is not called "elevation" a warning is raised ###
            if self.DefRasters[0] != "elevation":
                self._log_.warning("Warning: The first raster layer is not called elevation, but used to determine elevation!")
            x = abs(x - int(self.Properties["sim/west"])) * (self.Raster[0].width - 1) # -1 from widht required, because pixels cover also boundaries of dsf lon/lat grid
            y = abs(y - int(self.Properties["sim/south"])) * (self.Raster[0].height - 1) # -1 from height required, because pixels cover also boundaries of dsf lon/lat grid
            if self.Raster[0].flags & 4: #when bit 4 is set, then the data is stored post-centric, meaning the center of the pixel lies on the dsf-boundaries, rounding should apply
                x = round(x, 0)
                y = round(y, 0)
            x = int(x) #for point-centric, the outer edges of the pixels lie on the boundary of dsf, and just cutting to int should be right
            y = int(y)
            return self.Raster[0].data[x][y]   ############## REMOVED () around returned value here !!!!!! ######################


    def getPolys(self, type): #returns all polygons of one type (numbered as in DefPolys) in a list and for each poly parameter following all vertices as reference [poolId, index]
        l = [] #list of polygons to be returned
        for p in self.Polygons[type]:
            l.append([p[2]]) #start a new poly which starts with its parameter
            if p[1] == 13: #POLYGON RANGE
                for i in range(p[3], p[4]):
                    l[-1].append([p[0], i])
        ####### TO INCLUDE ALSO OTHER POLYGON COMMAND IDs between 12 and 15 !!!!!!!!!!!!!!!!!!!!!!!!!!
        return l

    def BoundingRectangle(self, vertices): #returns 4-tuple of (latS, latN, lonW, lonE) building the smallest rectangle to include all vertices in list as pairs of [lon, lat]
        minx = 181  #use maximal out of bound values to be reset by real coordinates from patch
        maxx = -181
        miny = 91
        maxy = -91
        for v in vertices: #all vertexes of each triangle in patch
            if v[0] < minx:
                minx = v[0]
            if v[0] > maxx:
                maxx = v[0]
            if v[1] < miny:
                miny = v[1]
            if v[1] > maxy:
                maxy = v[1]
        return miny, maxy, minx, maxx


    def TriaVertices(self, t): #returns 3 vertices of triangle as list of [lon, lat] pairs  ############# TBD: CHECK WHERE NEEDED !!! ################
        return [  [self.V[t[0][0]][t[0][1]][0], self.V[t[0][0]][t[0][1]][1]], [self.V [t[1][0]] [t[1][1]][0], self.V[t[1][0]][t[1][1]][1]],  [self.V[t[2][0]][t[2][1]][0], self.V[t[2][0]][t[2][1]][1]] ]

    
    def PatchesInArea (self, latS, latN, lonW, lonE):  #################### NEW UPDATE: Function uses now trias list of trias instead trias() build funciton  ### NEW #######################
    #
    # returns a list of all patches in dsf, where the rectangle bounding of the patch intersects
    # the rectangle area defined by coordinates
    # Note: it could be the case that there is not part of the patch really intersecting the area but this
    #       functions reduces the amount of patches that have then carefully to be inspected e.g. for real intersections
    #
        self._log_.info("Start to find patches intersecting area SW ({}, {}) and NE ({}, {}).".format(latS, lonW, latN, lonE))
        l = [] # list of patch-ids intersecting area
        count = 0
        for p in self.Patches:
            ## v = []    #### NEW: Not needed any more ###### NEW ###
            ## for t in p.triangles(): #all triangles in each patch
            ##    v.extend(self.TriaVertices(t)) #so all vertices of the patch
            ## miny, maxy, minx, maxx = self.BoundingRectangle(v)
            if self._DEBUG_: self._log_.debug("Checking patch {} of {} which lies in SW ({}, {}) and NE ({}, {}).".format(count, len(self.Patches), miny, minx, maxy, maxx))
            if not (p.minx < lonW and p.maxx < lonW): #x-range of box is not completeley West of area      #### NEW ### use p.minx instead of minx - also below ##### NEW ###
                if not (p.minx > lonE and p.maxx > lonE): #x-range of box is not completele East of area
                    if not (p.miny < latS and p.maxy < latS): #y-range is not completele South of area
                        if not (p.miny > latN and p.maxy > latN): #y-range is not conmpletele North of ares
                            l.append(p)  #so we have an intersection of box with area and append the patch index
                            if self._DEBUG_: self._log_.debug("Patch {} in SW ({}, {}) and NE ({}, {}) does intersect.".format(count, miny, minx, maxy, maxx))
            count += 1
        self._log_.info("{} intersecting patches of {} patches found.".format(len(l), len(self.Patches)))
        return l
           
                    
    def read(self, file):  
        self.__init__("_keep_logger_","_keep_statusfunction_") #make sure all values are initialized again in case additional read
        if not path.isfile(file):
            self._log_.error("File does not exist!".format(file))
            return 1
        flength = stat(file).st_size #length of dsf-file   
        self._progress_ = [0, 0, flength] #initilize progress start for reading
        with open(file, "rb") as f:    ##Open Tile as binary fily for reading
            self._log_.info("Opend file {} with {} bytes.".format(file, flength))
            start = f.read(12)
            if start.startswith(b'7z\xBC\xAF\x27\x1C'):
                if PY7ZLIBINSTALLED:
                    f.seek(0)
                    archive = py7zlib.Archive7z(f)
                    filedata = archive.getmember(archive.getnames()[0]).read()
                    self._log_.info("File is 7Zip archive. Extracted and read file {} from archive with decompressed length {}.".format(archive.getnames()[0], len(filedata)))
                    f.close()
                    f = BytesIO(filedata)
                    self._progress_ = [0, 0, len(filedata)] #set progress maximum to decompressed length
                    flength = len(filedata) #also update to decompressed length
                    start = f.read(12)
                else:
                    self._log_.error("File is 7Zip encoded! py7zlib not installed to decode.")
                    return 2
            identifier, version = unpack('<8sI',start)
            if identifier.decode("utf-8") != "XPLNEDSF" or version != 1:
                self._log_.error("File is no X-Plane dsf-file Version 1 !!!")  
                return 3            
            while f.tell() < flength - 16: #read chunks until reaching last 16 bytes hash value
                bytes = f.read(8)
                atomID, atomLength = unpack('<4sI',bytes)        
                atomID = atomID.decode("utf-8")
                if atomID in self._AtomStructure_.keys():
                    if self._DEBUG_: self._log_.debug("Reading top-level atom {} with length of {} bytes.".format(atomID, atomLength))
                else:
                    if self._DEBUG_: self._log_.debug("Reading atom {} with length of {} bytes.".format(atomID, atomLength))
                if atomID in self._AtomOfAtoms_:  
                    self._Atoms_[atomID] = [] #just keep notice in dictonary that atom of atoms was read
                elif atomID in self._AtomList_:
                    bytes = f.read(atomLength-8)   ##Length includes 8 bytes header
                    if atomID in self._Atoms_: #subatom already exists in dictionary
                        self._Atoms_[atomID].append(bytes) #append string to existing list
                    else:
                        if atomID in self._MultiAtoms_:
                            self._Atoms_[atomID] = [bytes] #create new list entry, as for multiple atoms more can follow to be appended
                        else:
                            self._Atoms_[atomID] = bytes #for single atoms there is just this string
                else:
                    self._log_.warning("Jumping over unknown Atom ID (reversed): {} with length {}!!".format(atomID, atomLength))
                    bytes = f.read(atomLength-8)
                    #return 4               
            if self._DEBUG_: self._log_.debug("Reached FOOTER with Hash-Value: {}".format(f.read(16)))
        self._log_.info("Finished pure file reading.")
        self._unpackAtoms_()
        return 0 #file successfull read

    
    def write(self, file): #writes data to dsf file with according file-name
        self._progress_[0] = 0
        self._progress_[1] = 0 #keep original file length as goal to reach in progress[2]
        self._packAtoms_() #first write values of Atom strings that below will written to file   
        m = md5() #m will at the end contain the new md5 checksum of all data in file
        with open(file, "w+b") as f:    ##Open Tile as binary fily for writing and allow overwriting of existing file
            self._log_.info("Write now DSF in file: {}".format(file)  )
            s = pack('<8sI', 'XPLNEDSF'.encode("utf-8"),1)  #file header
            m.update(s)
            f.write(s)
            for k in self._Atoms_.keys():
                if k in self._AtomOfAtoms_:
                    if self._DEBUG_: self._log_.debug("Writing atom of atoms {} with length {} bytes.".format(k, self._TopAtomLength_(k)))
                    s = pack('<4sI', k.encode("utf-8"), self._TopAtomLength_(k)) # for AtomsOfAtoms just store header (id + length including sub-atoms)
                    m.update(s)
                    f.write(s)
                elif k in self._MultiAtoms_:
                    for a in self._Atoms_[k]:
                        if self._DEBUG_: self._log_.debug("Writing multi atom {} with length {} bytes.".format(k, len(a) + 8))
                        s = pack('<4sI', k.encode("utf-8"), len(a) + 8) + a # add 8 for atom header length (id+length)
                        m.update(s)
                        f.write(s)
                        self._updateProgress_(len(s))
                else: #just single instance atom with plane data
                        if k in self._AtomStructure_.keys():
                            if self._DEBUG_: self._log_.debug("Writing top-level atom {} with length {} bytes.".format(k, len(self._Atoms_[k]) + 8))
                        else:
                            if self._DEBUG_: self._log_.debug("Writing single atom {} with length {} bytes.".format(k, len(self._Atoms_[k]) + 8))
                        s = pack('<4sI', k.encode("utf-8"), len(self._Atoms_[k]) + 8) + self._Atoms_[k] # add 8 for atom header length (id+length)
                        m.update(s)
                        f.write(s) 
                        self._updateProgress_(len(s))
            if self._DEBUG_: self._log_.debug("New md5 value appended to file is: {}".format(m.digest()))
            f.write(m.digest())
        self._log_.info("Finshed writing dsf-file.")
        return 0

