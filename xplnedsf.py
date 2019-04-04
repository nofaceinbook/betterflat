#******************************************************************************
#
# xplnedsf.py        Version 0.13
# ---------------------------------------------------------
# Python module for reading and writing X_Plane DSF files.
#   (zipped DSF files have to be unzipped with 7-zip first!)
#
#
# WARNIG: This code is still under development and may still have some errors.
#         In case you use it be very careful!
#
# Copyright (C) 2019 by schmax (Max Schmidt)
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



import os #required to retrieve length of dsf-file
import struct #required for binary pack and unpack
import hashlib #required for md5 hash in dsf file footer


class XPLNEpatch:
    def __init__(self, flag, near, far, poolIndex, defIndex):
        self.flag = flag
        self.near = near
        self.far = far
        self.defIndex = defIndex
        self.cmds = []
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
                    l.append( [ [p, c[v-2]], [p, c[v-1]], [p, c[v]] ] )
            elif c[0] == 27: # PATCH TRIANGLE STRIP CROSS POOL
                for v in range (6, len(c), 2): #skip first values - including command id are not needed any more
                    l.append( [ [c[v-5], c[v-4]], [c[v-3], c[v-2]], [c[v-1], c[v]] ] )                                
            elif c[0] == 28:  # PATCH TRIANGLE STRIP RANGE
                for v in range(c[1], c[2] - 2): # last index has one added so -2 as we go 2 up and last in range is not counted
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
            else:
                print ("Error: Patch command not know/supported. Skipped it....") 
        return l


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
    def __init__(self):
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

        
    def _extractPools_(self, bit = 16, log = 1):
        if bit == 32:
            if log:
                print ("Start to unpack and extract", len(self._Atoms_['23OP']),  "pools (", bit,  "bit)...", flush = True)
            atomstring = self._Atoms_['23OP']
            V = self.V32
            ctype = "<L"
            size = 4 #bytes read per coordinate in 32 bit pool
        else: #assuming the standard 16bit Pool case
            if log:
                print ("Start to unpack and extract", len(self._Atoms_['LOOP']),  "pools (", bit,  "bit)...")
            atomstring = self._Atoms_['LOOP']
            V = self.V
            ctype = "<H"
            size = 2 #bytes read per coordinate in 16 bit pool
        for s in atomstring: #goes through all Pools read; string s has to be unpacked
            nArrays, nPlanes = struct.unpack('<IB', s[0:5])
            if log > 1:
                print ("Pool number",len(V) ,"has", nArrays, "Arrays (vertices) with", nPlanes, "Planes (coordinates per vertex)! ")
            V.append([]) #the current pool starts empty
            for i in range(nArrays): ## span up multi-dimensional array for the new pool of required size (number of vertices in pool)
                V[-1].append([])
            pos = 5 #position in string s
            for n in range(nPlanes):
                encType, = struct.unpack('<B', s[pos : pos + 1])
                pos += 1
                if log > 1:
                    print ("    Plane", n, "is encoded:", encType)
                if encType < 0 or encType > 3: #encoding not defined
                    print ("Error: Stopping reading pool because not known encoding of plane!!!")
                    return [] ##This means we return empty pool, which can be used to detect error        
                i = 0  #counts how many arrays = vertices have been read in plane n
                while i < nArrays:
                    if encType >= 2: #this means data is run-length encoded
                        runLength, = struct.unpack('<B', s[pos : pos + 1])
                        pos += 1
                    else: #no run-length encoding (not tested yet!!!!!)
                        runLength = 1 #just read single values until end of this plane  
                    if runLength > 127: #means the following value is repeated
                        v, = struct.unpack(ctype, s[pos : pos + size]) #repeated value
                        pos += size
                        runLength -= 128 #only value without 8th bit gives now the number of repetitions
                        while runLength > 0:
                            V[-1][i].append(v)
                            runLength -= 1
                            i += 1
                    else:
                        while runLength > 0: #now just reading individual values
                            v, = struct.unpack(ctype, s[pos : pos + size])
                            pos += size
                            V[-1][i].append(v)
                            runLength -= 1
                            i += 1
                if encType == 1 or encType == 3: #values are also stored differenced
                    for i in range (1, nArrays):  
                        V[-1][i][n] = (V[-1][i][n] + V[-1][i-1][n]) % 65536  #undo differntiation and modulo for wrapping two byte unsigned integer 


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
        
    
    def _encodePools_(self, log = 1): #overwrites current Pool atom with actual values of all vertices
        ###### FOR THE MOMENT ONLY WORKING FOR 16 BIT POOLS !!!!!!!!! ################################
        if log:
            print ("Start to encode all pools... (read values of vertices are overwritten)", flush = True)
        self._Atoms_['LOOP'] = [] #start new (future version also think of creating Pool atom in case it new dsf file will be created !!!!!!!!!!)
        ####### This version only stores pools in differentiated run-length encoding !!! #############
        ## Start with differentiation of all values ##
        for p in self.V: # go through all pools
            encpool = struct.pack("<IB",len(p),len(p[0])) ## start string of binary encoded pool number of arrays and number of planes (taken from first vertex)
            for n in range(len(p[0])): # go through all planes; number of planes just taken from first vertex in pool
                plane = [p[0][n]] #set plane to first value (for starts with second value to calculete differences)
                for i in range(1, len(p)): #go through all values of a plane for differntiation and append to current plane
                    plane.append((p[i][n] - p[i-1][n]) % 65536)  #Calculate difference to previous value AND take care of wrapping for two byte unsigned integer
                encpool += struct.pack('<B',3) #plane will be encoded differntiated + runlength
                struct.pack('<H',p[0][n])
                ## Now perform run-length encoding ##
                for rlpair in self._encodeRunLength_(plane):
                    encpool += struct.pack('<B', rlpair[0]) #encode runlength value
                    for v in rlpair[1]:
                        encpool += struct.pack('<H', v)
            self._Atoms_['LOOP'].append(encpool)
 

    def _extractScalings_(self, bit = 16, log = 1): #extract scaling values from atoms and writes them to according lists
        if log:
            print ("Start to unpack and extract all saclings of", bit, "bit pools...", flush = True)
        if bit == 32:
            atomstring = self._Atoms_['23CS']
            Scalings = self.Scal32
        else: #assuming the standard 16bit Pool case
            atomstring = self._Atoms_['LACS']
            Scalings = self.Scalings
        for s in atomstring: #goes through all Scalings read; string s has to be unpacked
            Scalings.append([]) #new scaling entry for next pool   
            for i in range(int(len(s) / 8) ):  #for each scaling tuple - 4 bytes multiplier and 4 bytes offset
                m, o = struct.unpack('<ff', s[ i*8 : i*8 + 8] )  #retrieve multiplier and offset
                Scalings[-1].append([m, o])

                
    def _scaleV_(self, bit = 16, reverse = False, log = 1): #applies scaling to the vertices stored in V and V32
        if log:
            if reverse:
                print("Start to de-scale all", bit, "bit pools...", flush = True)
            else:
                print ("Start to scale all", bit, "bit pools...", flush = True)
        if bit == 32:
            V = self.V32
            Scalings = self.Scal32
        else: #assuming the standard 16bit Pool case
            V = self.V
            Scalings = self.Scalings
            
        if len(V) != len(Scalings):
            print ("Error: Amount of Scale atoms does not equal amount of Pools!!")
            return 1
        for p in range(len(V)): #for all Pools
            if V[p] == []: ###There can exist empty pools that have to be skipped for scaling!!!
                if log:
                    print("Info: empty pool number", p, "not scaled!")
                break
            if len(V[p][0]) != len(Scalings[p]): #take first vertex as example to determine number of coordinate planes in current pool
                print ("Error: Amount of scale values for pool", p, "does not equal the number of coordinate planes!!!")
                return 2
            for n in range(len(Scalings[p])): #for all scale tuples for that pool = all coordinate planes in pool
                if log > 1:
                    print("Will now scale pool", p, "plane", n, "with multiplier:", Scalings[p][n][0], "and offset:", Scalings[p][n][1])                
                if float(Scalings[p][n][0]) == 0.0:
                    if log > 1:
                        print ("    No! Plane will not be scaled because scale is 0 !!!!")
                    break
                for v in range(len(V[p])): #for all vertices in current plane
                    if reverse: #de-scale vertices
                        V[p][v][n] = round((V[p][v][n] - Scalings[p][n][1]) * 65535 / Scalings[p][n][0])   #de-scale vertex v in pool p for plane n by subtracting offset and dividing by multiplyer
                    else:  #scale vertices
                        V[p][v][n] = (V[p][v][n] * Scalings[p][n][0] / 65535) + Scalings[p][n][1]  #scale vertex v in pool p for plane n with multiplyer and offset
               
                
    def _extractRaster_(self, log = 2):  #extracts alll rasters from atoms and stores them in list
        if log:
            print("Extracting", len(self._Atoms_['IMED'])," raster layers...", flush = True)
        if len(self._Atoms_['IMED']) != len(self._Atoms_['DMED']):
            print("Error: Number of raster info atoms not equal to number of raster data atoms!!!")
            return 1           
        for rn in range(len(self._Atoms_['IMED'])):
            R = XPLNEraster()
            R.ver, R.bpp, R.flags, R.width, R.height, R.scale, R.offset = struct.unpack('<BBHLLff', self._Atoms_['IMED'][rn])
            if log > 1:
                print ("Info of new raster layer:", R.ver, R.bpp, R.flags, R.width, R.height, R.scale, R.offset)
            if R.flags & 1:  #signed integers to be read
                if R.bpp == 1:
                    ctype = "<b"
                elif R.bpp == 2:  ##### this is the only case tested so far !!!!!!!
                    ctype = "<h"
                elif R.bpp == 4:
                    ctype = "<i"
                else:
                    print ("Error: not allowed bytes per pixel in Raster Definition!!!")
                    return 2
            elif R.flags & 2: #unsigned integers to be read
                if R.bpp == 1:
                    ctype = "<B"
                elif R.bpp == 2:
                    ctype = "<H"
                elif R.bpp == 4:
                    ctype = "<I"
                else:
                    print ("Error: not allowed bytes per pixel in Raster Definition!!!")
                    return 3
            elif not (R.flags & 1) and not (R.flags & 2): #neither first nor second bit set means that 4 byte float has to be read
                if R.bpp == 4:
                    ctype ="<f"
                else:
                    print ("Error: not allowed bytes per pixel in Raster Definition!!!")
                    return 4

            for x in range(0, R.bpp * R.width, R.bpp): #going x-wise from east to west just the bytes per pixes
                line = []
                for y in range(0, R.bpp * R.height * R.width, R.bpp * R.width): #going y-wise from south to north, always jumping over the width of each x-line
                    v, = struct.unpack(ctype, self._Atoms_['DMED'][rn][y + x : y + x + R.bpp]) #unpack R.bpp bytes at position x+y of raster data atom number rn
                    v = v * R.scale + R.offset # APPLYING SCALE + OFFSET 
                    line.append(v) #the pixel appended is to a line from south to north (y-line)
                R.data.append(line) #south to north lines are appended to each other
            self.Raster.append(R) #so raster list of list is returned to be indexed by [x][y]
   
    def _unpackCMDS_(self, log=1):
        i = 0 #position in CMDS String
        if log:
            print ("Start unpacking of Commands", flush = True)
        while i < len(self._Atoms_['SDMC']):
            id, = struct.unpack('<B', self._Atoms_['SDMC'][i : i+1])
            self.CMDS.append([id])
            i += 1
            if id in self._CMDStructure_:
                l = self._CMDStructLen_[id][0] #length of bytes to read
                if l > 0:
                    y = struct.unpack('<' + self._CMDStructure_[id][0], self._Atoms_['SDMC'][i : i + l])
                    self.CMDS[-1].extend(y)
                    i += l
                if len(self._CMDStructLen_[id]) == 3: #read command with variable length
                    l = self._CMDStructLen_[id][1] #length of repeating value n to read
                    n, = struct.unpack('<' + self._CMDStructure_[id][1], self._Atoms_['SDMC'][i : i + l]) #number n of repetitions
                    if id == 15:
                        n += 1 #id = 15 seems a special case that there is one index more than windings  ########??????????
                    i += l
                    l = self._CMDStructLen_[id][2] #length of repeated bytes
                    for j in range(n): 
                       y = struct.unpack('<' + self._CMDStructure_[id][2], self._Atoms_['SDMC'][i : i + l])
                       self.CMDS[-1].extend(y)
                       i += l
            else:
                if id == 14: #special double packed case, which is explicetly treated separate and has special format returning lists inside commands !!!
                             ###### not tested yet !!!!!!!! ########
                    parameter, windings = struct.unpack('<HB', self._Atoms_['SDMC'][i : i + 3])
                    i += 3
                    self.CMDS[-1].extend(parameter)
                    for w in range(windings):
                        indices, = struct.unpack('<B', self._Atoms_['SDMC'][i : i + 1])
                        i += 1
                        windinglist = []
                        for m in range(indices):
                            y, = struct.unpack('<H', self._Atoms_['SDMC'][i : i + 2])
                            i += 2
                            windinglist.extend(y)
                        self.CMDS[-1].append(windinglist)
                else: #command id not tretated here until now
                    print ("Warning: ID", id, "not in CMDStrucuture --> not stored in list CMDS!!")
                    self.CMDS[-1] #delete already written id
            if log > 2:    
                print("CMD id", self.CMDS[-1][0], ":", self.CMDS[-1][1:], "(string pos next:", i, ")" )
        if log:
            print (len(self.CMDS) , "commands haven been unpacked.")


    def _extractCMDS_(self, log = 2): # extract CMDS and stores it as Mesh-Patches, Polygons, ...

        for i in range(len(self.DefPolygons)):
            self.Polygons.append([]) #span list of empty lists for all defined poygon types
        for i in range(len(self.DefObjects)): ########## OBJECTS NOT TESTED YET !!!!!!!!!!!!!!!!!!!!!!!!!!######
            self.Objects.append([]) #span list of empty lists for all defined poygon types

        patchPoolIndex = None #poolIndex currently used in current patch; if different from current in CMDS then change command is written to cmds of patch

        flag_physical = None #1 if physical, 2 if overlay
        nearLOD = None  
        farLOD = None
        poolIndex = 0
        defIndex = 0
        subroadtype = 0
        junctionoffset = 0

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
            elif 9 <= c[0] <= 11: #Network Commands
                if self.Networks == []: #first network command, so start with first entry
                    self.Networks.append([subroadtype, junctionoffset, poolIndex])
                elif self.Networks[-1][0] != subroadtype or self.Networks[-1][1] != junctionoffset or self.Networks[-1][2] != poolIndex: #chang of relevant base settings
                    self.Networks.append([subroadtype, junctionoffset, poolIndex]) #sp new entry with new base-settings
                self.Networks[-1].extend(c) #append complete command to build this network part on current base settings
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
                if patchPoolIndex != poolIndex:
                    self.Patches[-1].cmds.append( [1, poolIndex] ) #change pool with patch via according command
                    patchPoolIndex = poolIndex #now within patch also current pool will be used
                if self.Patches[-1].defIndex != defIndex:
                    print("Error: defIndex changed within patch")
                    return(1)
                self.Patches[-1].cmds.append(c) 

        if log:
            print(len(self.Patches), "Patches extracted from Commands.")
            print(len(self.Polygons), "different Polygon types including there definitions extracted from Commands.")
            print(len(self.Objects), "different Objects with placements coordinates extracted from Commands.")
            print(len(self.Networks), "different Network subtypes extracted from Commands (could include double count).", flush = True)


    def _unpackAtoms_(self, log = 1): #starts all functions to unpack and extract data froms strings in Atoms
        if log:
            print("Extracting properties and definitions ...", flush = True)
        self._extractProps_()
        self._extractDefs_()
        if 'IMED' in self._Atoms_:
            self._extractRaster_(log)
        else:
            print("Info: This dsf file has no raster layers.")
        self._extractPools_(16, log)
        self._extractScalings_(16, log)
        self._scaleV_(16, False, log) #False that scaling is not reversed
        if '23OP' in self._Atoms_:
            self._extractPools_(32, log)
            self._extractScalings_(32, log)
            self._scaleV_(32, False, log) #False that scaling is not reversed
        else:
            print("Info: This dsf file has no 32-bit pools.")
        self._unpackCMDS_(log)
        self._extractCMDS_(log)
        return 0
        

    def _packAtoms_(self, log = 1): #starts all functions to write all variables to strings (for later been written to file)
        if log:
            print("Preparing data to be written to file.", flush = True)
            print("    Info: This version only applies changes to POOL atoms (all other atoms are written as read)!", flush = True)
        self._scaleV_(16, True, log) #de-scale again
        self._encodePools_(log)
        return 0

            
    def getElevation(self, x, y, z = -32768): #gets Elevation at point (x,y) from raster grid Elevation
        if int(z) != -32768: #in case z vertex is different from -32768 than this is the right height and not taken from raster
            return z
        if "sim/west" not in self.Properties:
            print("Error: Properties like sim/west not defined!!!")
            return None
        if not (int(self.Properties["sim/west"]) <= x <= int(self.Properties["sim/east"])):
            print ("Error: x coordinate is not within boundaries!!!")
            return None
        if not (int(self.Properties["sim/south"]) <= y <= int(self.Properties["sim/north"])):
            print ("Error: y coordinate is not within boundaries!!!")
            return None
        ### THIS VERSION IS ASSUMING THAT ELEVATION RASTER IS THE FIRST RASTER-LAYER (index 0), if it is not called "elevation" a warning is printed ###
        if self.DefRasters[0] != "elevation":
            print ("Warning: The first raster layer is not called elevation, but used to determine elevation!")
        x = abs(x - int(self.Properties["sim/west"])) * (self.Raster[0].width - 1) # -1 from widht required, because pixels cover also boundaries of dsf lon/lat grid
        y = abs(y - int(self.Properties["sim/south"])) * (self.Raster[0].height - 1) # -1 from height required, because pixels cover also boundaries of dsf lon/lat grid
        if self.Raster[0].flags & 4: #when bit 4 is set, then the data is stored post-centric, meaning the center of the pixel lies on the dsf-boundaries, rounding should apply
            x = round(x, 0)
            y = round(y, 0)
        x = int(x) #for point-centric, the outer edges of the pixels lie on the boundary of dsf, and just cutting to int should be right
        y = int(y)
        return (self.Raster[0].data[x][y])  


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


    def TriaVertices(self, t): #returns 3 vertices of triangle as list of [lon, lat] pairs
        return [  [self.V[t[0][0]][t[0][1]][0], self.V[t[0][0]][t[0][1]][1]], [self.V [t[1][0]] [t[1][1]][0], self.V[t[1][0]][t[1][1]][1]],  [self.V[t[2][0]][t[2][1]][0], self.V[t[2][0]][t[2][1]][1]] ]

    
    def PatchesInArea (self, latS, latN, lonW, lonE, log=1):
    #
    # returns a list of all patches in dsf, where the rectangle bounding of the patch intersects
    # the rectangle area defined by coordinates
    # Note: it could be the case that there is not part of the patch really intersecting the area but this
    #       functions reduces the amount of patches that have then carefully to be inspected e.g. for real intersections
    #
        if log:
            print("Start to find patches intersecting area SW (", latS, lonW, ") and NE (", latN, lonE, ") ...")
        l = [] # list of patch-ids intersecting area
        count = 0
        for p in self.Patches:
            v = []
            for t in p.triangles(): #all triangles in each patch
                v.extend(self.TriaVertices(t)) #so all vertices of the patch
            miny, maxy, minx, maxx = self.BoundingRectangle(v)
            if log > 2:
                print("Checking patch", count, "of", len(self.Patches), "which lies SW (", miny, minx, ") and NE (", maxy, maxx, ")" )
            if not (minx < lonW and maxx < lonW): #x-range of box is not completeley West of area
                if not (minx > lonE and maxx > lonE): #x-range of box is not completele East of area
                    if not (miny < latS and maxy < latS): #y-range is not completele South of area
                        if not (miny > latN and maxy > latN): #y-range is not conmpletele North of ares
                            l.append(p)  #so we have an intersection of box with area and append the patch index
                            if log > 1:
                                print ("Patch", count, "SW (", miny, minx, ") and NE (", maxy, maxx, ") does intersect.")
            count += 1
        if log:
            print(len(l), "intersecting patches of", len(self.Patches), "patches found.")
        return l
           
                    
    def read(self, file, log=1):   ###### NEXT STEP: Also read 7-ZIP FILES ###########
        if not os.path.isfile(file):
            print("Error: File", file, "does not exist!")
            return 1
        flength = os.stat(file).st_size #length of dsf-file   ## Error condition if file not existing --> return 1
        with open(file, "rb") as f:    ##Open Tile as binary fily for reading
            if log:
                print ("Opend file", file, "with", flength, "bytes.", flush = True)
            bytes = f.read(12)
            identifier, version = struct.unpack('<8sI',bytes)
            if identifier.decode("utf-8") != "XPLNEDSF" or version != 1:
                print ("File is no X-Plane dsf-file Version 1 !!!")  
                return 2            
            while f.tell() < flength - 16: #read chunks until reaching last 16 bytes hash value
                bytes = f.read(8)
                atomID, atomLength = struct.unpack('<4sI',bytes)        
                atomID = atomID.decode("utf-8")
                if log > 1 and atomID in self._AtomStructure_.keys():
                    print ("Reading top-level atom", atomID, "with length of", atomLength, "bytes.")
                elif log > 2:
                    print ("Reading atom", atomID, "with length of", atomLength, "bytes.")
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
                    print ("WARNING: Jumping over unknown Atom ID (reversed):", atomID, "with length", atomLength) ####SHOULD IN NEXT VERSION RAISE EXCEPTION
                    bytes = f.read(atomLength-8)
                    #return 3               
            if log > 1:
                print("Reached FOOTER with Hash-Value:",f.read(16), flush = True)
        self._unpackAtoms_(log)
        return 0 #file successfull read

    
    def write(self, file, log=1): #writes data to dsf file with according file-name
        self._packAtoms_(log) #first write values of Atom strings that below will written to file   
        m = hashlib.md5() #m will at the end contain the new md5 checksum of all data in file
        with open(file, "w+b") as f:    ##Open Tile as binary fily for writing and allow overwriting of existing file
            if log:
                print ("Write now DSF in file:", file, flush = True)  
            s = struct.pack('<8sI', 'XPLNEDSF'.encode("utf-8"),1)  #file header
            m.update(s)
            f.write(s)
            for k in self._Atoms_.keys():
                if k in self._AtomOfAtoms_:
                    if log > 1:
                        print("Writing atom of atoms:", k, "with length:", self._TopAtomLength_(k), "bytes.")
                    s = struct.pack('<4sI', k.encode("utf-8"), self._TopAtomLength_(k)) # for AtomsOfAtoms just store header (id + length including sub-atoms)
                    m.update(s)
                    f.write(s)
                elif k in self._MultiAtoms_:
                    for a in self._Atoms_[k]:
                        if log > 2:
                            print("Writing multi atom:", k, "with length:", len(a) + 8, "bytes.")
                        s = struct.pack('<4sI', k.encode("utf-8"), len(a) + 8) + a # add 8 for atom header length (id+length)
                        m.update(s)
                        f.write(s)
                else: #just single instance atom with plane data
                        if log > 1 and k in self._AtomStructure_.keys():
                            print("Writing top-level atom:", k, "with length:", len(self._Atoms_[k]) + 8, "bytes.")
                        elif log > 2:
                            print("Writing single atom:", k, "with length:", len(self._Atoms_[k]) + 8, "bytes.")
                        s = struct.pack('<4sI', k.encode("utf-8"), len(self._Atoms_[k]) + 8) + self._Atoms_[k] # add 8 for atom header length (id+length)
                        m.update(s)
                        f.write(s)                    
            if log > 1:
                print ("New md5 value appended to file is:", m.digest())
            f.write(m.digest())
        if log:
            print("  ...finshed writing dsf-file.", flush = True)
        return 0

