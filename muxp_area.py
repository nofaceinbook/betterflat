# -*- coding: utf-8 -*-
#******************************************************************************
#
# muxp_area.py    Version: 0.0.2 exp
#        
# ---------------------------------------------------------
# Python Class for adapting mesh in a given area of an XPLNEDSF
# Used by Mesh Updater X-Plane (muxp)
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

from logging import getLogger
from muxp_math import *
from copy import deepcopy
    
class muxpArea:
    
    def __init__(self, dsf, logname):
        self.dsf = dsf #dsf from which area is extracted
        self.log = getLogger(logname + "." + __name__) #logging based on pre-defined logname
        self.apatches = set() #set of patches that are relevant for that area
        self.atrias = [] #area as array, of triangles where
                  # [0 - 2] list of all coordinates of tria vertices including s/t ##### ONLY REFERENCE NO DEEPCOPY / BETTER JUST DEEPCOPY of x,y corrds ?????
                  # [3 - 5] list of pool and vertex id in pool in dsf -- they stay even, if coordinates are changed or trias split; allows reference to original trias/values/scaling
                  # [6] index to patch tria was in dsf

    def extractMeshArea(self, latS, latN, lonW, lonE):
        """
        Extracts an area [latS, latN, lonW, lonE] from mesh in self.dsf and stores it in self.trias
        Important: This function really removes trias from dsf.
                   This function will add trias to existing ones.
        """
        self.log.info("Extracting area latS: {}, latN: {}, lonW: {}, lonE:{} from dsf.".format(latS, latN, lonW, lonE))
        triaCount = 0 #counts all trias in dsf
        for p in self.dsf.Patches:
            trias = p.triangles()
            tremoved = [] #trias that should be removed
            for t in trias:
                # get for each tria the bounding rectangle in miny, maxy, minx, maxxx
                minx = min(self.dsf.V[t[0][0]][t[0][1]][0], self.dsf.V[t[1][0]][t[1][1]][0], self.dsf.V[t[2][0]][t[2][1]][0])
                maxx = max(self.dsf.V[t[0][0]][t[0][1]][0], self.dsf.V[t[1][0]][t[1][1]][0], self.dsf.V[t[2][0]][t[2][1]][0])
                miny = min(self.dsf.V[t[0][0]][t[0][1]][1], self.dsf.V[t[1][0]][t[1][1]][1], self.dsf.V[t[2][0]][t[2][1]][1])
                maxy = max(self.dsf.V[t[0][0]][t[0][1]][1], self.dsf.V[t[1][0]][t[1][1]][1], self.dsf.V[t[2][0]][t[2][1]][1])
                #now check if bounding rectangle of tria intersects with area
                if not (minx < lonW and maxx < lonW): #x-range of box is not completeley West of area     
                    if not (minx > lonE and maxx > lonE): #x-range of box is not completele East of area
                        if not (miny < latS and maxy < latS): #y-range is not completele South of area
                            if not (miny > latN and maxy > latN): #y-range is not conmpletele North of ares
                                self.atrias.append([ self.dsf.V[t[0][0]][t[0][1]], self.dsf.V[t[1][0]][t[1][1]], self.dsf.V[t[2][0]][t[2][1]] ])
                                self.atrias[-1].extend(t)  #so we have an intersection of tria box with area and append the tria
                                self.atrias[-1].append(self.dsf.Patches.index(p))
                                self.apatches.add(self.dsf.Patches.index(p))
                                tremoved.append(t)
                                self.log.debug("Tria {} with latS: {}  latN: {}  lonW: {} and lonE: {} in area.".format(t, miny, maxy, minx, maxx))
                triaCount += 1
            for t in tremoved:
                trias.remove(t)
            p.trias2cmds(trias) #updates dsf patch with trias not including removed ones
        self.log.info("  ... dsf has {} trias and {} trias from {} different patches are now in the extracted area.".format(triaCount, len(self.atrias), len(self.apatches)))
        return
    
    def cutEdges(self, v, w, accuracy = 10):
        """
        Cuts all edges of trias that intersect with segment from vertex v to w in self.atrias.
        If existing vertices/edges are closer than accuracy in m, they will be used.
        Function returns all points on the cutting line (except endpoints, if they do not cut).
        Note: If edge is completely within a tria then nothing is cut.
              In that case the segment has to be inserted, then by just inserting the vertices of the edge.
        """
        self.log.info("Cutting area with edge from: {} to: {}.".format(v, w))
        cPs = [] #functions returns a list of vertices that are lying on intersection line v to w (either created or used within accuracy)
        new_trias = [] #trias to be added when cutting in patch
        old_trias = [] #old trias to be removed when cutting in patch
        for t in self.atrias: #go through all trias in area
            iv = [] #list of intersection vertices between line vw and edges of triangle t, could be between 0 and 3
            for edge in range(3): # go through edges of tria by index numbers
                cuttingPoint = intersection(v, w, t[edge][0:2], t[(edge+1)%3][0:2]) # modulo 3 returns to first vertex to close tria
                if cuttingPoint:
                    existing_vertex_close = False
                    for i in [edge, (edge+1)%3]: #check if v is too close to vertex of tria ONLY FOR ADJECENT vertices on edge  
                        if distance(t[i][0:2], cuttingPoint) < accuracy:
                            self.log.debug("   Cutting Point {} not inserted as too close to vertex of the triangle it is in.".format(cuttingPoint))
                            cPs.append(t[i][0:2]) ### Attention: adds all coordinates of t that are within accuracy, not just one
                            existing_vertex_close = True
                    if not existing_vertex_close:
                            cuttingPoint = (cuttingPoint[0], cuttingPoint[1], edge) #store number of edge cut as third element
                            iv.append(cuttingPoint)
                            cPs.append((cuttingPoint[0], cuttingPoint[1]))
            if len(iv) == 2:
                if iv[1][0] > iv[0][0] or (iv[1][0] == iv[0][0] and iv[1][1] > iv[0][1]):
                    iv[0], iv[1] = iv[1], iv[0] #make sure to always start with most western (or most southern if both have same west coordinate) cutting Point in order to always get same mesh for overlays
                self.log.debug("   Two intersections found at {}.".format(iv))
                
                edge = iv[0][2] ## start with edge with west/southern cutting point having index v_Index-1 (was inserted first)
                if iv[1][2] == (iv[0][2] + 1) % 3: #second cutting point lies on next edge
                    new_v0 = createFullCoords(iv[0][0], iv[0][1], t)
                    new_v1 = createFullCoords(iv[1][0], iv[1][1], t)
                    new_trias.append([ new_v0, t[(edge+1)%3], new_v1, t[3], t[4], t[5], t[6]]) 
                    new_trias.append([ new_v0, new_v1, t[(edge+2)%3], t[3], t[4], t[5], t[6]])
                    new_trias.append([ new_v0, t[(edge+2)%3], t[edge], t[3], t[4], t[5], t[6]])
                else: #second cutting point must lie on previous edge
                    new_v0 = createFullCoords(iv[0][0], iv[0][1], t)
                    new_v1 = createFullCoords(iv[1][0], iv[1][1], t)
                    new_trias.append([ new_v0, t[(edge+1)%3], t[(edge+2)%3], t[3], t[4], t[5], t[6]])
                    new_trias.append([ new_v0, t[(edge+2)%3], new_v1, t[3], t[4], t[5], t[6]])
                    new_trias.append([ new_v0,  new_v1, t[edge], t[3], t[4], t[5], t[6]])
                old_trias.append(t)                   
            elif len(iv) == 1:
                edge = iv[0][2]
                self.log.debug("   One intersection found at {}.".format(iv))
                new_v0 = createFullCoords(iv[0][0], iv[0][1], t)
                new_trias.append([ t[edge], new_v0, t[(edge+2)%3], t[3], t[4], t[5], t[6]])
                new_trias.append([ new_v0, t[(edge+1)%3], t[(edge+2)%3], t[3], t[4], t[5], t[6]])
                old_trias.append(t)
        for nt in new_trias:
            self.atrias.append(nt) #update area trias by appending new trias
        for ot in old_trias:
            self.atrias.remove(ot) #update area trias by by removing trias that are replaced by new ones
        self.log.info("  ... this cut returned {} cutting points.".format(len(cPs)))
        return cPs 
    
    def createDSFVertices(self, elevscal=1):
        """
        For the area this function will find for all new vertices created in the area that
        have not yet assigned a vertex in dsf vertex pool an according vertex in the pool.
        It allows also submeter elevations when elevscal (value in meter) is below 1.
        When elevscal=1 either the raster reference -32768.0 is used or existing elevations per vertex (which might also allow submeter)
        Important: All trias require reference to trias they are inside in the orginal dsf. This is required to get the correct scalings for the area.
        Important 2: This function sets direct elevation for vertices in case elevscal is below 1.
        OPTION FOR FUTURE: use area limits for scaling lat/lon values.
        """
        ############## TBD: Also adapt the vertex normals (perhaps only in case dsf has no raster) ################
        newPools = [] # create new pools for the new vertices, these are the according indices  ##### TBD: check for existing pools first that might be used ########
        for t in self.atrias:
            for vt in range(3): #all vertices of tria t 
                count = 0 #counting how many coordinates are still the same as the referenced vertex in dsf
                for vti in range(len(t[vt])):
                    if t[vt][vti] != self.dsf.V[t[vt+3][0]][t[vt+3][1]][vti]:
                        self.log.info("Reference for vertex {} not correct any more. {} differs to {}!".format(self.atrias.index(t), t[vt][vti], self.dsf.V[t[vt+3][0]][t[vt+3][1]][vti]))
                        break
                    count += 1
                if len(t[vt]) != count: #the reference for this vertex to the pool is not correct any more, new vertex needs to be inserted in dsf
                    v = deepcopy(t[vt]) #v has now deepcopy of all coordinates of that vertex that is to be inserted in the dsf
                    poolID4v = None #Searching for the pool ID for vertex v
                    if elevscal < 1 and v[2] < -32765: #we want to get submeter elevations but have vertex referencing raster
                        v[2] = self.dsf.getVertexElevation(v[0], v[1], v[2]) #make sure that v[2] is real elevation and not using raster (with raster no submeter); only works because raster is there
                    for i in newPools:
                        if len(v) == len(self.dsf.Scalings[i]): #does size of scaling vector / number of planes for pool fit
                            counter = 0
                            for j in range(len(v)): #check for each plane j if v could fit between minimum and maximum reach of scale
                                if v[j] >= self.dsf.Scalings[i][j][1] and v[j] <= self.dsf.Scalings[i][j][1] + self.dsf.Scalings[i][j][0]:
                                    counter += 1
                                else:
                                    break
                            if len(v) == counter:  #all values of v are within range of scale for pool i
                                poolID4v = i #existing pool id found fulfilling all requirements                            
                                break
                    if poolID4v != None: #existing pool was found
                        matchfound = False
                        for ev in self.dsf.V[poolID4v]: #check for all existing vertices ev if they match v already
                            counter = 0
                            for i in range(len(v)): #check for all planes/coordinates whether they are nearly equal
                                if abs(ev[i] - v[i]) >= self.dsf.Scalings[poolID4v][i][0] / 65535: #if difference is lower than scale multiplier both coordinates would end up same after endcoding, so they match
                                    break
                                counter +=1
                            if len(v) == counter: #matching vertex found
                                t[vt+3] = [poolID4v, self.dsf.V[poolID4v].index(ev)]
                                self.log.info("  Vertex {} equals vertex {} in existing Pool with index {} .".format(v, ev, poolID4v))
                                matchfound = True
                                break
                            if counter >= 2:
                                self.log.warning("  Vertex {} is at same location as {} but different higher coordinate {}!!".format(v, ev, counter))
                        if not matchfound:       
                            self.dsf.V[poolID4v].append(v)
                            t[vt+3] = [poolID4v, len(self.dsf.V[poolID4v])-1] #change reference in area tria vertex to the last vertex in existing pool
                            self.log.info("  Vertex {} inserted in existing pool no. {}.".format(v, poolID4v))
                    else: #no existing pool fulfilling requirements was found, so new pool has to be created and added with v
                        if len(self.dsf.V) >= 65535: #reached maximum number of pools, no pool could be added any more
                            self.log.error("DSF File already has maximum number of point pools. Addtional pools required for change can not be added!!!")
                            return -1
                        self.dsf.V.append([v])
                        self.dsf.Scalings.append(deepcopy(self.dsf.Scalings[t[vt+3][0]])) #get scalings from the original tria vertex the new vertex is inside
                        if elevscal < 1: #for given submeter elevation pool-scaling has to be adapted
                            self.dsf.Scalings[-1][2][0] = 65535 * elevscal #new multiplier based on required scaling for elevation defined for new pool
                            self.dsf.Scalings[-1][2][1] = int(-500 + int((v[2]+500)/(65535 * elevscal)) * (65535 * elevscal)) #offset for elevation of v   ######## -500m is deepest vaule that can be defined with this routine ####
                        poolID4v = len(self.dsf.Scalings) - 1 #ID for the pool is the last one added
                        self.log.info("  New pool with index {} and scaling {} added to insert vertex {}.".format(poolID4v, self.dsf.Scalings[poolID4v], v))
                        newPools.append(poolID4v)
                        t[vt+3] = [len(self.dsf.V)-1, 0] #change reference in area tria vertex to the first vertex in new pool
                else:
                    self.log.info("Vertex {} unchanged. Index to pool {} reused.".format(t[vt], t[vt+3][0]))
        return 0
    
    def insertMeshArea(self):
        """
        Insert mesh area in dsf again from which it was extracted.
        Important: The shape of extraction must still be the same. Also no new vertices on borders are allowed.
                   Changes are only allowed inside the area.
        """
        patchTrias = {} #dictionary that stores for each patch the list of trias that are in the area in this patch
        for t in self.atrias:
            if t[6] in patchTrias:
                patchTrias[t[6]].append([t[3], t[4], t[5]])
            else:
                patchTrias[t[6]] = [ [t[3], t[4], t[5]] ]
        for p in patchTrias:
            self.log.info("For patch no. {} trias will be added: {}".format(p, patchTrias[p]))
            dsftrias = self.dsf.Patches[p].triangles()
            dsftrias.extend(patchTrias[p])
            self.dsf.Patches[p].trias2cmds(dsftrias)
                
