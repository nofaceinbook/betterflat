

from math import sin, cos, sqrt, atan2, radians # for distance calculation etc.

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



 class xlpnedsfLefOvers:   
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
        #Now insert vertices on buttom in first tria of runway ### NEW 4 ###
        self._log_.info("INSERT VERTICES ON BUTTOM OF RUNWAY:")
        rwyPatch.trias.remove([rwyApv[0], rwyApv[1], rwyBpv[0]]) #remove first tria which is replaced below by new inserted trias for button vertices
        previousVonShoulder = rwyBpv[0]
        for k in sorted(dictShoulder[0].keys()):
            if k >= 0.01: #overjumps first corner of rwy on buttom, this is already in previousVonShoulder
                rwyPatch.trias.append([ previousVonShoulder, dictShoulder[0][k], rwyApv[1] ])
                previousVonShoulder = dictShoulder[0][k]
                #assumes that last k is key to second corner of runway button
        #Now insert vertices on top in last tria of runway ### NEW 4 ###
        self._log_.info("INSERT VERTICES ON TOP OF RUNWAY:")
        rwyPatch.trias.remove([rwyBpv[-2], rwyApv[-1], rwyBpv[-1]]) #remove last tria which is replaced below by new inserted trias for button vertices
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
    
