from xplnedsf import *
import os
import sys

def printInfo():
    print ("------------------------------------------------------------------------------")
    print ("             X-PLANE betterflat   version 0.12 (test)  by schmax")
    print ("This script adapts an X-Plane dsf-file to flatten height at an airport.")
    print ("Usage: betterflat <dsf-filename> <apt.dat including airport>")
    print ("               <opt: 4 letter uppercase ICAO Code> <opt: new heigth in m>")
    print ("Area to flatten is selected from airport boundary in apt.dat")
    print ("If apt.dat includes several airports, select airport by ICAO-Code.")
    print ("Works for the moment only for single boundary without holes.")
    print ("The new heigt is either given as optional argument in meter or")
    print ("it is derived from average height of vertices at boundary.")
    print ("CAUTION: This is still a test-version. No warrenties/guaranties!")
    print ("         Existing .org files are overwritten!!!!!")
    print ("-------------------------------------------------------------------------")


def _linsolve_(a1, b1, c1, a2, b2, c2):
    divisor = (a1 * b2) - (a2 * b1)
    if divisor == 0:
        return (-99999, -99999)  #actually None but with negative values calling intersection function returns None
    return (round(((c1 * b2) - (c2 * b1)) / divisor, 8), round(((a1 * c2) - (a2 * c1)) / divisor, 8)) #ROUNDING TO ALWAYS GET SAME CONCLUSION
 
def intersection(p1, p2, p3, p4): #checks if segment from p1 to p2 intersects segement from p3 to p4 
    s0, t0 = _linsolve_(p2[0] - p1[0], p3[0] - p4[0], p3[0] - p1[0], p2[1] - p1[1], p3[1] - p4[1], p3[1] - p1[1])
    if s0 >= 0 and s0 <= 1 and t0 >= 0 and t0 <= 1:
        return( round((p1[0] + s0 * (p2[0] - p1[0])), 8), round(p1[1] + s0 * (p1[1] - p1[1]), 8) )  ### returns the cutting point as tuple; ROUNDING TO ALWAYS GET SAME POINT
    else:
        return(None)   

def PointInPoly(p, poly): #test wether a point p with [lat, lon] coordinates lies in polygon (list of [lat, lon] pairs
# counts number of intersections from point outside poly to p on same y-coordinate, if it is odd the point lies in poly
# to avoid intersection at vertex of poly on same y-coordinate, such points are shifte about 1mm above for testing intersection
    count = 0
    for i in range(len(poly) - 1): #for all segments in poly
        epsilon0, epsilon1 = (0, 0) #added to p's y coordinate in case p is on same y-coordinate than according vertex of segment
        if poly[i][1] < p[1] and poly[i+1][1] < p[1]: #if vertices of segment below y-coordinate of p, no intersection
            continue
        if poly[i][1] > p[1] and poly[i+1][1] > p[1]: #if vertices of segment above y-coordinate of p, no intersection
            continue
        if poly[i][1] == p[1]:
            epsilon0 = 0.00000001
        if poly[i+1][1] == p[1]:
            epsilon1 = 0.00000001
        x = intersection([poly[i][0], poly[i][1] + epsilon0], [poly[i+1][0], poly[i+1][1] + epsilon1], [181, p[1]], p )
        if x:
            count += 1
    if count % 2:
        return True #odd number of intersections, so p in poly
    else:
        return False #even number of intersection, so p outside poly
            
            

def readBoundaryFromAPT(filename, icao_id = ""):
#
# Reads coordinates from airport boundary in given apt.dat file and returns them in a list of [lon, lat] values
# In case apt.dat contains severla airports the right airport is selected by giving the correct ICAO Code
# Currently only works for single airports with single boundary defintion without holes!!!
# For Bezier nodes of the boundary only the node without Bezier definitions is considered!!!
#
    if icao_id == "":
        Airport = True #if no icao code is given the first boundary in file is selected
    else:
        Airport = False #else first the correct icoa id has to be found befor Airport becomes true
    BoundarySection = False  #is set true when we are in boundary section in .apt file
    l = [] #list of boundary coordinates that will be returned
    if not os.path.isfile(filename):
        print("Error: Airport File", filename, "does not exist!")
        sys.exit(1)
    with open(filename, encoding="utf8") as f:
        for line in f:
            v = line.split()
            if len(v) == 0: #don't consider empty lines
                continue
            if len(v) > 2:
                if v[0] == '1302' and v[1] == 'icao_code' and v[2] == icao_id:
                    Airport = True
            if Airport == True and v[0] == '130':
                BoundarySection = True
            elif BoundarySection:
                if v[0] == '111' or v[0] == '112':
                    l.append([float(v[2]), float(v[1])]) #Note: Bezier definitins are not considered, just the base point
                elif v[0] == '113' or v[0] == '114':
                    l.append([float(v[2]), float(v[1])]) #Note: Bezier definitins are not considered, just the base point
                    BoundarySection = False
                    Airport = False #stop after read first boundary!!
    if l != []:
        l.append(l[0]) #form closed loop by adding again first vertex
    return l
                    

### MAIN ###
newheight = None  #new height that should be applied not defined yet
icao_id = ""      #icao code in case of searching apt.dat including several airports

try:
    dsffile = sys.argv[1]
    aptfile = sys.argv[2]
    if len(sys.argv) == 4:
        if sys.argv[3].isdigit():
            newheight = int(sys.argv[3])
        else:
            icao_id = sys.argv[3]
    elif len(sys.argv) == 5:
        icao_id = sys.argv[3]
        newheight = int(sys.argv[4])
except:
    printInfo()
    sys.exit()
    
dsffile = os.fspath(dsffile)
aptfile = os.fspath(aptfile)


print ("Reading", aptfile, "to extract boundary of airport", icao_id, flush = True)
poly = readBoundaryFromAPT(aptfile, icao_id) 
if len(poly):
    print ("Extracted boundary with", len(poly), "vertices.", flush = True)
else:
    print ("Error: No boundary found in ", aptfile)
    sys.exit(2)

dsf = XPLNEDSF()
if dsf.read(dsffile): #returns value > 0 in case of errors
    sys.exit(3)

print("Reading of dsf-file completed. Checking now for mesh triangles intersecting airport boundary.......", flush = True)
miny, maxy, minx, maxx = dsf.BoundingRectangle(poly)
s = set([])
for p in dsf.PatchesInArea(miny, maxy, minx, maxx): 
    for t in p.triangles():
        TriaV = dsf.TriaVertices(t)
        TriaV.append(TriaV[0]) #append first vertex to get a closed connection
        for i in range(3): #3 sides of triangle
            for j in range(len(poly)-1): #sides of poly
                if intersection(TriaV[i], TriaV[i+1], poly[j], poly[j+1]): #current triangle intersects with an poly line
                    s.add((t[0][0], t[0][1])) #add all 3 vertices of tria to set as tuples
                    s.add((t[1][0], t[1][1]))
                    s.add((t[2][0], t[2][1]))
        if PointInPoly(TriaV[0], poly): #Check also that not complete Tria lies in poly by checking for first vertex of Tria
            s.add((t[0][0], t[0][1])) #add all 3 vertices of tria to set as tuples
            s.add((t[1][0], t[1][1]))
            s.add((t[2][0], t[2][1]))        
        if PointInPoly(poly[0], TriaV): #Check also that not complete poly lies in current tria by checking for first vertex of poly
            s.add((t[0][0], t[0][1])) #add all 3 vertices of tria to set as tuples
            s.add((t[1][0], t[1][1]))
            s.add((t[2][0], t[2][1]))
print(len(s), "vertices of mesh found that belong to mesh triangles intersecting boundary.")

count = 0
sum = 0
for v in s:
    count += 1
    sum += dsf.getElevation(dsf.V[v[0]][v[1]][0], dsf.V[v[0]][v[1]][1])
averageheight = round(sum / count)
print ("Average height of boundary vertexes:", averageheight)

if newheight == None:
    newheight = averageheight
    print("Take this average height", newheight,"(in m) as new height for according vertices.")
else:
    print("Take given height", newheight,"(in m) as new height for according vertices.")
    
for v in s:
    print (dsf.V[v[0]][v[1]], "with height", dsf.getElevation(dsf.V[v[0]][v[1]][0], dsf.V[v[0]][v[1]][1]), "changed to: ", newheight )
    dsf.V[v[0]][v[1]][2] = newheight

print("Renaming original dsf-file to:", dsffile)
print("CAUTION: If you run this tool again this original file might be overwritten!")
try:
   os.replace(dsffile, dsffile + ".org")
except:
   print('Error:', dsffile, 'can not be replaced!')
   sys.exit(1)
print("Applied changes will now be written to:", dsffile)
dsf.write(dsffile)
