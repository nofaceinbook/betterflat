# -*- coding: utf-8 -*-
#******************************************************************************
#
# bflatKMLexport.py
#        
bflatKMLexport_VERSION = "0.4.0"
# ---------------------------------------------------------
# Python module for exporting mesh area to be flattened to KML-file.
# This module is called by bflat.py (Tool for flattening X-Plane Mesh)
#
# For more details refert to GitHub: https://github.com/nofaceinbook/betterflat

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

from xplnedsf import *
import os

def kmlExport(dsf, boundaries, vertices, filename):
    
    filename = os.fspath(filename) #encode complete filepath as required by os
    with open(filename + ".kml", "w") as f:
        f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        f.write("<kml xmlns=\"http://www.opengis.net/kml/2.2\" >\n")
        f.write("<Document>\n")
        
        ############## Style definitions for Polygons in kml ###############
        f.write("<Style id=\"Water\"><LineStyle><width>1</width></LineStyle><PolyStyle><color>40ff0000</color></PolyStyle></Style>\n")
        f.write("<Style id=\"grass\"><LineStyle><width>1</width></LineStyle><PolyStyle><color>407fffaa</color></PolyStyle></Style>\n")
        f.write("<Style id=\"FLAT\"><LineStyle><color>ffff00aa</color><width>3</width></LineStyle><PolyStyle><color>40f7ffff</color></PolyStyle></Style>\n")
        f.write("<Style id=\"ELSE\"><LineStyle><width>1</width></LineStyle><PolyStyle><color>4000aaaa</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Area\"><LineStyle><color>ff0000ff</color><width>4</width></LineStyle><PolyStyle><fill>0</fill></PolyStyle></Style>\n")
        
        ########### Show boundaries as Areas ################
        all_bounds = []
        for boundary in boundaries:     
            f.write("    <Placemark><name>Selected Area</name><styleUrl>#Area</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n")
            for p in boundary:
                f.write("        {},{},0\n".format(p[0], p[1]))
            f.write("    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n")
            all_bounds.extend(boundary)
        
        miny, maxy, minx, maxx = dsf.BoundingRectangle(all_bounds)
        pcount = 0
        for p in dsf.PatchesInArea(miny, maxy, minx, maxx):
            if p.flag == 1:
                flag = "PYS"
            else:
                flag = "OVL"
            terrain = dsf.DefTerrains[p.defIndex]
            if "Water" in terrain:
                style = "Water"
            elif "grass" in terrain:
                style = "grass"
            else:
                style = "ELSE"
                
            f.write("<Folder><name>Patch {} ({}): {}</name>\n".format(pcount, flag, terrain))
            tcount = 0
            saved_style = style #saved as style might be temporarely changed
            for t in p.triangles():
                if (t[0][0], t[0][1]) in vertices and (t[1][0], t[1][1]) in vertices and (t[2][0], t[2][1]) in vertices:
                    style = "FLAT" #if all vertices of current triangle are in vertices where height is adapted, then triangle is flattened
                else:
                    style = saved_style #set back to original style
                upcoords = "{} - ".format(t)  ########### Following NEW /EXPERIMENTAL to get upper coordinates shown in Google Earth
                upi = 5
                while upi < len(dsf.V[t[0][0]][t[0][1]]):
                    upxa = ", {0:.6f}".format(dsf.V[t[0][0]][t[0][1]][upi])
                    upcoords += upxa
                    upi += 1
                upi = 5
                upcoords += "/"
                while upi < len(dsf.V[t[1][0]][t[1][1]]):
                    upxa = ", {0:.6f}".format(dsf.V[t[1][0]][t[1][1]][upi])
                    upcoords += upxa
                    upi += 1
                upi = 5
                upcoords += "/"
                while upi < len(dsf.V[t[2][0]][t[2][1]]):
                    upxa = ", {0:.6f}".format(dsf.V[t[2][0]][t[2][1]][upi])
                    upcoords += upxa
                    upi += 1                    
                f.write("    <Placemark><name>T{} {}</name><styleUrl>#{}</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n".format(tcount, upcoords, style))  ##upcords NEW/Experimental
                v = dsf.TriaVertices(t)
                h = [] #stores heigth of vertices in triangles
                h.append(int(dsf.getElevation(v[0][0], v[0][1], dsf.V[t[0][0]][t[0][1]][2])))  #3rd Value is height from Vertex to be consideredn in case differnet from -32xxx
                h.append(int(dsf.getElevation(v[1][0], v[1][1], dsf.V[t[1][0]][t[1][1]][2])))
                h.append(int(dsf.getElevation(v[2][0], v[2][1], dsf.V[t[2][0]][t[2][1]][2])))
                f.write("        {0},{1},{2} {3},{4},{5} {6},{7},{8} {0},{1},{2}\n".format(v[0][0], v[0][1], h[0], v[1][0], v[1][1], h[1], v[2][0], v[2][1], h[2]))
                f.write("    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n")
                tcount += 1
            f.write("</Folder>\n")
            pcount += 1
        
        f.write("</Document></kml>\n")
        
