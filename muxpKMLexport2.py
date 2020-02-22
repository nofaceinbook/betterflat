# -*- coding: utf-8 -*-
#******************************************************************************
#
# muxpKMLexport2.py   for muxp
#        
muxpKMLexport2_VERSION = "0.0.2"
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

import os

def kmlExport2(dsf, boundaries, extract, filename):
 
    patchTrias = {} #dictionary that stores for each patch the list of trias that are in the area in this patch
    for t in extract:
        if t[6] in patchTrias:
            patchTrias[t[6]].append(t)
        else:
            patchTrias[t[6]] = [ t ]

    
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
        for boundary in boundaries:     
            f.write("    <Placemark><name>Selected Area</name><styleUrl>#Area</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n")
            for p in boundary:
                f.write("        {},{},0\n".format(p[0], p[1]))
            f.write("    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n")        

        for p_num in patchTrias:
            p = dsf.Patches[p_num]
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
                
            f.write("<Folder><name>Patch {} ({}): {}</name>\n".format(p_num, flag, terrain))
            tcount = 0

            for t in patchTrias[p_num]:
                upcoords = "{} - ".format(t[3:6])  ########### Following NEW /EXPERIMENTAL to get upper coordinates shown in Google Earth
                upi = 5
                while upi < len(t[0]):
                    upxa = ", {0:.6f}".format(t[0][upi])
                    upcoords += upxa
                    upi += 1
                upi = 5
                upcoords += "/"
                while upi < len(t[1]):
                    upxa = ", {0:.6f}".format(t[1][upi])
                    upcoords += upxa
                    upi += 1
                upi = 5
                upcoords += "/"
                while upi < len(t[2]):
                    upxa = ", {0:.6f}".format(t[2][upi])
                    upcoords += upxa
                    upi += 1                    
                f.write("    <Placemark><name>T{} {}</name><styleUrl>#{}</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n".format(tcount, upcoords, style))  ##upcords NEW/Experimental
                h = [] #stores heigth of vertices in triangles
                h.append(int(dsf.getVertexElevation(t[0][0], t[0][1], t[0][2])))  #3rd Value is height from Vertex to be consideredn in case differnet from -32xxx
                h.append(int(dsf.getVertexElevation(t[1][0], t[1][1], t[1][2])))
                h.append(int(dsf.getVertexElevation(t[2][0], t[2][1], t[2][2])))
                f.write("        {0},{1},{2} {3},{4},{5} {6},{7},{8} {0},{1},{2}\n".format(t[0][0], t[0][1], h[0], t[1][0], t[1][1], h[1], t[2][0], t[2][1], h[2]))
                f.write("    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n")
                tcount += 1
            f.write("</Folder>\n")

        f.write("</Document></kml>\n")
        
