# -*- coding: utf-8 -*-
#******************************************************************************
#
# muxpKMLexport2.py   for muxp
#        
muxpKMLexport2_VERSION = "0.0.3"
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

from os import fspath
from muxp_math import *

def kmlExport2(dsf, boundaries, extract, filename):
    
    #### Sort exported trias according to their patches, and get all vertices
    patchTrias = {} #dictionary that stores for each patch the list of trias that are in the area in this patch
    all_vertices = []
    for t in extract:
        all_vertices.extend(t[0:3]) #get coordinates of tria vertices
        if t[6] in patchTrias:
            patchTrias[t[6]].append(t)
        else:
            patchTrias[t[6]] = [ t ]
    
    #### Get the bounding rectangle over all boundaries and all tria vertices to be exported for area which is relevant for raster export
    for bounds in boundaries:
        all_vertices.extend(bounds)
    latS, latN, lonW, lonE = BoundingRectangle(all_vertices)

    #### Get index for raster pixel SW (yS, xW) for area to be exported
    xW = abs(lonW - int(dsf.Properties["sim/west"])) * (dsf.Raster[0].width - 1) # -1 from widht required, because pixels cover also boundaries of dsf lon/lat grid
    yS = abs(latS - int(dsf.Properties["sim/south"])) * (dsf.Raster[0].height - 1) # -1 from height required, because pixels cover also boundaries of dsf lon/lat grid
    if dsf.Raster[0].flags & 4: #when bit 4 is set, then the data is stored post-centric, meaning the center of the pixel lies on the dsf-boundaries, rounding should apply
        xW = round(xW, 0)
        yS = round(yS, 0)
    xW = int(xW) #for point-centric, the outer edges of the pixels lie on the boundary of dsf, and just cutting to int should be right
    yS = int(yS) 
    
    #### Get index for raster pixel NE (yN, xE) for area to be exported
    xE = abs(lonE - int(dsf.Properties["sim/west"])) * (dsf.Raster[0].width - 1) # -1 from widht required, because pixels cover also boundaries of dsf lon/lat grid
    yN = abs(latN - int(dsf.Properties["sim/south"])) * (dsf.Raster[0].height - 1) # -1 from height required, because pixels cover also boundaries of dsf lon/lat grid
    Rcentricity = "point-centric"
    if dsf.Raster[0].flags & 4: #when bit 4 is set, then the data is stored post-centric, meaning the center of the pixel lies on the dsf-boundaries, rounding should apply
        xE = round(xE, 0)
        yN = round(yN, 0)
        Rcentricity = "post-centric"
    xE = int(xE) #for point-centric, the outer edges of the pixels lie on the boundary of dsf, and just cutting to int should be right
    yN = int(yN) 
    
    #### Define relevant info for raster to be used later ####
    Rwidth = dsf.Raster[0].width
    xstep = 1 / (Rwidth - 1)  ##### perhaps only -1 when post-centric ---> also above !!! ########################################
    xbase = int(dsf.Properties["sim/west"])
    Rheight = dsf.Raster[0].height
    ystep = 1 / (Rheight -1)  ##### perhaps only -1 when post-centric ---> also above !!! ########################################
    ybase = int(dsf.Properties["sim/south"])



    
    filename = fspath(filename) #encode complete filepath as required by os
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
        
        ############## Style definitions for Raster Pixels ################
        f.write("<Style id=\"Raster0\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>80998066</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster1\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>80edefbe</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster2\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>80efff00</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster3\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>808fffff</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster4\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>8000beff</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster5\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>806596ff</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster6\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>806060ff</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster7\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>802c1ad3</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster8\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>8020206b</color></PolyStyle></Style>\n")
        f.write("<Style id=\"Raster9\"><LineStyle><color>ff000000</color><width>1</width></LineStyle><PolyStyle><color>80131340</color></PolyStyle></Style>\n")
        minelev = 10000
        maxelev = -500
        for x in range(xW, xE+1):
            for y in range (yS, yN+1):
                if dsf.Raster[0].data[x][y] < minelev:
                    minelev = dsf.Raster[0].data[x][y]
                if dsf.Raster[0].data[x][y] > maxelev:
                    maxelev = dsf.Raster[0].data[x][y]
        elevsteps = (maxelev - minelev) / 10  + 0.01 #add small value that the maxvalue is in last elevstep included
        
        ########### Show boundaries as Areas ################
        for boundary in boundaries:     
            f.write("    <Placemark><name>Selected Area</name><styleUrl>#Area</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n")
            for p in boundary:
                f.write("        {},{},0\n".format(p[0], p[1]))
            f.write("    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n")        

        ######### Export Raster #################
        f.write("<Folder><name>Raster {}x{} ({})from {}m to {}m </name>\n".format(Rwidth, Rheight, Rcentricity, minelev, maxelev))
        elevfolders = [[], [], [], [], [], [], [], [], [], []]
        if Rcentricity == "post-centric": #if post-centricity we have to move dem pixel half width/hight to left/down in order to get pixel center on border of dsf tile
            cx = 0.5 * xstep 
            cy = 0.5 * ystep
        else:
            cx = 0
            cy = 0
        for x in range(xW, xE+1):
            for y in range (yS, yN+1):
                folder = int((dsf.Raster[0].data[x][y] - minelev)/elevsteps)
                elevfolders[folder].append("    <Placemark><name>Pixel {}:{} at {} m</name><styleUrl>#Raster{}</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n".format(x, y, dsf.Raster[0].data[x][y], folder ))
                elevfolders[folder].append("        {},{},{}\n".format(xbase + x*xstep - cx, ybase + y*ystep - cy, dsf.Raster[0].data[x][y] ))
                elevfolders[folder].append("        {},{},{}\n".format(xbase + x*xstep - cx, ybase + (y+1)*ystep - cy, dsf.Raster[0].data[x][y] ))
                elevfolders[folder].append("        {},{},{}\n".format(xbase + (x+1)*xstep - cx, ybase + (y+1)*ystep - cy, dsf.Raster[0].data[x][y] ))
                elevfolders[folder].append("        {},{},{}\n".format(xbase + (x+1)*xstep - cx, ybase + y*ystep - cy, dsf.Raster[0].data[x][y] ))
                elevfolders[folder].append("        {},{},{}\n".format(xbase + x*xstep - cx, ybase + y*ystep - cy, dsf.Raster[0].data[x][y] ))
                elevfolders[folder].append("    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n")
        for folder in range(10):
            f.write("<Folder><name>Raster from {}m to {}m</name>\n".format(int(minelev + folder*elevsteps), int(minelev + (folder+1)*elevsteps)))
            for line in elevfolders[folder]:
                f.write(line)
            f.write("</Folder>\n")
        f.write("</Folder>\n")

        ########### Export Trias per Patch ##############
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
