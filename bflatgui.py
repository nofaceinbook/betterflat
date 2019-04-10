# -*- coding: utf-8 -*-
"""
GUI for betterflat   version 0.0.1
by schmax 
"""
############### IMPROVE FILE HEADER !!!!! ############


############ ADD LOGGING MODULE !!!! #################
######### Handle several boundaries !!!!!!! ###########

from xplnedsf import *
import os
from tkinter import *
from tkinter.filedialog import askopenfilename


def displayHelp(win):
    newwin = Toplevel(win)
    Label(newwin, text="Display some information on the tool.....").grid(row=0, pady=10, padx=10)
    ############# IMPROVE HELPT TEXT !!!! ###########################
    

def readBoundaryFromAPT(filename, icao_id = ""):
#
# Reads coordinates from airport boundary in given apt.dat file and returns them in a list of [lon, lat] values
# In case apt.dat contains severla airports the right airport is selected by giving the correct ICAO Code
# Currently only works for single airports with single boundary defintion without holes!!!
# For Bezier nodes of the boundary only the node without Bezier definitions is considered!!!
#
    print ("Reading boundary from:", filename)
    if icao_id == "":
        Airport = True #if no icao code is given the first boundary in file is selected
    else:
        Airport = False #else first the correct icoa id has to be found befor Airport becomes true
    BoundarySection = False  #is set true when we are in boundary section in .apt file
    l = [] #list of boundary coordinates that will be returned
    if not os.path.isfile(filename):
        print("Error: Airport File", filename, "does not exist!")
        return [], "Error: Airport file does not exist!"
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
    if len(l) == 0:
        return l, "Warning: No valid boundary found in file!"
    else:
        return l, None


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
            

def vertices_of_boundary_intersecting_trias(dsf, poly):
#
# Identifies all mesh triangles of dsf file having an intersection with an edge of the boundary (poly with [lon, lat] coordinates)
# or lie inside boundary (or single triangle in which boundary completely lies)
# Returns a set with all vertices (tuple of pool_id and vertex number in pool) that belong to such intersecting triangles
#
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
    print(len(s), "vertices of mesh found that belong to mesh triangles intersecting or within boundary.")
    changedVertices = {} #set up dictonary with all coords of changed vertices to find additional vertices at such coords for trias out of bound
    for v in s:
        changedVertices[(round(dsf.V[v[0]][v[1]][0], 7), round(dsf.V[v[0]][v[1]][1], 7))] = True #coords round to range of cm
    miny, maxy, minx, maxx = dsf.BoundingRectangle(changedVertices)
    delta = 0.000001 #check also for vertices 0.1m outside boundary
    for p in dsf.PatchesInArea(miny - delta, maxy + delta, minx - delta, maxx + delta): 
        for t in p.triangles():
            for v in t:
                if (round(dsf.V[v[0]][v[1]][0], 7), round(dsf.V[v[0]][v[1]][1], 7)) in changedVertices:
                    s.add((v[0], v[1])) #add tuple of vertex at coords to make sure to get all vertices at the coords
    print("References to vertices to be changed:", s)
    print("Unique coords to be changed:", changedVertices)
    return s, changedVertices


def averageheight(dsf, vertices):
    count = 0
    sum = 0
    for v in vertices:
        count += 1
        sum += dsf.getElevation(dsf.V[v[0]][v[1]][0], dsf.V[v[0]][v[1]][1], dsf.V[v[0]][v[1]][2])
    if count == 0:
        return None
    else:
        return round( sum / count)


class bflatGUI:
    def __init__(self):
        
        self.dsf = XPLNEDSF()  #includes all information about the read dsf-file
        self.boundary = []     #includes vertices of boundary
        self.vertices = set([])#set of vertices that should get new height
        self.dsfreadfile = ""  #name of the file that is stored in self.dsf
        
        self.window = Tk()
        self.window.title("X-Plane betterflat v0.1")
        
        self.header = Label(self.window, text="Work with copies here to protect your original files!!")
        self.header.grid(row=0, column=0, columnspan=2)
        self.help_button = Button(self.window, text=' Help ', fg='red', command = lambda: displayHelp(self.window))
        self.help_button.grid(row=0, column=3, sticky=W, pady=4, padx=10)
        
        self.aptfile_label = Label(self.window, text="Airport File (apt.dat):")
        self.aptfile_label.grid(row=1, column=0,  sticky=W)
        self.apt_entry = Entry(self.window, width = 60)
        self.apt_entry.grid(row=1, column=1, columnspan=2, sticky=W)
        self.apt_select = Button(self.window, text='Select', command = lambda: self.select_file(self.apt_entry))
        self.apt_select.grid(row=1, column=3, sticky=W, pady=4, padx=10)
        
        self.aptid_label = Label(self.window, text="Airport ICAO Id:")
        self.aptid_label.grid(row=2, column=0,  sticky=W)
        self.aptid_entry = Entry(self.window, width = 6)
        self.aptid_entry.grid(row=2, column=1, sticky=W)
        self.aptid_info = Label(self.window, text="")
        self.aptid_info.grid(row=2, column=2, sticky=W)
            
        self.apt_read = Button(self.window, text='Read Boundary', command = lambda: self.read_apt(self.apt_entry.get(), self.aptid_entry.get()))
        self.apt_read.grid(row=3, column=0, sticky=E, pady=4)
        self.apt_status_label = Label(self.window, text="ICOA Id only required if file contains several airports!")
        self.apt_status_label.grid(row=3, column=1, columnspan=2, sticky=W)

        self.dsffile_label = Label(self.window, text="DSF File:")
        self.dsffile_label.grid(row=4, column=0,  sticky=W)
        self.dsf_entry = Entry(self.window, width = 60)
        self.dsf_entry.grid(row=4, column=1, columnspan=2, sticky=W)
        self.dsf_select = Button(self.window, text='Select', command = lambda: self.select_file(self.dsf_entry)).grid(row=4, column=3, sticky=W, pady=4, padx=10)

        self.dsf_read = Button(self.window, text='  Read DSF   ', state=DISABLED, command = lambda: self.read_dsf(self.dsf_entry.get()))
        self.dsf_read.grid(row=5, column=0, sticky=E, pady=4)
        self.dsf_status_label = Label(self.window, text="   this can take some minutes....")
        self.dsf_status_label.grid(row=5, column=1, sticky=W)
        
        self.result_text = Label(self.window, text="Result:")
        self.result_text.grid(row=6, column=0, sticky=E, pady=8)
        self.result_label = Label(self.window, text=" not yet")
        self.result_label.grid(row=6, column=1, sticky=W, pady=8)
        
        self.height_label = Label(self.window, text="New height (in m):")
        self.height_label.grid(row=7, column=0, pady=4, sticky=E)        
        self.height_entry = Entry(self.window, width = 7)
        self.height_entry.grid(row=7, column=1, sticky=W)
        
        self.newdsf_label = Label(self.window, text="Updated DSF file:")
        self.newdsf_label.grid(row=8, column=0,  sticky=E, pady=4)
        self.newdsf_entry = Entry(self.window, width = 60)
        self.newdsf_entry.grid(row=8, column=1, columnspan=2, sticky=W)
        self.newdsf_change = Button(self.window, text='Change', command = lambda: self.select_file(self.newdsf_entry))
        self.newdsf_change.grid(row=8, column=3, sticky=W, pady=4, padx=10)
        
        self.bakdsf_label = Label(self.window, text="Backup orig. DSF:")
        self.bakdsf_label.grid(row=9, column=0,  sticky=E, pady=4)
        self.bakdsf_entry = Entry(self.window, width = 60)
        self.bakdsf_entry.grid(row=9, column=1, columnspan=2, sticky=W)
        self.bakdsf_change = Button(self.window, text='Change', command = lambda: self.select_file(self.bakdsf_entry))
        self.bakdsf_change.grid(row=9, column=3, sticky=W, pady=4, padx=10)
        
        self.write_button = Button(self.window, text="Write DSF", state=DISABLED, command= lambda: self.write_dsf(self.height_entry.get(), self.newdsf_entry.get(), self.bakdsf_entry.get()))
        self.write_button.grid(row=10, column=0, sticky=E, pady=4)
        self.write_status_label = Label(self.window, text="Warning: existing files are overwritten!! Takes time...")
        self.write_status_label.grid(row=10, column=1, sticky=W, pady=4)
        
        mainloop()
        
    def read_apt(self, filename, icao):
        self.boundary, err = readBoundaryFromAPT(filename, icao)
        if err:
            self.apt_status_label.config(text=err)
        else:
            self.apt_status_label.config(text="Boundary with {} vertices read.".format(len(self.boundary)))
            self.dsf_read.config(state="normal")
            
    def read_dsf(self, filename):
        self.dsfreadfile = filename        
        err = self.dsf.read(filename) #returns value > 0 in case of errors
        if err == 1:
            self.dsf_status_label.config(text = "Error: File not found!")
        elif err == 2:
            self.dsf_status_label.config(text = "Error: File is not correct dsf format! Might be 7zipped.")
        else:
            self.dsf_status_label.config(text = "DSF with {} pools of vertices and {} patches read.".format(len(self.dsf.V), len(self.dsf.Patches)))
            self.vertices, coords = vertices_of_boundary_intersecting_trias(self.dsf, self.boundary)
            ah = averageheight(self.dsf, self.vertices)
            if ah == None:
                self.result_label.config(text = "No mesh intersection found. Check that files are correct.")
            else:
                self.write_button.config(state="normal")
                self.height_entry.delete(0, END)
                self.height_entry.insert(0, ah)
                self.result_label.config(text = "{} vertices at {} coords with average height {} to be adapted.".format(len(self.vertices), len(coords), ah))
                self.newdsf_entry.delete(0, END)
                self.newdsf_entry.insert(0, filename)
                self.bakdsf_entry.delete(0, END)
                self.bakdsf_entry.insert(0, filename + ".bak")
        
    
    def write_dsf(self, newheight, dsffile, bakfile):
        newheight = int(self.height_entry.get())
        print("Writing dsf file:", dsffile, "for height:", newheight, "and backup-file", bakfile)  ######### ADD CODE TO REALLY WRITE ##########
        try:
           os.replace(self.dsfreadfile, bakfile)
        except:
           print('Error:', self.dsfreadfile, 'can not be replaced!')
           self.write_status_label.config(text = "Error: Original file {} can not be replaced!".format(self.dsfreadfile))
        else:
            for v in self.vertices:
                print (self.dsf.V[v[0]][v[1]], "with height", self.dsf.getElevation(self.dsf.V[v[0]][v[1]][0], self.dsf.V[v[0]][v[1]][1], self.dsf.V[v[0]][v[1]][2]), "changed to: ", newheight )
                self.dsf.V[v[0]][v[1]][2] = newheight            
            self.write_status_label.config(text = "Wrtiting changes ...") ############### FLUSH missing !!!!!!!! ######
            self.dsf.write(dsffile)
            self.write_status_label.config(text = "Done.")
            print("Done.")
              
    def select_file(self, entry):
        filename = askopenfilename()
        entry.delete(0, END)
        entry.insert(0, filename)

   

main = bflatGUI()