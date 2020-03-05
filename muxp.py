# -*- coding: utf-8 -*-
#******************************************************************************
#
# muxp.py
#        
muxp_VERSION = "0.0.5 exp"
# ---------------------------------------------------------
# Python Tool: Mesh Updater X-Plane (muxp)
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
from muxp_math import *
from muxp_area import *
from xplnedsf2 import *
from muxpKMLexport2 import *
from os import path, replace
from shutil import copy2
from tkinter import *
from tkinter.filedialog import askopenfilename, askdirectory
from glob import glob
from sys import argv, exit  ### exit only for TESTING ###

def displayHelp(win):
    helpwin = Toplevel(win)
    Label(helpwin, anchor=W, justify=LEFT, text=
          "This program updates the mesh of X-Plane based on a configuration\n"
          "given in a text file (*.muxp). \n"
          "Via the config butten you set your X-Plane base folder and the folder\n"
          "where the updated dsf files are stored. Make sure that this folder\n"
          "has in the scenery.ini file higher priority as other dsf mesh files\n"
          "in order to make changes visible.\n\n"
          "MORE INFORMATION, source code and contact info are\n"
          "available at GitHub: https://github.com/nofaceinbook/muxp/\n\n"
          "Hope the tool helps you.    (c) 2020 by schmax (Max Schmidt)"   
          ).grid(row=0, pady=10, padx=10)



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
        return [], [], None, None, None, "Error: Airport file does not exist!"
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

class muxpGUI:
    def __init__(self, rf, muxpfile=None):
        self.runfile = rf
        self.xpfolder = ""
        self.muxpfolder = ""
        self.kmlExport = 0

        self.window = Tk()
        self.window.title("Mesh Updater X-Plane (muxp version: {})".format(muxp_VERSION))
 
        self.current_action = "none"  #e.g. set to read or write when operating on dsf file     

        self.dsf = XPLNEDSF(LogName, self.showProgress)  # includes all information about the read dsf-file


        self.header = Label(self.window, text="WARNING - This is still a test version.")
        self.header.grid(row=0, column=0, columnspan=2)
        self.config_button = Button(self.window, text=' Config ', fg='black', command=lambda: self.ConfigMenu())
        self.config_button.grid(row=0, column=2, sticky=W, pady=4, padx=10)        
        self.help_button = Button(self.window, text=' Help ', fg='red', command=lambda: displayHelp(self.window))
        self.help_button.grid(row=0, column=3, sticky=W, pady=4, padx=10)

        self.muxpfile_label = Label(self.window, text="muxp File (*.muxp):")
        self.muxpfile_label.grid(row=1, column=0, sticky=W)
        self.muxpfile_entry = Entry(self.window, width=60)
        self.muxpfile_entry.grid(row=1, column=1, columnspan=2, sticky=W)
        self.muxpfile_select = Button(self.window, text='Select', command=lambda: self.select_muxpfile(self.muxpfile_entry))
        self.muxpfile_select.grid(row=1, column=3, sticky=W, pady=4, padx=10)
        
        self.muxp_start = Button(self.window, text='  Start muxp   ', state=DISABLED, command=lambda: self.runMuxp(self.muxpfile_entry.get()))
        self.muxp_start.grid(row=6, column=0, sticky=E, pady=4)
        self.muxp_status_label = Label(self.window, text="   this can take some minutes....")
        self.muxp_status_label.grid(row=6, column=1, sticky=W)
        self.info_label = Label(self.window, text=" ")
        self.info_label.grid(row=7, column=1, sticky=W, pady=8)
        self.muxp_create = Button(self.window, text='  Create muxp   ', state=DISABLED, command=lambda: self.create_muxp())
        self.muxp_create.grid(row=9, column=0, sticky=E, pady=4)
        self.muxp_undo = Button(self.window, text='  Undo muxp   ', state=DISABLED, command=lambda: self.undo_muxp())
        self.muxp_undo.grid(row=9, column=2, sticky=E, pady=4)        
        
        self.getConfig(runfile)
        
        if muxpfile != None:
            self.select_muxpfile(self.muxpfile_entry, muxpfile)
            self.runMuxp(muxpfile) ## directly run muxpfile if it was given as argument
        
        log.info("GUI is set up.")
        mainloop()

    def getConfig(self, runfile):
        """
        Gets configuration from muxp.config from same directory as runfile.
        Opens ConfigMenu if file is not present and creates file.
        """
        filename = self.runfile[:self.runfile.rfind('.')]+'.config'
        log.info("Searching Config File: {}".format(filename))
        c, err = readMuxpFile(filename) ## config file self has same syntax as muxpfile
        if c==None:
            log.error("{}".format(err))
            return
        if c['muxpconfigversion'] != "1":
            log.error("Config file has wrong version ({} instead of 1)".format(c['muxpconfigversion']))
            return
        self.xpfolder = c["xpfolder"]
        self.muxpfolder = c["muxpfolder"]
        self.kmlExpoft = int(c["kmlExport"])
        ##### tbd: check that keys are all exist, correct type and open config menu if not ####################
        #if values_read != 4:
        #    log.info("Config not found, not complete or wrong version. Open Config Window!")
        #    self.ConfigMenu()
        
    def showProgress(self, percentage):
        if self.current_action == 'read':
            self.muxp_status_label.config(text = "read {} percent".format(percentage))
        elif self.current_action == 'write':
            self.muxp_status_label.config(text = "written {} percent".format(percentage))
        self.window.update()
        
    def select_muxpfile(self, entry, filename=None): #if file is set it is directly displayed
        if filename == None:
            filename = askopenfilename()
        entry.delete(0, END)
        entry.insert(0, filename)
        self.muxp_start.config(state="normal")
        self.muxpfile_select.config(state="disabled")
        
    def ConfigMenu(self):
        def select_file(entry): #if file is set it is directly displayed
            file = askdirectory()
            entry.delete(0, END)
            entry.insert(0, file)
        configwin = Toplevel(self.window)
        configwin.attributes("-topmost", True)
        toplabel = Label(configwin, anchor=W, justify=LEFT, text="Settings for muxp").grid(row=0, column=0, columnspan=2, pady=10, padx=10)
        xpfolder_label = Label(configwin, text="X-Plane base folder:")
        xpfolder_label.grid(row=1, column=0, pady=4, sticky=E)
        xpfolder_entry = Entry(configwin, width=60)
        xpfolder_entry.grid(row=1, column=1, sticky=W)
        xpfolder_entry.insert(0, self.xpfolder)
        xpfolder_select = Button(configwin, text='Select', command=lambda: select_file(xpfolder_entry))
        xpfolder_select.grid(row=1, column=3, sticky=W, pady=4, padx=10)
        muxpfolder_label = Label(configwin, text="Folder to updated dsf:")
        muxpfolder_label.grid(row=2, column=0, pady=4, sticky=E)
        muxpfolder_entry = Entry(configwin, width=60)
        muxpfolder_entry.grid(row=2, column=1, sticky=W)
        muxpfolder_entry.insert(0, self.muxpfolder)
        muxpfolder_select = Button(configwin, text='Select', command=lambda: select_file(muxpfolder_entry))
        muxpfolder_select.grid(row=2, column=3, sticky=W, pady=4, padx=10)
        kmlExportType = IntVar() # 1 if kml should be exported, 0 if not
        kmlExportType.set(self.kmlExport)
        kmlExportCB = Checkbutton(configwin, text="Export to kml ", variable=kmlExportType)
        kmlExportCB.grid(row=3, column=0, sticky=E, pady=4)
        save_button = Button(configwin, text='  Save  ', command=lambda: self.safeConfig(xpfolder_entry.get(), muxpfolder_entry.get(), kmlExportType.get()))
        save_button.grid(row=10, column=0, pady=4)
        
    def safeConfig(self, xf, mf, ke):
        self.xpfolder = xf
        self.muxpfolder = mf
        self.kmlExport = ke
        log.info("Saving config {}, {}, {}".format(xf, mf, ke))
        filename = self.runfile[:self.runfile.rfind('.')]+'.config'
        with open(filename, "w", encoding="utf8", errors="ignore") as f:
            f.write("muxpconfigversion:  1\n")
            f.write("xpfolder:  {}\n".format(self.xpfolder))
            f.write("muxpfolder:  {}\n".format(self.muxpfolder))
            f.write("kmlExport:  {}\n".format(self.kmlExport))


    def runMuxp(self, filename):
        """
        Updates the mesh based on the muxp file.
        """
        update, error = readMuxpFile(filename)
        if update == None:
            log.error("muxpfile does not exist!")
            ### tbd: show error in GUI ####
            return
        error = validate_muxp(update) ### tbd: check if all relevant values are in and transform all values
        if error != None:
            log.error(error) ### tbd: show error in GUI ###
        log.info("muxpfile id: {} version: {} for area:{} with {} commands read.".format(update["id"], update["version"], update["area"], len(update["commands"])))
        log.info("This muxp version only changes default XP mesh....")
        ### tbd: check for other dsf file as hd-mesh etc. in Custom scenery and use the dsf that is loaded to XPlane
        self.muxp_start.config(state="disabled")
        filename = self.xpfolder +  "/Global Scenery/X-Plane 11 Global Scenery/Earth nav data/" + get10grid(update["tile"][:3]) + get10grid(update["tile"][3:]) + "/" + update["tile"] +".dsf"
        log.info("Loading dsf file {}".format(filename))
        self.current_action = "read"
        self.dsf.read(filename)
        a = muxpArea(self.dsf, LogName)
        a.extractMeshArea(*update["area"])
        areabound = [(update["area"][2],update["area"][0]), (update["area"][2],update["area"][1]), (update["area"][3],update["area"][1]), (update["area"][3],update["area"][0]), (update["area"][2],update["area"][0]) ]
        commandbound1 = [update["commands"][0]["coordinates"][0], update["commands"][0]["coordinates"][1], update["commands"][0]["coordinates"][2], update["commands"][0]["coordinates"][3], update["commands"][0]["coordinates"][0] ]
        log.info("Command Parameters: {}".format(commandbound1))
        new_polys = a.cutPoly(commandbound1)
        new_polys.append(commandbound1)
        ############# WARNING: Commands use x, y and Area: y, x coordinates ---> should be unique !!!! #####
        kmlExport2(self.dsf, new_polys, a.atrias, "D:\\Programmierung\\muxp_0_1\\test_extract.dsf")
        ##### TBD: Perhaps start kmlExport directly with area instead of vertex coordinates ?? !! ###########
        ### go through commands, write dsf to muxp-folder

def get10grid(tile):
    """
    returns for a tile definition string like -122 or +47 the 10 rounded string -130 or +40 
    """
    if tile[0] == '-':
        s = 5
    else:
        s = 4.99
    if len(tile) == 3:
        return "{0:+03d}".format(round((int(tile)-s)/10)*10)
    else:
        return "{0:+04d}".format(round((int(tile)-s)/10)*10)

          
def readMuxpFile(filename):
    log.info("Reading muxp File: {}".format(filename))
    d = {} #dictionary returnig the read content
    d["commands"] = [] #in this dictionary entry are commands listed
    if not path.isfile(filename):
        return None, "Error: File not existent!"
    with open(filename, encoding="utf8", errors="ignore") as f:
        for line in f:
            line_indent = len(line) - len(line.lstrip(' ')) #only count spaces for indent
            line = line.lstrip() #remove new leading spaces, tabs etc
            if line.find('\n') >= 0:
                line = line[:line.find('\n')] #remove line feed
            if line.find('#') >= 0:
                line = line[:line.find('#')] #remove everything after # which is comment
            if line.find(':') >= 0:
                key = line[:line.find(':')]
                key = key.lstrip()
                value = line[line.find(':')+2:]
                value = value.lstrip()
                if line_indent == 0: #we are on toplevel
                    if len(value) == 0: #if no value given, this is a mesh command
                        log.info("Entering new command: {}".format(key))
                        new_command_dict = {}
                        new_command_dict["command"] = key
                        d["commands"].append(new_command_dict)
                    else:
                        log.info("Read key: {} assigned to: {} and indent: {}".format(key, value, line_indent))  
                        d[key] = value
                else: #we are inside a command
                    if len(value) == 0: #if no value given, this is a data list for last element of command
                        log.info("Entering new datalist: {}".format(key))
                        new_datalist = []
                        d["commands"][-1][key] = new_datalist
                    else:
                        log.info("Read key inside command: {} assigned to: {} and indent: {}".format(key, value, line_indent))  
                        d["commands"][-1][key] = value
            if line.find('-') == 0: #we have now data list element which should belong to data set in command ---> error checks should be done as weel
                new_datalist.append(line[2:])
                log.info("Data-element: {}".format(line[2:]))
    return d, None

def validate_muxp(d):
    """
    Validates values in read muxp dictionary d and turns them inside d to correc format
    or returs error in case walues don't match.
    """
    ### Validate muxp file version ###
    if float(d["version"]) < 0.01: ### IMPORTANT: This is current version for muxp file
        return "Error muxp file: version is too old."
    
    ### Extract and validate area defined ###    
    d["area"] = d["area"].split()
    for i in range(len(d["area"])):
        try:
            d["area"][i] = float(d["area"][i])
        except ValueError:
            return "Error muxp file: area argument {} is not a float.".format(i+1)
    if i != 3:
        return "Error muxp file: area has {} instead of 4 arguments.".format(i+1)
    ################ TBD: Check that 4 floats really define coordinates for an area #########################
    
    ### Extract and validate commands
    for i in range(len(d["commands"])):
        log.info("Command {}: {}".format(i, d["commands"][i]))
        ### Extract and validate cut_polygon ###
        if d["commands"][i]["command"] == "cut_polygon":
            for j in range(len(d["commands"][i]["coordinates"])):
                d["commands"][i]["coordinates"][j] = d["commands"][i]["coordinates"][j].split()
                d["commands"][i]["coordinates"][j][0] = float(d["commands"][i]["coordinates"][j][0])
                d["commands"][i]["coordinates"][j][1] = float(d["commands"][i]["coordinates"][j][1])
                ######## TBD: Check that reall two integers, and also extract elevation ############
    
    ####### TBD: Extract all other stuff ----> Better place this function and read muxp in an object??!!! ################
    
    return None
 
    
########### MAIN #############
log = defineLog('bflat', 'INFO', 'INFO') #no log on console for EXE version --> set first INFO to None
log.info("Started muxp Version: {}".format(muxp_VERSION))
runfile = argv[0]
if len(argv) > 1:
    muxpfile = argv[1]
else:
    muxpfile = None
main = muxpGUI(runfile, muxpfile)

