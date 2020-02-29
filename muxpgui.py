# -*- coding: utf-8 -*-
#******************************************************************************
#
# muxp.py
#        
muxp_VERSION = "0.0.4 exp"
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
        
        self.muxp_start = Button(self.window, text='  Start muxp   ', state=DISABLED, command=lambda: self.read_dsf(self.dsf_entry.get()))
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
            #### directly starte hier muxp-creation #############################
        
        log.info("GUI is set up.")
        mainloop()

    def getConfig(self, runfile):
        """
        Gets configuration from muxp.config from same directory as runfile.
        Opens ConfigMenu if file is not present and creates file.
        """
        filename = self.runfile[:self.runfile.rfind('.')]+'.config'
        log.info("Searching Config File: {}".format(filename))
        values_read = 0
        version = "no-def"
        if path.isfile(filename):
            log.info("Read config from file....")
            with open(filename, encoding="utf8", errors="ignore") as f:
                for line in f:
                    v = line.split()
                    if len(v) >= 2: #only consider value pairs
                        if v[0] == "muxpconfigversion:" and v[1] == "1":
                            version = v[1]
                            values_read += 1
                        if v[0] == "xpfolder:":
                            self.xpfolder = " ".join(v[1:]) #folder could include spaces so join it together ######## ISSUE with multiple spaces ####
                            values_read += 1
                        if v[0] == "muxpfolder:":
                            self.muxpfolder = " ".join(v[1:]) #folder could include spaces so join it together ######## ISSUE with multiple spaces ####
                            values_read += 1
                        if v[0] == "kmlExport:":
                            self.kmlExport = int(v[1])
                            values_read += 1
                log.info("Configfile version {} defined: {}, {}, {}".format(version, self.xpfolder, self.muxpfolder, self.kmlExport))
        if values_read != 4:
            log.info("Config not found, not complete or wrong version. Open Config Window!")
            self.ConfigMenu()
        
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

        
########### MAIN #############
log = defineLog('bflat', 'INFO', 'INFO') #no log on console for EXE version --> set first INFO to None
log.info("Started muxp Version: {}".format(muxp_VERSION))
runfile = argv[0]
if len(argv) > 1:
    muxpfile = argv[1]
else:
    muxpfile = None
main = muxpGUI(runfile, muxpfile)

