# -*- coding: utf-8 -*-
#******************************************************************************
#
# muxp.py
#        
muxp_VERSION = "0.0.2 exp"
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
   
            
def showProgress(percentage):
    if percentage%10 == 0:
        log.info("dsf processing at {}%".format(percentage))



dsf_filename = "D:\\Programmierung\\muxp_0_1\\test.dsf"
log = defineLog('muxp', 'INFO', 'INFO') #no log on console for EXE version
log.info("Started muxp Version: {}".format(muxp_VERSION))

dsf = XPLNEDSF(LogName, showProgress)
log.info("Start reading dsf file: {}".format(dsf_filename))
dsf.read(dsf_filename)  # returns value > 0 in case of errors
log.info("dsf file with {} pools of vertices and {} patches read.".format(len(dsf.V), len(dsf.Patches)))


umesh = muxpArea(dsf, LogName)
umesh.extractMeshArea(-7.9790, -7.9614, -14.4083, -14.3783)
cps = umesh.cutEdges([-14.4026, -7.9650], [-14.3841, -7.9745])
kmlExport2(dsf, [[[-14.4083, -7.9790], [-14.4083, -7.9614], [-14.3783, -7.9614], [-14.3783, -7.9790], [-14.4083, -7.9790]]], umesh.atrias, "D:\\Programmierung\\muxp_0_1\\test_extract.dsf")
log.info("================")
umesh.createDSFVertices()
umesh.insertMeshArea()
dsf.write(dsf_filename+"new")
dsf.read(dsf_filename+"new")



#exarea = extractMeshArea(dsf, 50.4085, 50.4106, -125.1370, -125.1262)  ### CANADA
#cps = cutEdges(exarea, [-125.147, 50.412], [-125.119, 50.407]) ### CANADA

