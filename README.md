# Better Flattening for X-Plane Airports (version 0.4)

betterflat (short **bflat**) flattens the area for a X-Plane airport better – meaning much more close to the airport boundary 
than the default X-Plane solution using the flatten flag. 
Only the triangles of the mesh that intersect with the airport boundary or that are lying inside the boundary 
are flattened to a height you can define. You also have the option to just flatten the mesh triangles under the runways.
This version also offers cutting of triangles. When this option "cut" is ticked, then triangles are cut at runway or airport boundary. 
In the CONFIGURATION menu you can define the terrain that replaces the terrain within the boundary. Thus you can define airports how they will look like at the time LR has adapted the mesh for the airport based on your boundary for a later version of X-Plane.
In case SCIPY packe is installed blfat can also generate a runway profile. Refer to details below.

## Advantages of bflat
* User friendly GUI
* No need to convert to dsf-text files. The tool directly runs on the dsf binary file.
* Easy to install and run (executable file for Windows or Python 3 scripts for other platforms only requiring standard libraries).
* No need to handle coordinates for mesh editing. All required information is directly retrieved from the apt.dat of the airport.
* Directly adapts the apt.dat file in case a flatten flag was defined for that airport.
* New: Cut runway shape, runway profile or boundary in the mesh and change terrain if wanted.
* Depending on installed version the tool also operates on zipped dsf-files.
* Supports export of airport mesh area to kml-file. So you can check which mesh triangles will be flattened e.g. with Google Earth.

## Installation
For Windows you can just download and run the bflat.exe that includes all you need. This version will also read 7zipped dsf-fils directly. 
For other systems you need Python 64bit version installed (tested with version 3.6./3.7) and download all the Pyhton 3 scripts: bfalt.py, xlpnedsf.py, bflatKMLexport.py. 
The scripts all run with standard libraries. Optional you can install pylzma (https://github.com/fancycode/pylzma)  
library in python. This allows you to directly read zipped dsf-files.  As this requires some more IT knowledge 
you can skip this and just manually unzip the 7zipped dsf files you want to adapt e.g. by using https://www.7-zip.org/
If you want to generate runway profiles you need the "interpolate" function from SCIPY (www.SciPy.org). 

## Running the tool
Start the exe-file or run the bflat.py pyhton script. The GUI will open and steps should be self-explainable. 
**Important:** In order to avoid data loss make sure you always have copies of the original files you will modify.
The tool overwrites existing files without warning.

## Configuration Settings
In the main window you can click on "Config" to enter the configuration menu giving you the following options:
* Accuracy: If there are already other vertices or edges within this range they will be used instead of generating new ones. Thus the lower the value is the more accurate is your cut but the more triangles will be generated. The value 0.1 m is already very accurate.
* TerrainFile: When you use the replace function or the profile function then this terrain file is used for the triangles inside the cutted area. Note: Currently only terrain files without s/t coordinates are supported.
* Replace terrain: If this option is used the exisiting terrain will be replaced inside cutted area.
* Cut Runway Profile: Generates a profile for the runway. If now defintion is given below the profile is based on the raster inside the dsf file. 
* Profile Definition: Here you can define your own profile for the runway like -100@97 0@97 100@97 600@101 800@99 where the number before @ is the distance from runway start and the number after the height (both in meter). Values between are interpolated based on these numbers. If you define a profile, you need to give definition for every runway of the airport seperated with ";".

## Issues
Check out the FAQs included.
Cutting runway profile is not working very stable, especially if you have a HD mesh with many triangles. OFten it helps to modify parameters like higher accuracy or move boundaries/runay.
The tools creates a log file bflat.log. Open the file with a text editor. You might find additional information there about your problem.
Please report errors or give comments.
