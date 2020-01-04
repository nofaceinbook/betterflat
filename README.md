# ATTENTION: IN DEVELOPMENT

## Requires 64 bit Python 3.7 (currently issue with Python 3.8)
## Branch to allow sloped runways

# Better Flattening for X-Plane Airports

betterflat (short **bflat**) flattens the area for a X-Plane airport better – meaning much more close to the airport boundary 
than the default X-Plane solution using the flatten flag. 
Only the triangles of the mesh that intersect with the airport boundary or that are lying inside the boundary 
are flattened to a height you can define. You also have the option to just flatten the mesh triangles under the runways.

## Advantages of bflat
* User friendly GUI
* No need to convert to dsf-text files. The tool directly runs on the dsf binary file.
* Easy to install and run (executable file for Windows or Python 3 scripts for other platforms only requiring standard libraries).
* No need to handle coordinates for mesh editing. All required information is directly retrieved from the apt.dat of the airport.
* Directly adapts the apt.dat file in case a flatten flag was defined for that airport.
* Depending on installed version the tool also operates on zipped dsf-files.
* Supports export of airport mesh area to kml-file. So you can check which mesh triangles will be flattened e.g. with Google Earth.

## Installation
For Windows you can just download and run the bflat.exe that includes all you need. This version will also read 7zipped dsf-fils directly. 
For other systems you download all the Pyhton 3 scripts: bfalt.py, xlpnedsf.py, bflatKMLexport.py. 
The scripts all run with standard libraries. Optional you can install pylzma (https://github.com/fancycode/pylzma)  
library in python. This allows you to directly read zipped dsf-files.  As this requires some more IT knowledge 
you can skip this and just manually unzip the 7zipped dsf files you want to adapt e.g. by using https://www.7-zip.org/

## Running the tool
Start the exe-file or run the bflat.py pyhton script. The GUI will open and steps should be self-explainable. 
**Important:** In order to avoid data loss make sure you always have copies of the original files you will modify.
The tool overwrites existing files without warning.

## Issues
Check out the FAQs included.
The tools creates a log file bflat.log. Open the file with a text editor. You might find additional information there about your problem.
Please report errors or give comments.
