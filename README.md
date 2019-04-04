# betterflat for X-Plane
Python 3 script adapts an X-Plane dsf-file to flatten height at an airport.
Needs the included file xplnedsf.py for reading/writing dsf-files.
Note: dsf-files might be 7-zipped. In that case you need to unzip it first with some 7zip tool.

CAUTION: Work in progress.... Use with are at own risk.
         Run scripts always on copies of your original dsf / apt files.
         
Usage: betterflat <dsf-filename> <apt.dat including airport> <opt: 4 letter uppercase ICAO Code> <opt: new heigth in m>

Area to flatten is selected from airport boundary in 'apt.dat'. If apt.dat includes several airports, select airport by ICAO-Code.
Works for the moment only for single boundary without holes.
The new heigt is either given as optional argument in meter or it is derived from average height of vertices at boundary.
