# betterflat for X-Plane
Python 3 script that flattens a X-Plane dsf-file at the location of an airport.
Needs the included file xplnedsf.py for reading/writing dsf-files.

**Note:** dsf-files might be 7-zipped. In that case you need to unzip it first with some 7zip tool.

**CAUTION:** Work in progress.... Use with it with care.
         Run scripts always on copies of your original dsf / apt files.
         
**Usage:** betterflat <dsf-filename> <apt.dat including airport> <opt: 4 letter uppercase ICAO Code> <opt: new heigth in m>

Area to flatten is selected from airport boundary in 'apt.dat'. If apt.dat includes several airports, select airport by ICAO-Code.
Works for the moment only for single boundary without holes.
The new heigt is either given as optional argument in meter or it is derived from average height of all vertices found inside or near boundary that will be set to the same height.

**Installation:** Just download the two pyhton files and run it (python 3 needs to be installed on your computer)

        **As I'm still improving the tool I would be happy for your comments**


## Why should I use this tool
* If you created a new airport where the underground is not flat
* If you use a new mesh that is not prepared for that airport
* Much better than using flattenig option of X-Plane as only mesh triangles of airport area (depending on it's boundary) are changed and not hughe area as it is performed by X-Plane flattening
* No need to convert to dsf-text files - this script runs directly on the binary dsf-file
* Only standard libraries used meaning: download and run it directly


## FAQ

### How it works?
The scripts extracts the boundary polygon defined for that airport. Now it searches all mesh triangles of the dsf-file that intersect the boundary or that lie inside the boundary. All vertices of these triangles (and other vertices that have the same coordinates) will be set to the same new height. Which results in a flat airport.

### What to do if my airport has several boundaries?
You could run the tool several times always using a new boundary. Just copy the boundaries to seperate files and run the tool for them (not giving an icao-code).

### Can I also flatten other areas of the dsf?
Yes. Just create your own boundary for the area you want to have flattened and use this filename. Refer to example below.

## How should a boundary file look like?
Just copy a boundary from an apt.dat file or create your own file that looks like:
```
130 Airport Boundary 1
111  50.41021599 -125.13696606
111  50.40996632 -125.13692890
111      ...          ...
113  50.41052015 -125.13682289
```
The text behind 130 first line can be choosen free.
The last line needs to be 113. This must not be same coordinates as first ones. The polygon is automatically closed from the last coordinates to the first.
Instead of 111 you could have 112 and instead of 113 you can have 114 for Bezier type vertices. However the tool will not consider Bezier style. It will just take the first coordinates given.
If you use the tool for such a created file, please don't assign an icao-code.



