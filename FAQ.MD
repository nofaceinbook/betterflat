## How does bflat work?
The tool extracts the boundary polygon defined for that airport or the runway defintions from the airport defintion file apt.dat. 
Now it searches all mesh triangles of the dsf-file that intersect the boundary or that lie inside the boundary. 
All vertices of these triangles (and other vertices that have the same coordinates) will be set to the same new height. 
This will result in a flat airport.

## Why are the written dsf files so much bigger?
The tool does not write zipped archive. Just use a 7zip tool to zip the dsf file in case you like to save disk space.

## Why does the tool does not find a mesh intersection?
Usually you have the wrong dsf file used. Make sure that the dsf-file is really for the 1 by 1 grid where the airport is on. 
Also make sure you used a dsf file including the mesh. The dsf file that is delivered together with an airport usually
includes the placement of buildings but no mesh.

## What to do if my airport goes over several dsf files?
In case your airport is not within a one by one coordinate grid you have to flatten all the according dsf files with bflat.

## How to identify an airport without ICAO code ?
Identify the airport with the code that is in front of the airport name inside the apt.dat file. Such line starts  with 1. 
In the following example use "CA446":
´´´
1 321 0 0 CA446 Stuart Island Airstrip 
´´´

## Better choose airport boundary or runways for flattening?
Choosing runways will lead to a much more strict flattened area. However you than have issues with taxi ways and buildings.
For small airports or airports without boundary, runway might be a better choice. But usually you should use the airport
boundary to be on the safe side.

## How are Bezier type vertices handled?
Bezier curves are not part of calculation. Just the first coordinates of Bezier types are used. So the strict lines 
from point to point of the boundary will be considered but not the curve.

## What to do if my airport has holes defined for the boundaries?
Holes will not be considered as holes but as part of the airport. But having holes is very rarely used.

## Can I flat areas outside an airport?
Yes. Just create your own boundary for the area you want to have flattened and use this filename for reading the boundary. 
Refer to example below.

## How should a boundary file look like?
Just copy a boundary from an apt.dat file or create your own file that looks like:
´´´
1 <airport height in feet> 0 0 <icoa-id> <airport name>
130 Airport Boundary 1
111  50.41021599 -125.13696606
111  50.40996632 -125.13692890
111      ...          ...
113  50.41052015 -125.13682289
´´´
The first line defines the airport/area. The line starting with 130 defines the beginning of th boundary. 
The last line needs to be 113 to end the boundary. This must not be same coordinates as first ones. 
The polygon is automatically closed from the last coordinates to the first. Instead of 111 you could have 112 
and instead of 113 you can have 114 for Bezier type vertices. However the tool will not consider Bezier style. 
It will just take the first coordinates given.
