# Update branch

This branchs is used to do my step by step update for the next version of betterflat.

## Warning: this code is in develpment and probably contains errors!

## Changes
* Updates to xplnedsf in order to insert vertex in mesh, all changes are marked with NEW
* In order to write changes functions for writing Scalings and CMDS to binary have been written
* Correct error when reading STRIP COMMANDS. Order of strip 1,2,3,4, 5 is not triangles 123 234 345 but 123 243 345!!!!!
* Trias are stored explicetley with patch, not only the building cmds
* Function for inserting vertex in simple tria (only 5 planes for vertex)

## Next Steps
* Correct error when decoding 32-Bit Pools (uses still modulo of 16 bit pools) and include also function to pack 32-Bit scaling!
* Continue functions for insterting vertex: general for all type of triangles
* Think of reducing needed storage amount. E.g. trias are stored as binary atom, in full CMDS, in comds with patch and as trias array with patch. 
