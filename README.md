# Update branch

This branchs is used to do my step by step update for the next version of betterflat.

## Warning: this code is in develpment and probably contains errors!

## Changes
Updates to xplnedsf in order to insert vertex in mesh, all changes are marked with NEW
In order to write changes functions for writing Scalings and CMDS to binary have been written

## Next Steps
Correct error when reading STRIP COMMANDS. Order of strip 1,2,3,4, 5 is not triangles 123 234 345 but 123 243 345!!!!!
Set trias and PatchBoundary not inside extractCMDS function but in an own step after extraction
Correct error when decoding 32-Bit Pools (uses still modulo of 16 bit pools) and include also function to pack 32-Bit scaling!
Continue functions for insterting vertex.
