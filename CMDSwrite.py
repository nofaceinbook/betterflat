from xplnedsf import *
import copy

dsf = XPLNEDSF()
dsf.read("filname.dsf")

with open("CMDS_org.txt", "w+") as fc:
    for d in dsf.DefObjects:
        print(d, ':', dsf.DefObjects[d], ':', dsf.Objects[d])
    for d in dsf.DefPolygons:
        print(d, ':', dsf.DefPolygons[d], ':', dsf.Polygons[d])
    for d in dsf.DefNetworks:
        print(d, ':', dsf.DefNetworks[d])        
    for c in dsf.CMDS:
        fc.write(", ".join(str(n) for n in c))
        fc.write("\n")
        
with open("CMDS_new.txt", "w+") as fcn:
    flag_physical = None #1 if physical, 2 if overlay
    nearLOD = None  
    farLOD = None
    poolIndex = None
    defIndex = None
    subroadtype = None
    junctionoffset = None
    enccmds = struct.pack('<BB', 3, 0) ###### TO BE CHANGED: JUST DEFINE EMPTY STRING HERE #####
    for d in dsf.DefObjects:
        fcn.write("3, %s\n" % str(d)) ##### TO BE CHANGED WHEN index > 255 or > 16 Bit Value   --> CMD ID 4 or 5!!!
        enccmds += struct.pack('<BB', 3, d)
        for c in dsf.Objects[d]:
            if c[0] != poolIndex:
                fcn.write("1, %s\n" % str(c[0]))
                enccmds += struct.pack('<BH', 1, c[0])
                poolIndex = c[0]
            fcn.write(", ".join(str(n) for n in c[1:]))
            fcn.write("\n")
            if c[1] == 7:
                enccmds += struct.pack('<BH', 7, c[2])
            elif c[1] == 8:
                enccmds += struct.pack('<BHH', 8, c[2], c[3])
            else:
                print("ERROR!!! Unknown COMMAND ID!!!!!!!!!!!!!!!!")
                
    for d in dsf.DefPolygons:
        fcn.write("3, %s\n" % str(d)) ##### TO BE CHANGED WHEN index > 255 or > 16 Bit Value   --> CMD ID 4 or 5!!!
        enccmds += struct.pack('<BB', 3, d)
        for c in dsf.Polygons[d]:
            if c[0] != poolIndex:
                fcn.write("1, %s\n" % str(c[0]))
                enccmds += struct.pack('<BH', 1, c[0])
                poolIndex = c[0]
            fcn.write(", ".join(str(n) for n in c[1:]))
            fcn.write("\n")
            if c[1] == 13:
                enccmds += struct.pack('<BHHH', 13, c[2], c[3], c[4])
            elif c[1] == 15:
                enccmds += struct.pack('<BHB', 15, c[2], len(c)-4) ### seems that n is number of windings and not indices, so one less #### ERROR in case length > 255
                for i in range(3, len(c)-3):
                    enccmds += struct.pack('<H', c[i])
            else:
                print("ERROR!!! Unknown COMMAND ID, ", c[1], "!!!!! CMDS ID 12 an 14  not yet implemented!!!!!!!!!!!")
    
    fcn.write("3, 0\n") ### Currently there is only one Road-Defintion --> DefIndex set to 0
    for d in dsf.Networks:
        if d[0][0] != subroadtype:
            fcn.write("6, %s\n" % str(d[0][0]))
            subroadtype = d[0][0]
        if d[0][1] != junctionoffset:
            fcn.write("2, %s\n" % str(d[0][1]))
            junctionoffset = d[0][1]
        if d[0][2] != poolIndex:
            fcn.write("1, %s\n" % str(d[0][2]))
            poolIndex = d[0][2]
        for c in d[1:]:           
            fcn.write(", ".join(str(n) for n in c))
            fcn.write("\n")
        
    ###print(dsf._Atoms_['SDMC'])
    print("-----")
    dsf._Atoms_['SDMC'] = enccmds
    ###print(dsf._Atoms_['SDMC'])
    
###dsf.write("CMDS_Test.dsf")
