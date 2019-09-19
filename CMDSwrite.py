from xplnedsf import *
import copy

dsf = XPLNEDSF()
dsf.read("CMDS_Test_in.dsf") ##+50-126.dsf   CMDS_Test_in.dsf

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
    enccmds = b'' 
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
                for i in range(3, len(c)): 
                    enccmds += struct.pack('<H', c[i])
            else:
                print("ERROR!!! Unknown COMMAND ID, ", c[1], "!!!!! CMDS ID 12 an 14  not yet implemented!!!!!!!!!!!")
    
    
    enccmds += struct.pack('<BL', 2, 0) #In X-Plane original dsf files command 2 comes first
    fcn.write("2, 0\n")
    junctionoffset = 0
    fcn.write("3, 0\n") ### Currently there is only one Road-Defintion --> DefIndex set to 0
    enccmds += struct.pack('<BB', 3, 0)
    for d in dsf.Networks:
        if d[0][1] != junctionoffset:
            fcn.write("2, %s\n" % str(d[0][1]))
            enccmds += struct.pack('<BL', 2, d[0][1])
            print(2,d[0][1])
            junctionoffset = d[0][1]
        if d[0][2] != poolIndex:
            fcn.write("1, %s\n" % str(d[0][2]))
            enccmds += struct.pack('<BH', 1, d[0][2])
            print(1,d[0][2])
            poolIndex = d[0][2]
        if d[0][0] != subroadtype:   ## Order in org X-Plane files is with 6 at last
            fcn.write("6, %s\n" % str(d[0][0]))
            enccmds += struct.pack('<BB', 6, d[0][0])
            print(6,d[0][0])
            subroadtype = d[0][0]
        for c in d[1:]:           
            fcn.write(", ".join(str(n) for n in c))
            fcn.write("\n")
            if c[0] == 9:
                enccmds += struct.pack('<BB', 9, len(c)-1) ## ToDo: error if this len is bigger than 255
                print (9, len(c)-1, c[1:])
                for i in range(1, len(c)):
                    enccmds += struct.pack('<H', c[i])  
            elif c[0] == 10:
                enccmds += struct.pack('<BHH', 10, c[1], c[2])
                print (10, c[1], c[2])
            elif c[0] == 11:
                enccmds += struct.pack('<BB', 11, len(c)-1) ## ToDo: error if this len is bigger than 255
                print("CMD 11 should not be here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                for i in range(1, len(c)):
                    enccmds += struct.pack('<L', c[i]) 
    patchcmds = bytearray()
    for d in dsf.Patches:
        fcn.write("3, %s\n" % str(d.defIndex))  ##### Only required when defIndex changes  ##### TO BE CHANGED WHEN index > 255 or > 16 Bit Value   --> CMD ID 4 or 5!!!
        #enccmds += struct.pack('<BB', 3, d.defIndex)
        patchcmds.extend(struct.pack('<BB', 3, d.defIndex)) #### IMPORTANT: Change of DefIndex needs to be directly before new patch, otherwise error
        fcn.write("1, %s\n" % str(d.cmds[0][1])) 
        patchcmds.extend(struct.pack('<BH', 1, d.cmds[0][1])) ## In X-Plane original files setting pool comes befor the patch definition
        fcn.write("18, %s, %s, %s\n" % (str(d.flag), str(d.near), str(d.far))) #### IF flags stay the same use CMD ID 16 or 17
        #enccmds += struct.pack('<BBff', 18, d.flag, d.near, d.far)
        patchcmds.extend(struct.pack('<BBff', 18, d.flag, d.near, d.far))
        ###### Currently poolIndex is always written at beginning of each patch as first CMD in d.cmds --> only requied if it changes
        for c in d.cmds[1:]:   ##skip first command as this is pool defintion written above
            fcn.write(", ".join(str(n) for n in c))
            fcn.write("\n")
            ecmd = b'' 
            id = c[0] #id of CMD stored as first value in list
            ecmd += struct.pack('<B', id) #encode CMD id
            #fcn.write("%s " % str(id))
            i = 1 #index in value list to be packed
            if id in dsf._CMDStructure_: ##### dsf. to be changed in self.  !!!!! ALSO BELOW !!!!!! #####
                for valtype in dsf._CMDStructure_[id][0]: #pack fixed length values of CMD
                    ecmd += struct.pack('<' + valtype, c[i])
                    #fcn.write("%s " % str(c[i]))
                    i += 1
                if len(dsf._CMDStructLen_[id]) == 3: #pack command with variable length 
                    vlength = int( (len(c) - i) / len(dsf._CMDStructure_[id][2]) ) #fixed values inclding CMD id are counting not for the variable length, so minus i #### TBD: Test that this value does not exceed 255!!!!!! ####
                    if id == 15:
                        vlength -= 1 #id = 15 seems a special case that there is one index more than windings  ########??????????
                    ecmd += struct.pack('<B', vlength) #count packed
                    #fcn.write("L:%s " % str(vlength))
                    while i < len(c): #now pack the variable length value list
                        for valtype in dsf._CMDStructure_[id][2]: #pack fixed length individual value
                            ecmd += struct.pack('<' + valtype, c[i])
                            #fcn.write("%s " % str(c[i]))
                            i += 1
                #fcn.write("\n")
            else:
                print("CMD ID not supported !!!!!!!!!!!!!!!!!!!")
            patchcmds.extend(ecmd)

        
    ###print(dsf._Atoms_['SDMC'])
    print("-----")
    dsf._Atoms_['SDMC'] = enccmds + patchcmds
    ###print(dsf._Atoms_['SDMC'])
    
dsf.write("CMDS_Test.dsf")
