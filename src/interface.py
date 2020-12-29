#!/usr/bin/env python
# coding: utf-8
import numpy as np
import copy
from funPDB import *
from geometricks import *
from PlumedCV import *
import warnings
from Bio.PDB import *
from glob import glob
import os



#####################################
#             Extract info            #
#####################################

def hydrogenBond(resname, res, tmp):
#######################################################################
#http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/#
    Output = []
    if resname == "ASN":
        Output = [tmp[res]["OD1"].get_serial_number()]
        Output.append(tmp[res]["ND2"].get_serial_number())
    if resname == "ASP":
        Output = [tmp[res]["OD1"].get_serial_number()]
        Output.append(tmp[res]["OD2"].get_serial_number())
    if resname == "GLN":
        Output = [tmp[res]["OE1"].get_serial_number()]
        Output.append(tmp[res]["NE2"].get_serial_number())
    if resname == "GLU":
        Output = [tmp[res]["OE1"].get_serial_number()]
        Output.append(tmp[res]["OE2"].get_serial_number())
    if resname == "SER":
        [tmp[res]["OG"].get_serial_number()]
    if resname == "THR":
        Output.append(tmp[res]["OG1"].get_serial_number())
    if resname == "ARG":
        Output = [tmp[res]["NE"].get_serial_number()]
        Output.append(tmp[res]["NH1"].get_serial_number())
        Output.append(tmp[res]["NH2"].get_serial_number())
    if resname == "HIS":
        Output = [tmp[res]["ND1"].get_serial_number()]
        Output.append(tmp[res]["NE2"].get_serial_number())
    if resname == "LYS":
        Output.append(tmp[res]["NZ"].get_serial_number())
    if resname == "TYR":
        Output.append(tmp[res]["OH"].get_serial_number())
    if resname == "CYS":
        Output.append(tmp[res]["SG"].get_serial_number())
    return Output



def Apolar(resname, res, tmp):
    #['VAL','ILE','LEU','MET','PHE','TRP']
    Output = []
    if resname == "VAL":
        Output = [tmp[res]["CB"].get_serial_number()]
        Output.append(tmp[res]["CG1"].get_serial_number())
        Output.append(tmp[res]["CG2"].get_serial_number())
        return Output
    if resname == "ILE":
        Output = [tmp[res]["CB"].get_serial_number()]
        Output.append(tmp[res]["CG1"].get_serial_number())
        Output.append(tmp[res]["CG2"].get_serial_number())
        Output.append(tmp[res]["CD"].get_serial_number())
        return Output
    if resname == "LEU":
        Output = [tmp[res]["CB"].get_serial_number()]
        Output.append(tmp[res]["CG"].get_serial_number())
        Output.append(tmp[res]["CD1"].get_serial_number())
        Output.append(tmp[res]["CD2"].get_serial_number())
        return Output
    if resname == "MET":
        Output = [tmp[res]["CB"].get_serial_number()]
        Output.append(tmp[res]["CG"].get_serial_number())
        Output.append(tmp[res]["CE"].get_serial_number())
        return Output
    if resname == "PHE":
        Output = [tmp[res]["CB"].get_serial_number()]
        Output.append(tmp[res]["CG"].get_serial_number())
        Output.append(tmp[res]["CD1"].get_serial_number())
        Output.append(tmp[res]["CZ"].get_serial_number())
        Output.append(tmp[res]["CE1"].get_serial_number())
        Output.append(tmp[res]["CE2"].get_serial_number())
        Output.append(tmp[res]["CD2"].get_serial_number())
        return Output
    if resname == "TRP":
        Output = [tmp[res]["CB"].get_serial_number()]
        Output.append(tmp[res]["CG"].get_serial_number())
        Output.append(tmp[res]["CD1"].get_serial_number())
        Output.append(tmp[res]["CE2"].get_serial_number())
        Output.append(tmp[res]["CZ2"].get_serial_number())
        Output.append(tmp[res]["CH2"].get_serial_number())
        Output.append(tmp[res]["CZ3"].get_serial_number())
        Output.append(tmp[res]["CE3"].get_serial_number())
        Output.append(tmp[res]["CD2"].get_serial_number())
        return Output


def getCommonAtoms(chainOne, chainTwo, atomNames, subset = range(-9999, 9999)):
    """
    Get common atoms
    """
    fixed = []
    moving = []
    for resID in subset:
        if chainOne.has_id(resID) and chainTwo.has_id(resID):
            resOne = chainOne[resID]
            resTwo = chainTwo[resID]
            if resOne.get_resname() == resTwo.get_resname():
                for atomName in atomNames:
                    if resOne.has_id(atomName) and resTwo.has_id(atomName):
                        fixed.append(resOne[atomName])
                        moving.append(resTwo[atomName])
            else:
                print("Skipped residue ", resID, " as amino acid type does not match: ", resOne.get_resname(), " vs. ", resTwo.get_resname())
    return (fixed, moving)


def getNeighbors(filename, maxRadius):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    p = PDBParser(PERMISSIVE=1)
    s = p.get_structure("complex", filename)
    firstModel = s[0]
    atom_list = Selection.unfold_entities(firstModel, 'A')
    neigh = NeighborSearch(atom_list)
    allNeigh = neigh.search_all(maxRadius, 'R')
    return allNeigh


def getInterfaceDic(allNeigh, receptorID, ligandID):
    neighbourDicA = {}
    neighbourDicB = {}
    for n in allNeigh:
        resOne = n[0].get_id()[1]
        resTwo = n[1].get_id()[1]
        resNameOne = n[0].get_resname()
        resNameTwo = n[1].get_resname()
        if (not (resNameOne in Polypeptide.standard_aa_names) or not (resNameTwo in Polypeptide.standard_aa_names)):
            continue
        chainOne = "X"
        chainTwo = "X"
        chainOne = n[0].get_full_id()[2]
        chainTwo = n[1].get_full_id()[2]
    #    if resOne in surfaceRes and resTwo in surfaceRes:
        if (chainOne == receptorID and chainTwo == ligandID):
            if resOne in neighbourDicA:
                if not(resTwo in neighbourDicA[resOne]):
                    neighbourDicA[resOne].append(resTwo)
            else:
                neighbourDicA[resOne] = []
                neighbourDicA[resOne].append(resTwo)
            if resTwo in neighbourDicB:
                if not(resOne in neighbourDicB[resTwo]) :
                    neighbourDicB[resTwo].append(resOne)
            else:
                neighbourDicB[resTwo] = []
                neighbourDicB[resTwo].append(resOne)
        elif (chainOne == ligandID and chainTwo == receptorID):
            if resOne in neighbourDicB:
                if not(resTwo in neighbourDicB[resOne]):
                    neighbourDicB[resOne].append(resTwo)
            else:
                neighbourDicB[resOne] = []
                neighbourDicB[resOne].append(resTwo)
            if resTwo in neighbourDicA:
                if not(resOne in neighbourDicA[resTwo]) :
                    neighbourDicA[resTwo].append(resOne)
            else:
                neighbourDicA[resTwo] = []
                neighbourDicA[resTwo].append(resOne)
    return (neighbourDicA, neighbourDicB)


def getInterfaceList(interfaceDic):
    neighbourDicA = interfaceDic[0]
    neighbourDicB = interfaceDic[1]
    #print "neighbours from chain A to chain B:"
    interfaceA = []
    for (resA, resB) in neighbourDicA.items():
        interfaceA.append(resA)
    #print "neighbours from chain B to chain A:"
    interfaceB = []
    for (resA, resB) in neighbourDicB.items():
        interfaceB.append(resA)
    return (interfaceA, interfaceB)



def getContactList(interfaceDic):
    neighbourDicA = interfaceDic[0]
    contactList = []
    for (resA, resBlist) in neighbourDicA.items():
        for resB in resBlist:
            contactList.append((resA, resB))
    return contactList


ResApolar = ['VAL','ILE','LEU','MET','PHE','TRP']
ResPolar = ["ASN","ASP","GLU","GLN","ARG","LYS","SER","THR", "TYR", "CYS"]



def new_model(tmpmodel, vector):
	tmpmodel.chains[lig] = tmpmodel.chains[lig]+vector
	return tmpmodel



def SelectAtoms(ref,interfaceList):
    CA = []
    H =[]
    ApolarList = []
    CoM = [] #Center of Mass
    for res in interfaceList:
        CA.append(ref[res]['CA'].get_serial_number())
        resname = ref[res].get_resname()
        coord = ref[res]['CA'].get_vector()
        CoM.append([coord[0],coord[1],coord[2]])
        if resname in ResPolar:
            H = H + hydrogenBond(resname, res, ref)
        if resname in ResApolar:
            ApolarList = ApolarList + Apolar(resname, res, ref)
    CA.sort()
    H.sort()
    ApolarList.sort()
    tmp = np.matrix(CoM)
    CoM = [np.mean(tmp[:,0]),np.mean(tmp[:,1]),np.mean(tmp[:,2])]
    interf = CA + H + ApolarList
    interf.sort()
    return CA, H, ApolarList, CoM, interf



def createPymolScript(interfaceList, filename, receptorID, ligandID):
    interfaceA = interfaceList[0]
    interfaceB = interfaceList[1]
    outStringA = ""
    for resID in interfaceA:
        outStringA += str(resID)
        outStringA += ","
    outStringB = ""
    for resID in interfaceB:
        outStringB += str(resID)
        outStringB += ","
    pymolFile = filename.replace(".pdb", ".pymol")
    f = open(pymolFile, "w")
    f.write("cmd.load('%s')\n" % os.path.basename(filename))
    f.write("cmd.select('interfaceReceptor', '(chain " + receptorID + " and resi " + outStringA + ")')\n")
    f.write("cmd.select('interfaceLigand', '(chain " + ligandID + " and resi " + outStringB + ")')\n")
    f.write("cmd.hide('everything')\n")
    f.write("cmd.show('cartoon')\n")
    f.write("cmd.bg_color('white')\n")
    f.write("cmd.color('gray', 'chain " + receptorID + "')\n")
    f.write("cmd.color('yellow', 'chain " + ligandID + "')\n")
    f.write("cmd.color('red', 'interfaceReceptor')\n")
    f.write("cmd.color('blue', 'interfaceLigand')\n")
    f.write("cmd.select('not all')\n")
    f.close()


def OutputList(list):
    tmp=str(list[0])
    for elmt in list[1:]:
        tmp=tmp+","+str(elmt)
    return tmp


def structure4pathCV(filename, receptorH,ligandH):
    """
    Export proteins's interface to refStructure.pdb
    Arguments:
        filename: reference PDB file to read
        receptorH: list atoms's number to export
        ligandH: list atoms's number to export
    """
    with open("refStructure.pdb", "w") as filout:
        with open(filename, 'r') as filin:
            for line in filin:
                if len(line.split()) < 9:
                    continue
                elif line.split()[0] == "ATOM" and int(line.split()[1]) in receptorH:
                    filout.write(line)
                elif line.split()[0] == "ATOM" and int(line.split()[1]) in ligandH:
                    filout.write(line)



if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='Generate pdb file for pathCV \
and compute lambda, S and Z score.\n\
Example: python pathCV_S_Z.py -g input.pdb -o ref_structures.pdb -distmax 3.2 -unit nm')
	parser.add_argument('-g', action="store", dest="g", type=str,
		            help="Topolgy file (PDB).\nExample: reference.pdb")
	parser.add_argument('-r', action="store", dest="r", type=str,
		            help="id of receptor chain, default = A", default = "A")
	parser.add_argument('-l', action="store", dest="l", type=str, default = "B",
		            help="id of ligand chain, default = B")
	parser.add_argument('-d', action="store", dest="d", type=float, default = 5, \
	help="Distance maximum between the chains, default = 5")
	parser.add_argument('-log', action="store", dest="log", type=str,\
	 default = "file.log", help="log file's name: default file.log")
	parser.add_argument('-o', action="store", dest="o", type=str, default = "pathMSD.pdb",\
	help="output path filename ")


	arg = parser.parse_args()
	print("Parameters used (∩｀-´)⊃━☆ﾟ\n{0}\nSee the help for more details.\n\n"\
.format(arg))
	print("Results summary is store in {0}".format(arg.log))

	filename = arg.g
	allNeigh = getNeighbors(filename, arg.d)
	interfaceDic_ref = getInterfaceDic(allNeigh, arg.r, arg.l)
	contactList_ref = getContactList(interfaceDic_ref)
	interfaceList_ref = getInterfaceList(interfaceDic_ref)
	createPymolScript(interfaceList_ref, filename, arg.r, arg.l)
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
	p = PDBParser(PERMISSIVE=1)
	reference = p.get_structure("reference", filename)[0]

	receptor = reference[arg.r]
	ligand = reference[arg.l]
	ResApolar = ['VAL','ILE','LEU','MET','PHE','TRP']
	ResPolar = ["ASN","ASP","GLU","GLN","ARG","LYS","SER","THR", "TYR", "CYS"]

    #For the receptor
	receptorCA, receptorH, receptorApo, CoMReceptor = SelectAtoms(\
                                            reference[arg.r],\
                                            interfaceList_ref[0])
    #For the ligand
	ligandCA, ligandH, ligandApo, CoMLigand = SelectAtoms(\
                                            reference[arg.l],\
                                            interfaceList_ref[1])


	structure4pathCV(filename, np.array(receptorH),np.array(ligandH))


    #Save the Center of Mass of the 2 interfaces
	f = open(filename[:-4]+"_CoM.pdb", "w")
	f.write("TITLE     "+filename[:-4]+"\nMODEL\n")
	f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}\
	   {:8.3f}{:8.3f}{:8.3f}\
	\n".format("ATOM",1,"O"," ","COM","X", 1, " ",\
	CoMReceptor[0],CoMReceptor[1],CoMReceptor[2], 1.0, 50.0))
	f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}\
	   {:8.3f}{:8.3f}{:8.3f}\
	".format("ATOM",1,"O"," ","COM","X", 2, " ",\
	CoMLigand[0],CoMLigand[1],CoMLigand[2], 1.0, 50.0))
	f.write("\nTER")
	f.close()


	plumed_common([receptorCA, ligandCA], [receptorH, ligandH], [receptorApo, ligandApo])
	plumed_dat()
