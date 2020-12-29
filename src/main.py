#!/usr/bin/env python
# coding: utf-8
import numpy as np
import copy
from Bio.PDB import *
from glob import glob
import os
import sys
import warnings
from funPDB import *
from geometricks import *
from PlumedCV import *
from pathCV_S_Z import *
from interface import *


#ResApolar = ['VAL','ILE','LEU','MET','PHE','TRP']
#ResPolar = ["ASN","ASP","GLU","GLN","ARG","LYS","SER","THR", "TYR", "CYS"]



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
    parser.add_argument('-d', action="store", dest="d", type=str, default = 5,
                    help="Distance maximum between the chains, default = 5")
    parser.add_argument('-vector', dest="v", type=str, default = "-1 0 0",
                    help="vector for the translation by default the translation is \
    made on X axis. It is possible to provide other vector. \
    Example: -v \"-1 0 0\", default = \"-1 0 0\"")
    parser.add_argument('-method', action="store", dest="m", type=str, default = "MSD", \
    help="The method use to evaluate the distance between frame (Default is MSD).\n\
RMSD/MSD are available and DRMSD: ALL/INTER-DRMSD/INTRA-DRMSD chain.\n\
    /!\ To increase speed with DRMSD use: INTER-DRMSD (faster)")
    parser.add_argument('-distmax', action="store", dest="distmax", type=float, default = 3.2, \
    help="Distance maximum between the chains, default = 3.2")
    parser.add_argument('-eps', action="store", dest="eps", type=float, default = 0.0001, \
    help="The epsilon precision range (value +/- accuracy)")
    parser.add_argument('-stride', action="store", dest="s", type=float, default = 0.04,
                    help="stride: Distance value between structures, default = 0.04")
    parser.add_argument('-log', action="store", dest="log", type=str,\
     default = "file.log", help="log file's name: default file.log")
    parser.add_argument('-o', action="store", dest="o", type=str, default = "pathMSD.pdb",\
    help="output path filename ")
    parser.add_argument('-unit', action="store", dest="u", type=str, default = "nm", \
    help="Length unit /!\ Use the same as your MD engine for a correct lambda, default = nm")


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
    #print(interfaceList_ref)
    #print(interfaceList_ref[0])
    #print(interfaceList_ref[1])
    receptorCA, receptorH, receptorApo, CoMReceptor, interfR = SelectAtoms(reference[arg.r],interfaceList_ref[0])

    #For the ligand
    ligandCA, ligandH, ligandApo, CoMLigand, interfL = SelectAtoms(\
                                            reference[arg.l],\
                                            interfaceList_ref[1])
    mdl = ParsePDB(filename, "nm'")
    receptorH, ligandH = contacPair(mdl, receptorH, ligandH)
    structure4pathCV(filename, np.array(interfR),np.array(interfL))

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


    #sys.exit()
    #pdb = arg.g
    lig = arg.l
    pdb = "refStructure.pdb"
    rec = arg.r
    vector = np.fromstring(arg.v, dtype=float, sep=' ')
    method = arg.m
    distmax= arg.distmax
    eps = arg.eps
    unit = arg.u
    stride = arg.s

    mdl = ParsePDB(pdb, unit)
    if method == "MSD":
        print("The epsilon precision range is {0} {1}²".format(eps, unit))
        print("vector\t{0} ({1}²)".format(method, unit))
    else:
        print("The epsilon precision range is {0} {1}".format(eps, unit))
        print("vector\t{0} ({1})".format(method, unit))

    model, drmsd, transI = Vector4DRMSD(mdl, lig, rec, vector, eps, method, stride, distmax)
    exportStructure(model, arg.o, method, unit)
    exortlog(arg, transI, drmsd)
    lam = CV_S(model,drmsd, method, lig, rec, arg.log)
    print("\n\nZ values estimation when the ligand references's structures are translated \
and rotated.\n\
It is an indication to limit exploration on the Z axis. ┐(￣ヮ￣)┌ ")
    print("\nTranslation\tRotate\tZ mean\tZ sd\n\t(A)\t(°)")
    for j in range(0,16,5):
        translate_mdl = Z_wall(drmsd, transI, unit, copy.deepcopy(model), lig, lam, j)
        Zlim_Estimation(model, translate_mdl, unit, lig, rec, lam, j)
    #Generate input for plumed
    #list  atoms compose reference structure
    #number of atoms in fragA
    nb_atm = np.shape(mdl.atoms[rec])[0]
    fragA=""
    for i in range(nb_atm):
        fragA += mdl.atoms[rec][i].split()[1]+","
    fragA = fragA[:-1]
    nb_atm = np.shape(mdl.atoms[lig])[0]
    fragB=""
    for i in range(nb_atm):
        fragB += mdl.atoms[lig][i].split()[1]+","
    fragB = fragB[:-1]
    plumed_common([receptorCA, ligandCA], [fragA,fragB], [receptorH, ligandH], [receptorApo, ligandApo], lam)
    plumed_dat()
