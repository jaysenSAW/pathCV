# coding: utf-8
from geometricks import *
import copy
import numpy as np

def OutputList(list):
    tmp=str(list[0])
    for elmt in list[1:]:
        tmp=tmp+","+str(elmt)
    return tmp

def DRMSD(model1, model2, compute, lig, rec):
	"""
	Compute the RMSD or MSD between two models.
	"""
	if compute == "INTER-DRMSD":
		dist = model1.MatrixInterDist(rec, lig) - model2.MatrixInterDist(rec, lig)
		dist2 = np.square(dist)
		N = len(model1.MatrixCoordinate())
		#dN = N*(N-1)
		dN = len(dist)
		return np.sqrt(np.sum(dist2)/dN)
	elif compute == "ALL":
		dist = model1.MatrixAllDist() - model2.MatrixAllDist()
		dN = len(dist)
		dist2 = np.square(dist)
		return np.sqrt(np.sum(dist2)/dN)
	elif compute == "INTRA-DRMSD":
		dist = model1.MatrixIntraDist(rec, lig) - model2.MatrixIntraDist(rec, lig)
		dN = len(dist)
		dist2 = np.square(dist)
		return np.sqrt(np.sum(dist2)/dN)
	elif compute == "MSD" or compute == "RMSD":
		#Do a fit on the system as function of the center of mass (CoM)
		CoM1 = np.mean(model1.MatrixCoordinate(),0)
		CoM2 = np.mean(model2.MatrixCoordinate(),0)
		centerMdl1 = model1.MatrixCoordinate() - CoM1
		centerMdl2 = model2.MatrixCoordinate() - CoM2
		msd = np.square(centerMdl1-centerMdl2)
		#without fit
		#rmsd = model1.MatrixCoordinate() - model2.MatrixCoordinate()
		#rmsd = np.square(rmsd)
		#rmsd= np.sum(rmsd,1)
		if compute == "MSD":
			return np.sum(msd)/len(msd)
		else:
			return np.sqrt( np.sum(msd)/len(msd) )
	else:
		print("Method was not found !")


def plumed_common(CA, frag, H, Apo, lamb):
    print("\n\n\ninput for plumed-common.dat")
    f = open("plumed-common.dat", "w")
    f.write("RESTART\nRANDOM_EXCHANGES\nFLUSH STRIDE=100\n")
    f.write("\n\n#CoM Corbon alpha\n")
    f.write("CoM1: COM ATOMS="+OutputList(CA[0]) )
    f.write("\nCoM2: COM ATOMS="+OutputList(CA[1]) )
    f.write("\n\n#Atoms for the pathCV")
    f.write("\nInterfFragA: GROUP ATOMS="+frag[0])
    f.write("\nInterfFragB: GROUP ATOMS="+frag[1])
    f.write("\n\n#Hydrophobic residues" )
    f.write("\nHydroPhobicFragA: GROUP ATOMS="+OutputList(Apo[0]) )
    f.write("\nHydroPhobicFragB: GROUP ATOMS="+OutputList(Apo[1]) )
    f.write("\n\n#Atoms involved in hydrogen bonds")
    f.write("\nSaltA: GROUP ATOMS="+OutputList(H[0]) )
    f.write("\nSaltB: GROUP ATOMS="+OutputList(H[1]) )
    f.write("\n\ncomdist: DISTANCE ATOMS=CoM1,CoM2\n")
    f.write("WHOLEMOLECULES ENTITY0=InterfFragA ENTITY1=InterfFragB STRIDE=1\n")
    f.write("p1: PATHMSD REFERENCE=pathMSD.pdb LAMBDA="+str(lamb)+"\n")
    f.write("uwall: UPPER_WALLS ARG=p1.zzz AT=ValGiven KAPPA=100000.0\n")
    f.write("CoordHydropphobe: COORDINATION GROUPA=HydroPhobicFragA GROUPB=HydroPhobicFragB R_0=0.55 NN=8 MM=12\n")
    f.write("CoordPolar: COORDINATION GROUPA=SaltA GROUPB=SaltB  R_0=0.35 NN=8 MM=12 PAIR\n")
    f.write("w1: COORDINATION GROUPA=SaltA GROUPB=3104-84310:3 SWITCH={RATIONAL R_0=0.35 NN=8 MM=10 D_MAX=1}\n\n\n")
    f.close()

def plumed_dat():
    print("input for plumed.0.dat")
    f = open("plumed.0.dat", "w")
    f.write("INCLUDE FILE=plumed-common.dat\n\
    be: METAD ARG=p1.sss HEIGHT=0.5 SIGMA=0.15 PACE=8000\n\
    PRINT STRIDE=100 ARG=p1.sss,p1.zzz,CoordHydropphobe,CoordPolar,w1,comdist,be.bias,uwall.bias FILE=COLVAR")
    f.close()
    ###############################
    print("input for plumed.1.dat")
    f = open("plumed.1.dat", "w")
    f.write("INCLUDE FILE=plumed-common.dat\n\
    be: METAD ARG=CoordHydropphobe HEIGHT=0.5 SIGMA=1.5 PACE=4000 INTERVAL=1.5,110.0\n\
    PRINT STRIDE=100 ARG=p1.sss,p1.zzz,CoordHydropphobe,CoordPolar,w1,comdist,be.bias,uwall.bias FILE=COLVAR")
    f.close()
    ###############################
    print("input for plumed.2.dat")
    f = open("plumed.2.dat", "w")
    f.write("INCLUDE FILE=plumed-common.dat\n\
    be: METAD ARG=CoordPolar HEIGHT=0.5 SIGMA=0.3 PACE=4000 INTERVAL=0.3,25.0\n\
    PRINT STRIDE=100 ARG=p1.sss,p1.zzz,CoordHydropphobe,CoordPolar,w1,comdist,be.bias,uwall.bias FILE=COLVAR")
    f.close()
    ###############################
    print("input for plumed.3.dat")
    f = open("plumed.3.dat", "w")
    f.write("INCLUDE FILE=plumed-common.dat\n\
    be: METAD ARG=w1 HEIGHT=0.5 SIGMA=5.0 PACE=4000 INTERVAL=60,300.0\n\
    PRINT STRIDE=100 ARG=p1.sss,p1.zzz,CoordHydropphobe,CoordPolar,w1,comdist,be.bias,uwall.bias FILE=COLVAR")
    f.close()
