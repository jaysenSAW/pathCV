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
import sys


def Vector4DRMSD(mdl, lig, rec, vector,  eps, compute ="MSD", stride = 2.5, Distmax = 3.2):
    """
    Give the vector translation to create a regular space DRMSD
    Arguments:
        _mdl: Model object
        _lig: chain we translate
        _compute: ALL/INTER-DRMSD/INTRA-DRMSD chain (MORE DETAILS)
        https://plumed.github.io/doc-v2.5/user-doc/html/_d_r_m_s_d.html
        _stride: DRMSD value between structures
        _Distmax: maximum value for the translation
        _eps: RMSD's tolerence
    Return:
        _model: list of Model
        _drmsd: list countaining drmsd values
        _transI: list countaining vector values for the translation
    """
    drmsd = []
    ind = 0
    model = [mdl]
    transI = []
    for    i in np.arange(0.2, Distmax,Distmax/15):
        tmp = copy.deepcopy(mdl)
        tmp.chains[lig] = tmp.chains[lig] + i*vector
        val = DRMSD(model[ind], tmp, compute, lig, rec)
        if val < stride - eps:
            #print "{0} A\t{1} A".format(i*vector, round(val,3))
            continue
        elif val < stride + eps:
            model.append(tmp)
            drmsd.append(val)
            transI.append(i*vector)
            print("Find a structure")
            print("{0}\t{1}".format(transI[ind], round(drmsd[ind],4)))
            ind += 1
        elif val > stride + eps:
            #print "{0} A\t{1} A".format(i*vector, round(val,3))
            while val > stride + eps:
                i = i - eps
                tmp = copy.deepcopy(mdl)
                tmp.chains[lig] = tmp.chains[lig] + i*vector
                val = DRMSD(model[ind], tmp, compute, lig, rec)
                #print "{0} A\t{1} A".format(i*vector, round(val,3))
            model.append(tmp)
            drmsd.append(val)
            transI.append(i*vector)
            print("Find a structure")
            print("{0}\t{1}".format(transI[ind], round(drmsd[ind],4)) )
            ind += 1
    return model, drmsd, transI


def RotateXYZ(model, angle):
    """
    Rotate a structure in function of X, Y or Z axis
    Arguments:
        model: coordinate matrix
        _angle: rotation  angle (degree)
    """
    rotate = []
    tmp = copy.deepcopy(model)
    tmp.chains[lig] = rotate_matrix(tmp.chains[lig], [angle,0,0])
    rotate.append(tmp)
    tmp = copy.deepcopy(model)
    tmp.chains[lig] = rotate_matrix(tmp.chains[lig], [0,angle,0])
    rotate.append(tmp)
    tmp = copy.deepcopy(model)
    tmp.chains[lig] = rotate_matrix(tmp.chains[lig], [0,0,angle])
    rotate.append(tmp)
    return rotate


def computeLambda(model, sumDist):
    """
    Arguments:
        _model: references structures
        _sumDist: sum of all the distances
    """
    return round(2.3*(len(model)-1)/sumDist,3)


def give_S(model, lam, sumDist):
    """
    Arguments:
        _model: references structures
        _lam: lambda factor
        _sumDist: sum of all the distances
    """
    indice = range(1,len(model)+1)
    numer = np.sum(indice*np.exp(-lam*np.array(sumDist)))
    denumer = np.sum(np.exp(-lam*np.array(sumDist)))
    return numer/denumer


def Z_wall(drmsd, transI, unit, model, lig, lam, j):
	"""
	Compute normal vector and translate reference structures to estimate Z score
	and proposal values for a upperwall
	Argument:
		_drmsd: List countaining drmsd values
		_transI: List countaining vector values for the translation
		_model: List of Model
		_unit: Unit use by the MD engine
		_lig: Ligand's chain
		_Lam: Lambda factor (float) to compute S and Z score
		_j: Multiply unit vector by j (float)
	Return:
		_translate_mdl: List of Model
	"""
	translate_mdl = []
	if unit == "nm":
		factor = 0.1
	else:
		factor = 1
	#Compute the normal vector
	for i in range(1,len(transI)):
		vec = transI[i] - transI[i-1]
		vec2, vec3 = NormalVector4Vector(vec)
		tmp = copy.deepcopy(model[i])
		tmp.chains[lig] = tmp.chains[lig]+vec2*j*factor
		translate_mdl.append(tmp)
		tmp = copy.deepcopy(model[i])
		tmp.chains[lig] = tmp.chains[lig]+vec3*j*factor
		translate_mdl.append(tmp)
	return translate_mdl



def Z_score(lam, model, translate_mdl, lig, rec, method="MSD"):
	"""
	Compute Z score
	arguments:
		_lam: lambda value
		_model: list of reference structures
		_translate_mdl: list of structures
		_method: How distance is computed
		_lig: Ligand's chain
		_rec: Receptor's chain
	"""
	Zscore = []
	for i in range(len(translate_mdl)):
		tmpZ = 0
		for j in range(len(model)):
			tmpZ = tmpZ + np.exp(-lam*DRMSD(model[j],
			translate_mdl[i],
			method, lig, rec)
			)
		Zscore.append(-1/lam * np.log(tmpZ))
	return Zscore


def Zlim_Estimation(model, translate_mdl, unit, lig, rec, lam, j, method="MSD"):
    """
    Estimate the Z values after translated and rotated the ligand references's\
    structures. It is an indication to limit exploration on the Z axis.
    Arguments:
        _model:
        _unit: unit use to compute the distance (nm/A)
        _lig: ligand's chain
        _lam= Lambda value
        _j: Multply vector by j and translate references structure to an orthogonal
        direction.
    """
    #print "Translation (A)\tRotation (°)\tZ mean\tZ sd"
    #tmp_model = copy.deepcopy(model)
    #translate_mdl = Z_wall(drmsd, transI, unit, tmp_model, lig, lam, j)
    for angle in range(10,35,10):
        rotate = []
        rotate2 = []
        tmp = 0
        for i in range(len(translate_mdl)):
#            rotate2.append(RotateXYZ(translate_mdl[i],angle))
            tmp = copy.deepcopy(translate_mdl[i])
            tmp.chains[lig] = rotate_matrix(tmp.chains[lig], [angle,0,0])
            rotate.append(tmp)
            tmp = copy.deepcopy(translate_mdl[i])
            tmp.chains[lig] = rotate_matrix(tmp.chains[lig], [0,angle,0])
            rotate.append(tmp)
            tmp = copy.deepcopy(translate_mdl[i])
            tmp.chains[lig] = rotate_matrix(tmp.chains[lig], [0,0,angle])
            rotate.append(tmp)
        #exportStructure(rotate, "z_trans"+str(j)+"ang_rotate_"+str(angle)+".pdb", method, unit)
        Zval = Z_score(lam, model, rotate, lig, rec, method)
        print("\t{0}\t{1}\t{2:.3f}\t{3:.3f}".format(j,angle,np.mean(Zval), np.std(Zval)))



def exortlog(parameters, transI, drmsd):
    """
    Generate a log file with the parameters used and summary results
    Arguments:
        _parameters: argparse arguments
        _tansI: list countaining all vectors for the translation
        _list countaining drmsd values
    """
    filelog = parameters.log
    with open(filelog, 'w') as output:
        output.write("Summary ( ͡° ͜ʖ ͡°)\n")
        output.write("Reference PDB {0}\n".format(parameters.g))
        output.write("Use chain {0} (receptor) and {1} (ligand) to compute \
{2}\n".format(parameters.r, parameters.l, parameters.m))
        output.write("All option given (See the help for more details)\n{0}"\
.format(parameters))
        output.write("\n\n################\nSUMMARY results (☞ﾟヮﾟ)☞\n")
        output.write("{0}\tvector applyed from the reference\
 structure\n".format(parameters.m))
        for i in range(len(transI)):
            output.write("{0}\t{1}\n".format(round(drmsd[i],3),transI[i]))




def CV_S(model,drmsd, method, lig, rec, filelog, unit = "nm"):
    """
    Compute and display the S value and Z value
    """
    if method == "MSD":
        unit = unit+"²"
    with open(filelog, 'a+') as output:
        denulam = 0
        #print "\n\n\nCompute the distance between each ref structure (method: \
#{0})".format(method)
        for i in range(len(model)-1):
            tmp = DRMSD(model[i], model[i+1], method, lig, rec)
            #print "Dist frame {0}-{1} = {2}".format(i,i+1,round(tmp,3))
            denulam += tmp
        lam = computeLambda(model, denulam)
        print("\n\n\n#############\nLambda = {0} ({1})^-1".format(lam, unit))
        print("\nframe\ts\tz")
        output.write("\n\n\n#############\nLambda value {0} ({1})^-1\n\
/!\ lambda depend on the length unit use by your MD engine\n"\
.format(lam, unit))
        output.write("\n\nframe\ts\tz\n")
        for i in range(len(model)):
            numer=0
            denumer=0
            tmp = []
            for j in range(len(model)):
                val = DRMSD(model[i], model[j], method, lig, rec)
                tmp.append(val)
            s = give_S(model, lam, tmp)
            indice = range(1,len(model)+1)
            numer = np.sum(indice*np.exp(-lam*np.array(tmp)))
            denumer = np.sum(np.exp(-lam*np.array(tmp)))
            z = (-1/lam)*np.log(denumer)
            print("{0}\t{1}\t{2}".format(i+1,round(s,2),round(z,2)))
            output.write("{0}\t{1}\t{2}\n".format(i+1,round(s,3),round(z,3)))
        if method == "RMSD" or method == "MSD":
            output.write("\n\n\n#############\n{0} Matrix distance {1}\n\
#############\n".format(method,unit))
            #header
            output.write("Frame")
            tmp=0
            for i in range(len(model)):
                output.write("\t{0}".format(str(i)))
            #Write distance matrix
            for i in range(len(model)):
                output.write("\n{0}".format(str(i)))
                for j in range(len(model)):
                    tmp = DRMSD(model[i], model[j], method, lig, rec)
                    output.write("\t{0}".format(round(tmp,4)))
    return lam



def new_model(tmpmodel, vector):
    tmpmodel.chains[lig] = tmpmodel.chains[lig]+vector
    return tmpmodel


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
    #pdb = arg.g
    lig = arg.l
    pdb = arg.g
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
