# coding: utf-8

import numpy as np



def normalize(tmpvec):
    """
    take a vector and return it with a norme = 1
    Argument:
        _tmpvec (vector)
    Return:
        vector
    """
    norme = np.sqrt(np.sum(np.square(tmpvec)))
    return tmpvec/norme

def NormalVector4Vector(vec):
    """
    take a vector and return its 2 normal vector
    Argument:
        _tmpvec (vector)
    Return:
        _vector
    """
    vec2 = []
    vec3 = []
    if np.sum([vec[1],-vec[0],0]) != 0:
        vec2 = normalize([vec[1],-vec[0],0])
    else:
        print(np.sum([vec[1],-vec[0],0]))
        vec2 = normalize([vec[0],-vec[2],0])
    if np.sum([vec[2],0,-vec[0]]) != 0:
        vec3 = normalize([vec[2],0,-vec[0]])
    else:
        print(np.sum([vec[2],0,-vec[0]]))
        vec3 = normalize([vec[0],0,-vec[1]])
    return vec2,vec3


def rotate_matrix(mat, angle):
    """
    The function center a matrix to zero, then applies a generic rotation
    Arguments:
        mat: coordinates
        angle: angle for the rotation (the right-hand rule)
        axis: which axis the rotation is made (1=x;2=y;3=z)
    Retrun:
        mat
    """
    thetaX = np.radians(angle[0])
    thetaY = np.radians(angle[1])
    thetaZ = np.radians(angle[2])
    newMat=[]
    CoM = np.mean(mat,0)
    #Rotate on X axis
    Rx = np.array((
                    (1, 0, 0),
                    (0,  np.cos(thetaX), -np.sin(thetaX)),
                    (0,  np.sin(thetaX),  np.cos(thetaX))
                    ))
    #Rotate on Y axis
    Ry = np.array((
                    (np.cos(thetaY), 0, np.sin(thetaY)),
                    (0,  1, 0),
                    (-np.sin(thetaY),  0,  np.cos(thetaY))
                    ))
    #Rotate on Z axis
    Rz = np.array((
                    (np.cos(thetaZ), -np.sin(thetaZ), 0),
                    (np.sin(thetaZ),  np.cos(thetaZ), 0),
                    (0,  0,  1)
                    ))
    #print('apply the rotation matrix to structure: r*v')
    #print( Rx.dot((np.transpose(mat)) ))
    #center to origin and apply translation
    #Rotate matrix must be made in two step
    #R np.dot(Rx,Ry,Rz) doesn't work
    R = np.dot(Rx, Ry)
    R = np.dot(R, Rz)
    return np.transpose( R.dot(np.transpose(mat-CoM)))+CoM



def computeDistance(elmt1, mat):
    return np.sqrt(np.sum(np.square(elmt1 - mat), 1))


def contacPair(mdl, receptorH, ligandH):
    """
    Compute distance between list atoms
    and return a list of contact pair
    Arguments:
        mdl: Model class
        receptorH: list1 atoms
        ligandH: list2 atoms
    """
    listRec = []
    listLig = []
    for indice in receptorH:
        tmp = computeDistance(mdl.MatrixCoordinate()[indice], mdl.MatrixCoordinate()[ligandH])
        #Keep atoms distant bellow 7 ang
        ind_pair = np.where(tmp < 7 )[0]
        if len(ind_pair) > 0:
            for elmt in ind_pair:
                listRec.append(indice)
                listLig.append(ligandH[elmt])
    return listRec, listLig
