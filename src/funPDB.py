# coding: utf-8
import numpy as np


class Model:
	"""Model for the complex
		Store the PDB's line and store the coordinate of
		each  chain's atoms
	"""
	def __init__(self, chain, elmt1, elmt2):
		self.List_Chain = chain
		self.chains = elmt1
		self.atoms = elmt2
	def getSerial(self, chain):
		#atoms's serial in chain
		tmp = []
		for i in range(len(self.atoms[chain])):
			tmp.append(int(self.atoms[chain][i].split()[1]))
		return tmp
	def MatrixCoordinate(self):
		keys = self.List_Chain
		tmp = self.chains[keys[0]]
		#concatenate all the coordinate
		for i in range(len(keys)-1):
			tmp = np.concatenate(
					(tmp, self.chains[keys[i+1]])
					)
		return tmp
	#################################################
	#			Methods for distance RMSD 			#
	#################################################
	#I.	Compute distances between all atoms
	def MatrixAllDist(self):
		flag = False
		keys = self.List_Chain
		tmp = self.chains[keys[0]]
		#concatenate all the coordinate
		for i in range(1,len(keys)):
			tmp = np.concatenate(
					(tmp, self.chains[keys[i]])
					)
		#compute the pairing distance
		for i in range(len(tmp)):
			for j in range(len(tmp)):
				if j == i:
					continue
				atm1 = tmp[i]
				atm2 = tmp[j]
				tmpdist = [np.sqrt(np.square(
						atm1[0]-atm2[0])+\
						np.square(atm1[1]-atm2[1])+\
						np.square(atm1[2]-atm2[2])\
					)]
				if flag is False:
					dist = tmpdist
					flag = True
				elif flag is True:
					dist = np.append(dist,tmpdist)
		return dist
	#II. Compute distance between atoms localizated in two distinc chains
	def MatrixInterDist(self, chain1, chain2):
		#Only for two chains distinc
		flag = False
		for i in range(len(self.chains[chain1])):
			for j in range(len(self.chains[chain2])):
				atm1 = self.chains[chain1][i]
				atm2 = self.chains[chain2][j]
				tmpdist = [np.sqrt(np.square(
						atm1[0]-atm2[0])+\
						np.square(atm1[1]-atm2[1])+\
						np.square(atm1[2]-atm2[2])\
					)]
				if flag is False:
					dist = tmpdist
					flag = True
				elif flag is True:
					dist = np.append(dist,tmpdist)
		return dist
	#III. Compute distance between atoms localizated in same chain
	def MatrixIntraDist(self, chain1, chain2):
		#Only for two chains distinc
		flag = False
		for i in range(len(self.chains[chain1])):
			for j in range(len(self.chains[chain1])):
				if i == j:
					continue
				atm1 = self.chains[chain1][i]
				atm2 = self.chains[chain1][j]
				tmpdist = [np.sqrt(np.square(
						atm1[0]-atm2[0])+\
						np.square(atm1[1]-atm2[1])+\
						np.square(atm1[2]-atm2[2])\
					)]
				if flag is False:
					dist = tmpdist
					flag = True
				elif flag is True:
					dist = np.append(dist,tmpdist)
		if chain1 == chain2:
			return dist
		else:
			for i in range(len(self.chains[chain2])):
				for j in range(len(self.chains[chain2])):
					if i == j:
						continue
					atm1 = self.chains[chain1][i]
					atm2 = self.chains[chain1][j]
					tmpdist = [np.sqrt(np.square(
							atm1[0]-atm2[0])+\
							np.square(atm1[1]-atm2[1])+\
							np.square(atm1[2]-atm2[2])\
						)]
					dist = np.append(dist,tmpdist)
		return dist
	def __str__(self):
		return "The model has {0} chain(s) and {1} atoms".format(\
				len(self.chains.keys()),  len( self.MatrixCoordinate() ) )


def ParsePDB(pdb, unit):
	"""
	Parse a PDB file and store the coordinate and atoms's id in the object mdl
	Argument
		_PDB: PDB file
	Return
		_modl: class Model
	"""
	if unit == "nm":
		factor = 0.1
	else:
		factor = 1
	flag = False
	with open(pdb, "r") as filin:
		for line in filin:
			if line == "ENDMDL\n":
				return modl
			elif len(line.split()) <=1 or line.split()[0] != "ATOM":
				continue
			key = line[21]
			if line.split()[0] == "ATOM" and flag is False:
				modl = Model(
							[key],
							{key: np.array([ line[30:54].split() ], dtype = np.float64)*factor},
							{key: [line[:30]]}
				)
				flag = True
			elif key not in modl.List_Chain and flag is True:
				modl.List_Chain.append(key)
				modl.chains[key] = np.array([ line[30:54].split() ], dtype = np.float64)*factor
				modl.atoms[key] = [line[:30]]
			elif key in modl.List_Chain:
				modl.chains[key] = np.concatenate(
							(
							modl.chains[key],
							np.array([ line[30:54].split() ], dtype = np.float64)*factor
							)
				)
				modl.atoms[key].append(line[:30])
	return modl


def XYZ(coord):
	"""
	Convert coordinate to string and check the good number of caracters
	Argument:
		_coord: float
	"""
	while len(str(coord).split(".")[1]) < 3:
		coord = str(coord)+"0"
	while len(str(coord).split(".")[1]) > 3:
		coord = str(coord)[:-1]
	while len(str(coord).split(".")[0]) < 3:
		coord = " "+str(coord)
	return str(coord)


def exportStructure(model, filout, method, unit="ang"):
	"""
	Export structures and save them into filout (PDB format)
	Argument:
		_model: list of Model
		_filout: pdb file's name
		_method: the method to compute the DRMSD
	"""
	if unit == "nm":
	#convert nm to angstrom for the pdb file
		factor = 10.0
	else:
		factor = 1.0
	with open(filout, 'w') as output:
		for i in range(len(model)):
			for chain in model[i].List_Chain:
				for ind in range(len(model[i].atoms[chain])):
					atomID = model[i].atoms[chain][ind]
					xyz = model[i].chains[chain][ind]*factor
					xyz = XYZ(xyz[0])+" "+XYZ(xyz[1])+" "+XYZ(xyz[2])
					output.write(atomID+" "+xyz+"  1.00  1.00\n")
				if method != "ALL" and method != "RMSD" and method != "MSD":
					output.write("TER\n")
			output.write("ENDMDL\n")
