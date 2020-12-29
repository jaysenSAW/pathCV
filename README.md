# Author: Jaysen SAWMYNADEN

#                               LICENCE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the [GNU General Public License](http://www.gnu.org/licenses/)
along with this program.

##                               Principle

It generates a pathMSD.pdb file countaining all the models. By default the Mean Square
Deviation (MSD) is computed from chain A (receptor) and B (ligand),
each model are distant of 0.04 nm² (4 A²) equivalent to a RMSD of 0.2 nm (2 A).



##                               Install and execute
To use it, first install all the requiered module:
```
numpy
```
By default, the id of receptor and ligand chain are respectively A 
and B. The option -r and -l allow change those default values.
The translation vector is [1,0,0], you can provide other vector with the argument -vector.
The list of all the options and their default values are descirbed below.
```
 -g Reference PDB
 
```
 Arguments bellow are optional
```
 -o: PDB with all the models conserved (pathMSD.pdb)
 -r: ID of receptor chain (default A)
 -l: ID of ligand chain (default B)
 -unit: Length unit (default nm)
 /!\ use the same as your MD engine for a correct lambda
 -vector: vector use to translate receptor ("1 0 0")
 -distmax: Maximum float use to multipliate the vector (3.2 nm)
 -method: How distance is computed RMSD, RMS or DRMSD: (default MSD)
 /!\ To increase speed with DRMSD use: INTER-DRMSD (faster)
 ```
 more details about the computation are available [here](https://plumed.github.io/doc-v2.5/user-doc/html/_d_r_m_s_d.html)
 ```
 -stride: metric value between models (0.04)
 -eps: The epsilon precision range (value +/- accuracy). By default use 0.001
 -log: Log file's name (file.log)
```

##                              Example
You can execute the script by using the command line:
```
python pathCV_S_Z.py -g reference.pdb
```
To change receptor and ligand respectively by A and E:
```
python pathCV_S_Z.py -g reference.pdb -r B -l E
```
To specify the translation vector
```
python pathCV_S_Z.py -g reference.pdb -vector "0 0 1"
```
##                               Output

2 files:
* Log file (arguments used and other info)
* PDB file contain the structures
Length unit is (and MUST be) ALWAYS in angstrom whatever your MD engine

#                               Warning

##                              RMSD/MSD

RMSD/MSD results can differt (a little bit) of Plumed. The main cause is the method to align structure [before](https://plumed.github.io/doc-v2.3/user-doc/html/_r_m_s_d.html).
The python script process only a SIMPLE alignment, whereas PATHMSD made a OPTIMAL alignment before to compute RMSD/MSD.

##                              Length unit

Be sure length unit used by your MD engine is the same specify by the flag -unit.
Indeed S and Z depends on a [lambda factor](https://plumed.github.io/doc-v2.5/user-doc/html/_p_a_t_h.html). This factor is obtain from the distance between reference frames.
By default the script uses nm unit (like GROMACS) and method MSD.
If your MD engine uses angstrom unit, you MUST change these defaults values:

Example with a MSD method:
* -stride 4
* -distmax 32
* -unit A

##                              UPCOMING FEATURES

* PEP 8
* Create curve path between state A and state B
* optimal alignment
* Add method rotation to class Model
* Compatibility with Python 3.X
