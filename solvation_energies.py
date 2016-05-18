"""
Perform single point energy calculations in solvent on a set of gas phase
TS geometries in a database.
"""
import sqlite3
import openbabel as ob
from subprocess import Popen
import os
import sys
from rmgpy.molecule import Molecule

if len(sys.argv)>1:
	i = int(sys.argv[-1])
elif os.getenv('LSB_JOBINDEX'):
	i = int(os.getenv('LSB_JOBINDEX'))
else:
	raise Exception("Specify a TS number!")

conn = sqlite3.connect('/scratch/westgroup/ts_data.db') # Table name is TSData
print "Opened database successfully"

family = 'intra_H_migration'
solvent = 'n-octane'

geometries = conn.execute("SELECT uniqueID, geometry FROM TSData WHERE rxnFamily=\'" + family + "\'") # returns cursor
geometries = list(geometries.fetchall()) # make into list
entry = geometries[i-1]
#for n, entry in enumerate(geometries):
# Convert the geometry to xyz
obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("cml", "xyz")
mol = ob.OBMol()
obConversion.ReadString(mol, entry[1].encode('utf8'))
xyz_geom = obConversion.WriteString(mol)
# if n==0:
#     print entry[0].encode('utf8')
#     print xyz_geom
rmg_mol = Molecule().fromSMILES(entry[0].encode('utf8').split('_')[0]) # for intra-H!
mult = rmg_mol.multiplicity
xyz_geom = "0 {0}\n".format(mult) + '\n'.join(xyz_geom.split('\n')[2:])
# if n==0:
#     print obConversion.WriteString(mol)
#     print rmg_mol.toAdjacencyList()
#     print xyz_geom

# Check for folder
reaction_folder = os.path.join("/home/slakman.b/Gaussian/SMD/", family, entry[0].encode('utf8'))
if not os.path.exists(reaction_folder):
    os.makedirs(reaction_folder)

# Write input file
options = "%mem=1GB\n%nprocshared=10"
keywords = "# m062x/6-311+G(2df,2p) scrf(smd, solvent=" + solvent +") int=ultrafine freq nosymm"
title = entry[0].encode('utf8')

input_file_path = os.path.join(reaction_folder, "ts_" + solvent)
input_file_ext = ".gjf"
input_file = open(input_file_path + input_file_ext, 'w')
input_file.write(options + "\n" + keywords + "\n\n" + title + "\n\n" + xyz_geom + "\n")
input_file.close()

# submits the input file to Gaussian
#import ipdb; ipdb.set_trace()
process = Popen(["/shared/apps/gaussian/G09_LINUX_LINDA/g09/g09", input_file_path + ".gjf", input_file_path + ".log"])
process.communicate() # necessary to wait for executable termination!
