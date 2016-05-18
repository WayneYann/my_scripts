"""
Perform single point energy calculations in solvent on a set of gas phase
reactant geometries (gets list from database and geometry from existing gas
phase calculations)
"""
import sqlite3
import openbabel as ob
from subprocess import Popen
import os
from shutil import copy
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

geometries = conn.execute("SELECT uniqueID, r1 FROM (SELECT * from TSData WHERE rxnFamily=\'" + family + "\') WHERE method='m062x'") # returns cursor
geometries = list(geometries.fetchall()) # make into list
entry = geometries[i-1]

# Check for folder
reaction_folder = os.path.join("/home/slakman.b/Gaussian/SMD/", family, entry[0].encode('utf8'))
if not os.path.exists(reaction_folder):
	os.makedirs(reaction_folder)

err_message = None
# Find reactant gas-phase output file in scratch and copy to my reaction folder
if not os.path.exists(os.path.join(reaction_folder, entry[1] + ".log")):
    try:
    	gas_phase_output = os.path.join('/scratch/bhoorasingh.p/QMscratch/Species', entry[1], 'm062x', entry[1] + ".log")
    	copy(gas_phase_output, reaction_folder)
    except:
	err_message = "Gas phase output file not found for {0}".format(entry[1])

gas_phase_output = os.path.join(reaction_folder, entry[1] + ".log")

# Find ts gas-phase output file in scratch and copy to my reaction folder
try:
	gas_phase_output_ts = os.path.join('/scratch/bhoorasingh.p/QMfiles/Reactions', family, entry[0], 'm062x', entry[0] + ".log")
	copy(gas_phase_output_ts, os.path.join(reaction_folder, "ts.log"))
except:
	err_message = "Gas phase output file not found for {0}".format(entry[0])

if err_message is None:
	# Extract geometry from gas-phase output
	xyz_geom = ""
	with open(gas_phase_output, 'r') as gpf:
            lines = gpf.read().split('\n')
            input_geoms = [i for i,x in enumerate(lines) if 'Input orientation' in x]
	    geom_start = input_geoms[-1] + 5
	    geom_end = [lines[geom_start:].index(l) for l in lines[geom_start:] if "---" in l][0] + geom_start

	for atom_line in lines[geom_start:geom_end]:
	    atomic_number = int(atom_line.split()[1])
	    coordinates = atom_line.split()[-3:]
	    if atomic_number == 6: element = "C"
	    elif atomic_number == 1: element = "H"
	    elif atomic_number == 8: element = "O"
            else: element = "x"
            xyz_geom += "{0}\t{1} {2} {3}\n".format(element, coordinates[0], coordinates[1], coordinates[2])

	# multiplicity
	rmg_mol = Molecule().fromSMILES(entry[1].encode('utf8')) # for intra-H!
	mult = rmg_mol.multiplicity
	xyz_geom = "0 {0}\n".format(mult) + xyz_geom

	# Write solvation input file
	options = "%mem=1GB\n%nprocshared=10"
	keywords = "# m062x/6-311+G(2df,2p) scrf(smd, solvent=" + solvent +") int=ultrafine freq nosymm"
	title = entry[1].encode('utf8')

	input_file_path = os.path.join(reaction_folder, entry[1] + "_" + solvent)
	input_file_ext = ".gjf"
	input_file = open(input_file_path + input_file_ext, 'w')
	input_file.write(options + "\n" + keywords + "\n\n" + title + "\n\n" + xyz_geom + "\n")
	input_file.close()

	# submits the input file to Gaussian
	#import ipdb; ipdb.set_trace()
	process = Popen(["/shared/apps/gaussian/G09_LINUX_LINDA/g09/g09", input_file_path + ".gjf", input_file_path + ".log"])
	process.communicate() # necessary to wait for executable termination!

else:
	err_file = open(os.path.join(reaction_folder, "err.txt"), 'w')
	err_file.write(err_message)
	err_file.close()
