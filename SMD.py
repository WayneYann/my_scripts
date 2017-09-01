#!/usr/bin/env python
# encoding: utf-8

# A script for gathering the SMD barrier correction values from Gaussian jobs and putting into the RMG database.
# Also, creates a .csv file that can be used for later processing, like making a decision tree regressor.

import os.path
import argparse
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
import pandas as pd

parser = argparse.ArgumentParser(description="""
Given solvent, will get the delEa values from the output files.
""")
parser.add_argument("-s", "--solvent", help="Name of solvent, i.e. \'water\'")
parser.add_argument("-f", "--family", help="Name of family, i.e. \'intra_H_migration\'")
args= parser.parse_args()

rmgDatabase = RMGDatabase()
rmgDatabase.load(settings['database.directory'], kineticsFamilies=[args.family], reactionLibraries=[])

db_folder = '/Users/belinda/Code/RMG-database/input/kinetics/families'

data_dir = os.path.join('/Users/belinda/Gaussian/SMD/', args.family)
data_file = 'output_' + args.solvent + '.txt'
if not os.path.exists(os.path.join(db_folder, args.family, args.solvent, 'training')):
    os.makedirs(os.path.join(db_folder, args.family, args.solvent, 'training'))
output_file = os.path.join(db_folder, args.family, args.solvent, 'training', 'solvationReactions.py')
output_dict = os.path.join(db_folder, args.family, args.solvent, 'training', 'solvationDictionary.txt')

out = open(output_file, 'w+')
out.write("""
#!/usr/bin/env python
# encoding: utf-8

name = \"""" + args.family + """\/training\"
shortDesc = u\"A list of solvation reactions used to train group additivity values\"
longDesc = u\"\"\"
\"\"\"
"""
)
out.close()

index = 1
unique = {}

if args.family == "H_Abstraction":
    data = [["X_rad_count",
"X_lp_count",
"X_element",
"X_bonds",
"Y_rad_count",
"Y_lp_count",
"Y_element",
"Y_bonds",
'X_CO',
 'X_CS',
 'X_Cb',
 'X_Cbf',
 'X_Cd',
 'X_Cdd',
 'X_Cs',
 'X_Ct',
 'X_H',
 'X_N3b',
 'X_N3d',
 'X_N3s',
 'X_N3t',
 'X_Oa',
 'X_Od',
 'X_Os',
 'X_Ot',
 'Y_CO',
 'Y_CS',
 'Y_Cb',
 'Y_Cbf',
 'Y_Cd',
 'Y_Cdd',
 'Y_Cs',
 'Y_Ct',
 'Y_H',
 'Y_N3b',
 'Y_N3d',
 'Y_N3s',
 'Y_N3t',
 'Y_Oa',
 'Y_Od',
 'Y_Os',
 'Y_Ot',
"barrier_correction"]]

if args.family == "intra_H_migration":
    data = [[]]

for subdir, dirs, files in os.walk(data_dir):
    reactionString = os.path.basename(os.path.normpath(subdir))
    #if reactionString in ('H_Abs','solv2'): continue
    if not data_file in files: continue
    try:
        reactants, products = reactionString.split("_")
    except:
        print("String {0} is wrong length".format(reactionString))
        raise
    reactants = reactants.split("+")
    products = products.split("+")

    # Create the real reaction, so we can get the labeled atoms for the dictionary
    if args.family in ['H_Abstraction', 'Disproportionation', 'Cl-Abstraction']:
    	rSpecies1, rSpecies2 = [Species(molecule=[Molecule().fromSMILES(r)]) for r in reactants]
    	pSpecies1, pSpecies2 = [Species(molecule=[Molecule().fromSMILES(p)]) for p in products]
    	rSpecies1.generateResonanceIsomers()
    	rSpecies2.generateResonanceIsomers()
    	pSpecies1.generateResonanceIsomers()
    	pSpecies2.generateResonanceIsomers()
    	testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
    	reactionList = []
    	for moleculeA in rSpecies1.molecule:
    		for moleculeB in rSpecies2.molecule:
    			tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[args.family])
    			for rxn0 in tempList:
    				reactionList.append(rxn0)
    elif args.family in ['intra_H_migration']:
    	rSpecies = Species(molecule=[Molecule().fromSMILES(reactants[0])])
    	pSpecies = Species(molecule=[Molecule().fromSMILES(products[0])])
    	rSpecies.generateResonanceIsomers()
    	pSpecies.generateResonanceIsomers()
    	testReaction = Reaction(reactants=[rSpecies], products=[pSpecies], reversible=True)
    	reactionList = []
    	for moleculeA in rSpecies.molecule:
    		tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA], [], only_families=[args.family])
    		for rxn0 in tempList:
    			reactionList.append(rxn0)
    elif args.family in ['R_Addition_MultipleBond', 'Silylene_Insertion']:
        rSpecies1, rSpecies2 = [Species(molecule=Molecule().fromSMILES(r)) for r in reactants]
    	pSpecies = [Species(molecule=Molecule().fromSMILES(p)) for p in products]
    	rSpecies1.generateResonanceIsomers()
    	rSpecies2.generateResonanceIsomers()
    	pSpecies.generateResonanceIsomers()
    	testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies], reversible=False)
    	reactionList = []
    	for moleculeA in rSpecies1.molecule:
    		for moleculeB in rSpecies2.molecule:
    			tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[args.family])
    			for rxn0 in tempList:
    				reactionList.append(rxn0)
    gotOne=False
    for reaction in reactionList:
    	# Check if any of the RMG proposed reactions matches the reaction in the mechanism
    	if testReaction.isIsomorphic(reaction):
    		# Now add the labeled atoms to the Molecule, and check all labels were added
    		atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
    		atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

    		for reactant in reaction.reactants:
    			#reactant = reactant.molecule[0]
    			reactant.clearLabeledAtoms()
    			for atom in reactant.atoms:
    				for atomLabel in reaction.labeledAtoms:
    					if atom==atomLabel[1]:
    						atom.label = atomLabel[0]
    						atLblsR[atomLabel[0]] = True
    		for product in reaction.products:
    			#product = product.molecule[0]
    			product.clearLabeledAtoms()
    			for atom in product.atoms:
    				for atomLabel in reaction.labeledAtoms:
    					if atom==atomLabel[1]:
    						atom.label = atomLabel[0]
    						atLblsP[atomLabel[0]] = True
    		if all( atLblsR.values() ) and all( atLblsP.values() ):
    			gotOne=True
                my_reaction = reaction
                break

    if not gotOne: continue

    for i, r in enumerate(reactants):
        if r in unique:
            unique[r] += 1
            reactants[i] = "{0}-{1}".format(reactants[i], unique[r])
        else: unique.update({r: 1})
    for j, p in enumerate(products):
        if p in unique:
            unique[p] += 1
            products[j] = "{0}-{1}".format(products[j], unique[p])
        else: unique.update({p: 1})

    f = open(os.path.join(subdir,data_file), 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
       if line.startswith("Difference"):
          text, barrierCorrection = line.split("= ")
          if abs(float(barrierCorrection)) > 50.0:
              print reactionString

    out = open(output_file, 'a')
    out.write('entry( \n')
    out.write('    index = ' + str(index) + ', \n')
    out.write('    reactants = [\'' + "\' , \'".join(reactants) + '\'], \n')
    out.write('    products = [\'' + "\' , \'".join(products) + '\'], \n')
    out.write('    solvent = \'' + args.solvent + '\', \n')
    out.write('    correction = BarrierCorrection(correction = (' + barrierCorrection + ', \'kJ/mol\')), \n')
    out.write('    shortDesc = u"""MO6-2X/MG3S calculations in g09 with SMD solvation model""", \n')
    out.write('    longDesc = \n u\"\"\" \n \"\"\" \n )')
    out.write('\n')
    out.close()

    index += 1

    # Create the dictionary
    if os.path.exists(output_dict):
        d = open(output_dict, 'r')
    else:
        d = open(output_dict, 'a+')
    entries = d.read().splitlines()
    d.close()
    for ind, reac in enumerate(reactants):
        if not reac in entries:
           d = open(output_dict, 'a')
           d.write(reac + '\n')
           d.write(my_reaction.reactants[ind].toAdjacencyList() + '\n')
           d.close()
    for indP, prod in enumerate(products):
        if not prod in entries:
            d = open(output_dict, 'a')
            d.write(prod + '\n')
            d.write(my_reaction.products[indP].toAdjacencyList() + '\n')
            d.close()

    x_atTypes = {
    'Cs' : 0,
    'Cd' : 0,
    'Cdd': 0,
    'Ct' : 0,
    'CO' : 0,
    'Cb' : 0,
    'Cbf': 0,
    'CS' : 0,
    'H'  : 0,
    'Os' : 0,
    'Od' : 0,
    'Oa' : 0,
    'Ot' : 0,
    'N3s': 0,
    'N3d': 0,
    'N3t': 0,
    'N3b': 0 }

    y_atTypes = x_atTypes.copy()

    if args.family == "H_Abstraction":
        for m in my_reaction.reactants:
            try:
                x_atom = m.getLabeledAtom('*1')
                x_neighbors = [atom for atom in m.getLabeledAtom('*1').bonds]
            except:
                pass
            try:
                y_atom = m.getLabeledAtom('*3')
                y_neighbors = [atom for atom in m.getLabeledAtom('*3').bonds]
            except:
                pass
        x_rad = x_atom.radicalElectrons
        x_lp = x_atom.lonePairs
        x_atomType = x_atom.atomType.label
        x_element = x_atom.element.symbol
        if len(x_atom.bonds) > 0:
            x_bonds = "".join(sorted([x_atom.bonds[b].order for b in x_atom.bonds]))
        else:
            x_bonds = "-"
        for bond in x_atom.bonds:
            atom_type = bond.atomType.label
            x_atTypes[atom_type] += 1
        x_neighbors = [a[1] for a in sorted(x_atTypes.items())]

        y_rad = y_atom.radicalElectrons
        y_lp = y_atom.lonePairs
        y_atomType = y_atom.atomType.label
        y_element = y_atom.element.symbol
        if len(y_atom.bonds) > 0:
            y_bonds = "".join(sorted([y_atom.bonds[b].order for b in y_atom.bonds]))
        else:
            y_bonds = "-"
        for bond in y_atom.bonds:
            atom_type = bond.atomType.label
            y_atTypes[atom_type] += 1
        y_neighbors = [a[1] for a in sorted(y_atTypes.items())]

        #data.append([x_rad, x_lp, x_atomType, x_neighborC, x_neighborO, x_neighborH, x_neighborN, y_rad, y_lp, y_atomType, y_neighborC, y_neighborO, y_neighborH, y_neighborN, barrierCorrection])
        #data.append([x_rad, x_lp, x_element, x_bonds, x_neighborC, x_neighborO, x_neighborH, x_neighborN, y_rad, y_lp, y_element, y_bonds, y_neighborC, y_neighborO, y_neighborH, y_neighborN, barrierCorrection])
        this_data = [x_rad, x_lp, x_element, x_bonds, y_rad, y_lp, y_element, y_bonds]
        this_data.extend(x_neighbors)
        this_data.extend(y_neighbors)
        this_data.append(barrierCorrection)
        data.append(this_data)
    if args.family == "intra_H_migration":
        #other stuff
        data.append([])
pd.DataFrame(data[1:], columns=data[0]).to_csv('data_h_abs_water.csv')
