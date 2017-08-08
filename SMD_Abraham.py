#!/usr/bin/env python
# encoding: utf-8

# A script that uses quantum chemistry outputs for kinetic solvent effects and
#creates a .csv file that can be used for later processing, like making a
# decision tree regressor. This particular script will save the Abraham
# parameters, and possibly some solvent parameters

import os.path
import argparse
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.data.solvation import SolvationDatabase
import pandas as pd

parser = argparse.ArgumentParser(description="""
Given solvent, will get the delEa values from the output files.
""")
parser.add_argument("-s", "--solvent", help="Name of solvent, i.e. \'water\'")
parser.add_argument("-f", "--family", help="Name of family, i.e. \'intra_H_migration\'")
args= parser.parse_args()

rmgDatabase = RMGDatabase()
rmgDatabase.load(settings['database.directory'], kineticsFamilies=[args.family], reactionLibraries=[], solvation=True)

db_folder = '/Users/belinda/Code/RMG-database/input/kinetics/families'

data_dir = os.path.join('/Users/belinda/Gaussian/SMD/', args.family)
data_file = 'output_' + args.solvent + '.txt'

if args.family == "H_Abstraction":
    data = [['r1S','r1B','r1E','r1L','r1A',
    'r2S','r2B','r2E','r2L','r2A',
    'p1S','p1B','p1E','p1L','p1A',
    'p2S','p2B','p2E','p2L','p2A','barrier_correction']]

# Since this is just for one solvent, leave out any solvent parameters, but this
# could be a consideration in the future.

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

    f = open(os.path.join(subdir,data_file), 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
       if line.startswith("Difference"):
          text, barrierCorrection = line.split("= ")
          if abs(float(barrierCorrection)) > 50.0:
              print reactionString

    if args.family == "H_Abstraction":
        r1_data = rmgDatabase.solvation.getSoluteData(Species(molecule=my_reaction.reactants[0].generateResonanceIsomers()))
        r2_data = rmgDatabase.solvation.getSoluteData(Species(molecule=my_reaction.reactants[1].generateResonanceIsomers()))
        p1_data = rmgDatabase.solvation.getSoluteData(Species(molecule=my_reaction.products[0].generateResonanceIsomers()))
        p2_data = rmgDatabase.solvation.getSoluteData(Species(molecule=my_reaction.products[1].generateResonanceIsomers()))
        this_data = [r1_data.S, r1_data.B, r1_data.E, r1_data.L, r1_data.A,
        r2_data.S, r2_data.B, r2_data.E, r2_data.L, r2_data.A,
        p1_data.S, p1_data.B, p1_data.E, p1_data.L, p1_data.A,
        p2_data.S, p2_data.B, p2_data.E, p2_data.L, p2_data.A, barrierCorrection]
        data.append(this_data)
    if args.family == "intra_H_migration":
        #other stuff
        data.append([])
pd.DataFrame(data[1:], columns=data[0]).to_csv('data_Abraham.csv')
