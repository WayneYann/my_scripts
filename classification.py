import os
import sys

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, KineticsDatabase
from rmgpy.data.rmg import RMGDatabase
import pandas as pd

def get_features(features, reaction):
	feature_dict = dict.fromkeys(features[:-1], False)
	if '[H][H]' in [r.molecule[0].toSMILES() for r in reaction.reactants]:
		feature_dict.update(dict.fromkeys(['H2'], True))
	else:
		for r in reaction.reactants:
			try:
				si_h = r.molecule[0].getLabeledAtom('*1')
			except:
				pass
		if si_h.radicalElectrons > 0:
			feature_dict.update(dict.fromkeys(['Sirad_H'], True))
		for atom, bond in si_h.bonds.items():
			if atom.element.symbol == 'Si':
				if atom.radicalElectrons > 0:
					feature_dict.update(dict.fromkeys(['Sirad_Si_H'], True))
				for atom2, bond2 in atom.bonds.items():
					if atom2 != si_h and bond2.order == 'D':
						feature_dict.update(dict.fromkeys(['Sid_Si_H'], True))
			if bond.order == 'D':
				feature_dict.update(dict.fromkeys(['Sid_H'], True))
	for p in reaction.products:
		try:
			sil_h = p.molecule[0].getLabeledAtom('*3')
		except:
			pass
	neighbor_atoms = [neighbor[0] for neighbor in sil_h.bonds.items()]
	if len(neighbor_atoms) == 2 and all([neighbor.atom.symbol == 'H' for neighbor in neighbor_atoms]):
		feature_dict.update(dict.fromkeys(['sil_H2'], True))
	else:
		if sil_h.radicalElectrons > 0:
			feature_dict.update(dict.fromkeys(['sil_rad'], True))
		for atom, bond in sil_h.bonds.items():
			if bond.order == 'D':
				feature_dict.update(dict.fromkeys(['sil_d'], True))

	return feature_dict

def get_success(family, possibleFolders):
	success = None
	my_folder = None
	for folder in possibleFolders:
		if os.path.exists(os.path.join('/gss_gpfs_scratch/slakman.b/QMfiles/Reactions', family, folder)):
			my_folder = folder
			break
	if my_folder is None:
		print "Couldn't find any of the following: {0}".format(','.join(possibleFolders))
		return

	if os.path.exists(os.path.join('false_positives', my_folder)):
		success = 'FP' # false positive
	else:
		if os.path.exists(os.path.join('/gss_gpfs_scratch/slakman.b/QMfiles/Reactions', family, my_folder, 'm062x', 'error.txt')):
			with open(os.path.join('/gss_gpfs_scratch/slakman.b/QMfiles/Reactions', family, my_folder, 'm062x', 'error.txt')) as error_file:
				for line in error_file: # There should only be one
					if line.strip().startswith('Success'):
						success = 'S'
						break
					if line.strip().startswith('TS not converged'):
						success = 'F1' # TS not converged failure
						break
					if line.strip().startswith('IRC failed'):
						success = 'F2' # IRC failure
						break
			if success is None:
				print "Strange error {0} in reaction {1}".format(line, my_folder)
		elif os.path.exists(os.path.join('/gss_gpfs_scratch/slakman.b/QMfiles/Reactions', family, my_folder, 'm062x', 'optDists.txt')):
			success = 'S'
		else:
			print "Not sure about reaction {0}!".format(my_folder)

	return success

def get_family_data(family):

	data = []

	print 'Loading RMG Database ...'
	rmgDatabase = RMGDatabase()
	rmgDatabase.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')), kineticsFamilies=[family], seedMechanisms=[], solvation=False)
	print 'Finished loading RMG Database ...'

	# Doesn't matter which family, it loads the entire species dict
	loadSpecies = rmgDatabase.kinetics.families[family]
	species_dict_file = 'sih4_mech/species_edge_dictionary.txt'
	species_dict = loadSpecies.getSpecies(species_dict_file)

	file_object = open('sih4_mech/chem_edge_annotated.inp', 'r')
	mechLines = file_object.readlines()

	rxnList = []
	gotit = []
	for k, line in enumerate(mechLines):
		if line.startswith('! Template reaction: {0}'.format(family)):
			for m in range(10):
				reaction = mechLines[k+m].split()[0]
				if not reaction.startswith('!'):
					break
			if reaction not in gotit:
				gotit.append(reaction)
				rxnList.append((family, mechLines[k+m]))

	features = ['H2', 'Sirad_H', 'Sirad_Si_H', 'Sid_H', 'Sid_Si_H', 'sil_H2', 'sil_rad', 'sil_d', 'Success Code']

	for reactionTuple in rxnList:
		rxnFamily, reactionLine = reactionTuple
		rxnFormula, A, n, Ea = reactionLine.split()
		reactants, products = rxnFormula.split('=')
		if rxnFamily in ['R_Addition_MultipleBond', 'Silylene_Insertion']:
			if '(+M)' in reactants:
				reactants = reactants.split('(+M)')[0]
				products = products.split('(+M)')[0]
			if len(reactants.split('+'))==2:
				rSpecies1, rSpecies2 = [species_dict[j] for j in reactants.split('+')]
				pSpecies = species_dict[products]
			else:
				rSpecies1, rSpecies2 = [species_dict[j] for j in products.split('+')]
				pSpecies = species_dict[reactants]
			rSpecies1.generateResonanceIsomers()
			rSpecies2.generateResonanceIsomers()
			pSpecies.generateResonanceIsomers()
			testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies], reversible=False)

			# Figure out whether we ran a calculation for this reaction or not. Based on SMILES.
			r1_smiles = rSpecies1.molecule[0].toSMILES()
			r2_smiles = rSpecies2.molecule[0].toSMILES()
			p_smiles = pSpecies.molecule[0].toSMILES()
			possibleFolders = ["{0}+{1}_{2}".format(r1_smiles, r2_smiles, p_smiles), "{0}+{1}_{2}".format(r2_smiles, r1_smiles, p_smiles)]
			exists = [os.path.exists(os.path.join('/gss_gpfs_scratch/slakman.b/QMfiles/Reactions', family, folder)) for folder in possibleFolders]
			if exists: # Make the RMG reaction
				reactionList = []
				for moleculeA in rSpecies1.molecule:
					for moleculeB in rSpecies2.molecule:
						tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[family])
						for rxn0 in tempList:
							reactionList.append(rxn0)
					gotOne=False
					for reaction in reactionList:
						# Check if any of the RMG proposed reactions matches the reaction in the mechanism
						if reaction.isIsomorphic(testReaction):
							# Now add the labeled atoms to the Molecule, and check all labels were added
							atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
							atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

							for reactant in reaction.reactants:
								reactant = reactant.molecule[0]
								reactant.clearLabeledAtoms()
								for atom in reactant.atoms:
									for atomLabel in reaction.labeledAtoms:
										if atom==atomLabel[1]:
											atom.label = atomLabel[0]
											atLblsR[atomLabel[0]] = True
							for product in reaction.products:
								product = product.molecule[0]
								product.clearLabeledAtoms()
								for atom in product.atoms:
									for atomLabel in reaction.labeledAtoms:
										if atom==atomLabel[1]:
											atom.label = atomLabel[0]
											atLblsP[atomLabel[0]] = True
							if all( atLblsR.values() ) and all( atLblsP.values() ):
								gotOne=True
								break

				if gotOne:
					features_for_reaction = get_features(features, reaction)
					features_for_reaction['Success Code'] = get_success(family, possibleFolders)
					data.append(features_for_reaction)

	return pd.DataFrame().from_dict(data)

get_family_data('Silylene_Insertion').to_csv('data_Silylene_Insertion.csv')
