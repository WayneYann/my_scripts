"""
Given a species dictionary and chemkin file, will modify the reaction barriers
for liquid phase based on the group additive scheme in RMG-Py, and will return
the chemkin file back with the barriers changed. This is for H-Abstraction for
now.
"""
import sys
import os
from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.solvation import DatabaseError, SoluteData, SolvationDatabase, SolvationKinetics
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.base import Database
#########################################################
database = Database()
species_dict = database.getSpecies('/home/slakman.b/Code/mech/MH_mechs1/prelim_mech_V3/RMG_Dictionary.txt')

reaction_list = []
family_list = ['H_Abstraction', 'intra_H_migration']

with open('/home/slakman.b/Code/mech/MH_mechs1/prelim_mech_V3/chem.inp', 'r') as mech_file:
    for line in mech_file:
        if line.strip().startswith('REACTIONS'): break
    for line in mech_file:
        if line.strip().startswith('!'): break
        if 'H_Abstraction' in line or 'intra_H_migration' in line: reaction_list.append(line.strip())

rmg_database = RMGDatabase()
rmg_database.load(settings['database.directory'], kineticsFamilies=family_list, reactionLibraries=[])
kinetics_database = rmg_database.kinetics
families = [kinetics_database.families[f] for f in family_list]

barrier_database = SolvationKinetics()
barrier_database.family = families[0]
barrier_database.load(os.path.join(settings['database.directory'], 'kinetics', 'families', "H_Abstraction"), None, None)

barrier_database2 = SolvationKinetics()
barrier_database2.family = families[1]
barrier_database2.load(os.path.join(settings['database.directory'], 'kinetics', 'families', "intra_H_migration"), None, None)

delEa_list = []
for rxn in reaction_list:
        if 'H_Abstraction' in rxn:
            bd = barrier_database
            reactants, products = rxn.split('=')[0], rxn.split('=')[1]
            r1_ID, r2_ID = reactants.split('+')
            p1_ID, p2_ID = products.split('+')[0], products.split('+')[1]
            p2_ID = p2_ID.split(')')[0] + ')'
            rSpecies1 = species_dict[r1_ID]
            rSpecies2 = species_dict[r2_ID]
            pSpecies1 = species_dict[p1_ID]
            pSpecies2 = species_dict[p2_ID]
            rSpecies1.generateResonanceIsomers()
            rSpecies2.generateResonanceIsomers()
            pSpecies1.generateResonanceIsomers()
            pSpecies2.generateResonanceIsomers()
            testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
            reactionList = []
            for moleculeA in rSpecies1.molecule:
                for moleculeB in rSpecies2.molecule:
        	    tempList = rmg_database.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=['H_Abstraction'])
                    for rxn0 in tempList:
                        reactionList.append(rxn0)
        else:
            bd = barrier_database2
            reactants, products = rxn.split('=')[0], rxn.split('=')[1]
            r1_ID = reactants
            p1_ID = products.split(')')[0] + ')'
            rSpecies1 = species_dict[r1_ID]
            pSpecies1 = species_dict[p1_ID]
            rSpecies1.generateResonanceIsomers()
            pSpecies1.generateResonanceIsomers()
            testReaction = Reaction(reactants=[rSpecies1], products=[pSpecies1], reversible=True)
            reactionList = []
            for moleculeA in rSpecies1.molecule:
                tempList = rmg_database.kinetics.generateReactionsFromFamilies([moleculeA], [], only_families=['intra_H_migration'])
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

        if not gotOne:
            continue

        barrierCorrection = bd.estimateBarrierCorrection(reaction)
        delEa_list.append(barrierCorrection.correction.value_si) # Value in J/mol

new_mech_file = open('/home/slakman.b/Code/mech/MH_mechs1/prelim_mech_V3/chem_modified.inp', 'a+')
num = 0
with open('/home/slakman.b/Code/mech/MH_mechs1/prelim_mech_V3/chem.inp', 'r') as mech_file:
    # write header and species list
    for line in mech_file:
        new_mech_file.write(line)
        if line.strip().startswith('REACTIONS'): break
    # write reactions, including the H_Abs one with corrected Ea
    for line in mech_file:
        if 'H_Abstraction' in line or 'intra_H_migration' in line:
           # get the correction, apply to appropriate place in line
           corr = delEa_list[num] / 4184 # kcal/mol
           Ea = float(line[72:77]) + corr
           Ea_string = str(Ea)
           # Make sure this string is 5 characters
           if len(Ea_string) < 5:
               i=0
               while i < 5-len(Ea_string):
                   Ea = Ea + " "
                   i += 1
           new_mech_file.write(line[:72] + str(Ea) + line[77:])
           num += 1
        else:
           new_mech_file.write(line)
           if line.strip().startswith('!'): break
    # write remainder of file
    for line in mech_file:
        new_mech_file.write(line)
new_mech_file.close()
