"""
Given two chemkin files and species dictionaries, will replace rates in the second
mechanism with the rate in the first.
"""
import sys
import os
from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.base import Database
#########################################################
database = Database()
species_dict_V3 = database.getSpecies('/home/slakman.b/Code/mech/MH_mechs1/prelim_mech_V3/RMG_Dictionary.txt')
species_dict_V4 = database.getSpecies('/home/slakman.b/Code/mech/MH_mechs1/final_mech_V4/RMG_Dictionary.txt')

reaction_list_V3 = []
reaction_list_V4 = []
family_list = ['H_Abstraction', 'intra_H_migration']

with open('/home/slakman.b/Code/mech/MH_mechs1/prelim_mech_V3/chem.inp', 'r') as mech_file:
    for line in mech_file:
        if line.strip().startswith('REACTIONS'): break
    for line in mech_file:
        if line.strip().startswith('!'): break
        if 'H_Abstraction' in line or 'intra_H_migration' in line: reaction_list_V3.append(line.strip())

with open('/home/slakman.b/Code/mech/MH_mechs1/final_mech_V4/chem.inp', 'r') as mech_file2:
    for line in mech_file2:
        if line.strip().startswith('REACTIONS'): break
    for line in mech_file2:
        if line.strip().startswith('!'): break
        if line.strip().startswith('END'): break
        reaction_list_V4.append(line.strip())

print reaction_list_V3
print "\n"
print reaction_list_V4

for i, rxn in enumerate(reaction_list_V4):
     h_abs = False
     if '+' in rxn.split('=')[0] and '+' in rxn.split('=')[1].split()[0]:
         h_abs = True # maybe, but it's at least bimolecular in both directions
         reactants, products = rxn.split('=')[0], rxn.split('=')[1]
         r1_ID, r2_ID = reactants.split('+')
         p1_ID, p2_ID = products.split('+')[0], products.split('+')[1]
         p2_ID = p2_ID.split(')')[0] + ')'
         rSpecies1 = species_dict_V4[r1_ID]
         rSpecies2 = species_dict_V4[r2_ID]
         pSpecies1 = species_dict_V4[p1_ID]
         pSpecies2 = species_dict_V4[p2_ID]
         rSpecies1.generateResonanceIsomers()
         rSpecies2.generateResonanceIsomers()
         pSpecies1.generateResonanceIsomers()
         pSpecies2.generateResonanceIsomers()
         testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
         # reactionList = []
 #             for moleculeA in rSpecies1.molecule:
 #                 for moleculeB in rSpecies2.molecule:
 #         	    tempList = rmg_database.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=['H_Abstraction'])
 #                     for rxn0 in tempList:
 #                         reactionList.append(rxn0)
     elif not '+' in rxn.split('=')[0] and not '+' in rxn.split('=')[1].split()[0]: # it's at least unimolecular in both directions
         reactants, products = rxn.split('=')[0], rxn.split('=')[1]
         r1_ID = reactants
         p1_ID = products.split(')')[0] + ')'
         rSpecies1 = species_dict_V4[r1_ID]
         pSpecies1 = species_dict_V4[p1_ID]
         rSpecies1.generateResonanceIsomers()
         pSpecies1.generateResonanceIsomers()
         testReaction = Reaction(reactants=[rSpecies1], products=[pSpecies1], reversible=True)
         # reactionList = []
 #             for moleculeA in rSpecies1.molecule:
 #                 tempList = rmg_database.kinetics.generateReactionsFromFamilies([moleculeA], [], only_families=['intra_H_migration'])
 #                 for rxn0 in tempList:
 #                     reactionList.append(rxn0)

     else:
         continue

     for r in reaction_list_V3:
         if h_abs:
             if 'H_Abstraction' in r:
                 reactants, products = r.split('=')[0], r.split('=')[1]
                 r1_ID, r2_ID = reactants.split('+')
                 p1_ID, p2_ID = products.split('+')[0], products.split('+')[1]
                 p2_ID = p2_ID.split(')')[0] + ')'
                 rSpecies1 = species_dict_V3[r1_ID]
                 rSpecies2 = species_dict_V3[r2_ID]
                 pSpecies1 = species_dict_V3[p1_ID]
                 pSpecies2 = species_dict_V3[p2_ID]
                 rSpecies1.generateResonanceIsomers()
                 rSpecies2.generateResonanceIsomers()
                 pSpecies1.generateResonanceIsomers()
                 pSpecies2.generateResonanceIsomers()
                 V3_reaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
                 if testReaction.isIsomorphic(V3_reaction):
                     reaction_list_V4[i] = r
                     break
         else:
             if 'intra_H_migration' in r:
                 reactants, products = r.split('=')[0], r.split('=')[1]
                 r1_ID = reactants
                 p1_ID = products.split(')')[0] + ')'
                 rSpecies1 = species_dict_V3[r1_ID]
                 pSpecies1 = species_dict_V3[p1_ID]
                 rSpecies1.generateResonanceIsomers()
                 pSpecies1.generateResonanceIsomers()
                 V3_reaction = Reaction(reactants=[rSpecies1], products=[pSpecies1], reversible=True)
                 if testReaction.isIsomorphic(V3_reaction):
                     reaction_list_V4[i] = r
                     break

new_mech_file = open('/home/slakman.b/Code/mech/MH_mechs1/final_mech_V4/chem_modified.inp', 'a+')
num = 0
with open('/home/slakman.b/Code/mech/MH_mechs1/final_mech_V4/chem.inp', 'r') as mech_file:
    # write header and species list
    for line in mech_file:
        new_mech_file.write(line)
        if line.strip().startswith('REACTIONS'): break
    for line in mech_file:
        new_mech_file.write(reaction_list_V4[num])
        num += 1
        if line.strip().startswith('!'): break
    # write remainder of file
    for line in mech_file:
        new_mech_file.write(line)
new_mech_file.close()
