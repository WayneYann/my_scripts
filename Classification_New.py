
# coding: utf-8

# In[1]:

# %load classification.py
import os
import sys

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, KineticsDatabase
from rmgpy.data.rmg import RMGDatabase

import pandas as pd


# In[2]:

def get_features(reaction):
    """
    Given a rmgpy.reaction.Reaction() 'reaction', calculate the value of 
    the features for this reaction and return as a dictionary
    
    Currently specific to Silylene Insertion reactions
    """
    features = ['H2', 'Sirad_H', 'Sirad_Si_H', 'Sid_H', 'Sid_Si_H', 'sil_H2', 'sil_rad', 'sil_d', 'Success Code']
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


# In[22]:

def get_success(reaction):
    """
    Given a rmgpy.reaction.Reaction() 'reaction', get whether the TS calculation succeeded or how it
    failed (success, false positive, failure1, failure2)
    
    Currently, will only work if reaction is bimolecular->unimolecular (2 reactants->1)
    """
    scratch_dir = os.getenv('SCRATCH') 
    if scratch_dir is None:
        scratch_dir = os.getenv('RMGpy')
    
    success = None
    family = reaction.family
    try:
        my_folder = find_reaction_folder(reaction)
    except TypeError:
        print "Couldn't find any folders for the reaction: {0}+{1}->{2}".format(
            reaction.reactant[0].molecule[0].toSMILES(), reaction.reactant[1].molecule[0].toSMILES,
            reaction.product[0].molecule[0].toSMILES())
    
    if my_folder:
        if os.path.exists(os.path.join('false_positives', my_folder)):
            success = 'FP' # false positive
        else:
            if os.path.exists(os.path.join(scratch_dir, 'QMfiles/Reactions', family, my_folder, 'm062x', 'error.txt')):
                with open(os.path.join(scratch_dir, 'QMfiles/Reactions', family, my_folder, 'm062x', 'error.txt')) as error_file:
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
            elif os.path.exists(os.path.join(scratch_dir, 'QMfiles/Reactions', family, my_folder, 'm062x', 'optDists.txt')):
                success = 'S'
            else:
                print "Not sure about reaction {0}, has contents {1}\n".format(my_folder, 
                            ' , '.join(os.listdir(os.path.join(scratch_dir, 'QMfiles/Reactions', 
                                                               family, my_folder, 'm062x'))))

    return success


# In[16]:

def find_reaction_folder(reaction, family=None):
    """
    Given rmgpy.reaction.Reaction() 'reaction', return the folder containing the TS calculation results as a string
    """
    scratch_dir = os.getenv('SCRATCH') 
    if scratch_dir is None:
        scratch_dir = os.getenv('RMGpy')
        
    if family is None:
        family = reaction.family
        
    r1_smiles = reaction.reactants[0].molecule[0].toSMILES()
    r2_smiles = reaction.reactants[1].molecule[0].toSMILES()
    p_smiles = reaction.products[0].molecule[0].toSMILES()
    possible_folders = ["{0}+{1}_{2}".format(r1_smiles, r2_smiles, p_smiles), "{0}+{1}_{2}".format(r2_smiles, r1_smiles, p_smiles)]
            
    my_folder = None
    for folder in possible_folders:
        if os.path.exists(os.path.join(scratch_dir, 'QMfiles/Reactions', family, folder)):
            my_folder = folder
            break
    return my_folder


# In[5]:

def get_rmg_reaction(database, species_dict, family, mech_line):
    """
    This function takes in:
    
    a loaded rmgpy.data.rmg.RMGDatabase() 'database'
    RMG reaction family 'family' as a string
    a reaction line from a chemkin file 'mech_line' as a string
    species dictionary 'species_dict' associated with the mechanism
    
    It returns the associated RMG reaction if it's possible 
    to make one. If not, return 'None'
    """
    rxnFormula, A, n, Ea = mech_line.split()
    reactants, products = rxnFormula.split('=')
    if family in ['H_Abstraction', 'Disproportionation', 'Cl-Abstraction']:
        rSpecies1, rSpecies2 = [species_dict[j] for j in reactants.split('+')]
        pSpecies1, pSpecies2 = [species_dict[j] for j in products.split('+')]
        rSpecies1.generateResonanceIsomers()
        rSpecies2.generateResonanceIsomers()
        pSpecies1.generateResonanceIsomers()
        pSpecies2.generateResonanceIsomers()
        testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
        reactionList = []
        for moleculeA in rSpecies1.molecule:
            for moleculeB in rSpecies2.molecule:
                tempList = database.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[family])
                for rxn0 in tempList:
                    reactionList.append(rxn0)
    elif family in ['intra_H_migration']:
        rSpecies = species_dict[reactants]
        pSpecies = species_dict[products]
        rSpecies.generateResonanceIsomers()
        pSpecies.generateResonanceIsomers()
        testReaction = Reaction(reactants=[rSpecies], products=[pSpecies], reversible=True)
        reactionList = []
        for moleculeA in rSpecies.molecule:
            tempList = database.kinetics.generateReactionsFromFamilies([moleculeA], [], only_families=[family])
            for rxn0 in tempList:
                reactionList.append(rxn0)
    elif family in ['R_Addition_MultipleBond', 'Silylene_Insertion']:
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
        reactionList = []
        for moleculeA in rSpecies1.molecule:
            for moleculeB in rSpecies2.molecule:
                tempList = database.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[family])
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
                return reaction
    
    return None


# In[46]:

def load_mechanism(mech_file, dict_file):
    """
    Loads the RMG database and processes a chemical mechanism. 
    
    Returns loaded rmgpy.data.rmg.RMGDatabase(), list of 
    reactions and their families, and species dictionary
    """
    print 'Loading RMG Database ...'
    families = ['Silylene_Insertion', 'H_Abstraction']
    rmgDatabase = RMGDatabase()
    rmgDatabase.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')), 
                     kineticsFamilies=families, seedMechanisms=[], solvation=False)
    print 'Finished loading RMG Database ...'

    loadSpecies = rmgDatabase.kinetics.families[families[0]] # any family will do
    species_dict = loadSpecies.getSpecies(dict_file)

    file_object = open(mech_file, 'r')
    mechLines = file_object.readlines()

    rxnList = []
    gotit = []
    for k, line in enumerate(mechLines):
        if line.startswith('! Template reaction:'):
            for m in range(10):
                reaction = mechLines[k+m].split()[0]
                if not reaction.startswith('!'):
                    break
            if reaction not in gotit:
                gotit.append(reaction)
                rxnList.append((line.split(': ')[1], mechLines[k+m]))

    return rmgDatabase, rxnList, species_dict


# In[47]:

def get_family_data(family):
    """
    Given a reaction family 'family' as string, returns a pandas.DataFrame that has the features for each
    reaction in a mechanism, in the family, that had a TS calculation done and what the success code is. 
    
    Features only will be relevant for silylene insertion
    """
    data = []

    rmg_database, reaction_list, species_dict = load_mechanism(
        'sih4_mech/chem_edge_annotated.inp', 'sih4_mech/species_edge_dictionary.txt')

    for reaction_tuple in reaction_list:
        rxn_family, reaction_line = reaction_tuple
        if rxn_family == family:
            reaction = get_rmg_reaction(rmg_database, species_dict, rxn_family, reaction_line)

            if reaction is not None:
                features_for_reaction = get_features(reaction)
                features_for_reaction['Success Code'] = get_success(reaction)
                data.append(features_for_reaction)
            else:
                print "No RMG reactions found for reaction {0}".format(reaction_line.split()[0])

    return pd.DataFrame().from_dict(data)


# In[48]:

get_family_data('Silylene_Insertion').to_csv('data_Silylene_Insertion.csv')


# In[ ]:



