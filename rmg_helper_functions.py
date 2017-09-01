
# coding: utf-8

# In[1]:

import os
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.data.rmg import RMGDatabase


# In[2]:

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


# In[3]:

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


# In[4]:

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

