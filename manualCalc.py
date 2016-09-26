from rmgpy.qm.reaction import QMReaction
from rmgpy.qm.main import QMSettings
from rmgpy.qm.gaussian import GaussianTSM062X
from rmgpy.reaction import Reaction
from rmgpy.molecule import Molecule, Atom
from rmgpy.species import Species
from rmgpy.data.rmg import RMGDatabase
import os

database = RMGDatabase()
database.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')),
kineticsFamilies=['Silylene_Insertion'], seedMechanisms=[], solvation=False)
tsDatabase = database.kinetics.families['Silylene_Insertion'].transitionStates

def fixLonePairMolecule(mol_string):
    """
    Fix an incorrect SMILES parsing for molecule that should have a lone pair
    We fix the first atom we encounter. Not sure what to do if there's more
    than one
    """
    if "(S)" in mol_string:
        smiles = mol_string.split("(S")[0]
        rmg_mol = Molecule().fromSMILES(smiles)
        for atom in rmg_mol.atoms:
            if atom.radicalElectrons >= 2:
                 atom.radicalElectrons -= 2
                 atom.lonePairs += 1
                 rmg_mol.update()
                 break
        return rmg_mol
    else: return Molecule().fromSMILES(mol_string)

rxnString = os.getcwd().split('/')[-2]
r1 = fixLonePairMolecule(rxnString.split('_')[0].split('+')[0])
r2 = fixLonePairMolecule(rxnString.split('_')[0].split('+')[1])
p1 = fixLonePairMolecule(rxnString.split('_')[1])#.split('+')[0])
#p2 = fixLonePairMolecule(rxnString.split('_')[1].split('+')[1])
rSpecies1 = Species(molecule=[r1])
rSpecies2 = Species(molecule=[r2])
pSpecies1 = Species(molecule=[p1])
#pSpecies2 = Species(molecule=[p2])
rSpecies1.generateResonanceIsomers()
rSpecies2.generateResonanceIsomers()
pSpecies1.generateResonanceIsomers()
#pSpecies2.generateResonanceIsomers()

testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1], reversible=True)
reactionList = []
for moleculeA in rSpecies1.molecule:
    for moleculeB in rSpecies2.molecule:
        tempList = database.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=['Silylene_Insertion'])
        for rxn0 in tempList:
            reactionList.append(rxn0)

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

qmSettings = QMSettings(
software='gaussian',
method='m062x',
fileStore=os.getcwd(),
scratchDirectory='/gss_gpfs_scratch/slakman.b/QMscratch'
)

qmReac = GaussianTSM062X(reaction, qmSettings, tsDatabase)
qmReac.generateTSGeometryDirectGuess()
