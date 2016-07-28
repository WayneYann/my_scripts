from rmgpy.qm.reaction import QMReaction
from rmgpy.qm.main import QMSettings
from rmgpy.reaction import Reaction
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.rmg import RMGDatabase
import os

database = RMGDatabase()
database.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')),
kineticsFamilies=['H_Abstraction'], seedMechanisms=[], solvation=False)

rxnString = os.getcwd().split('/')[-1]
r1 = Molecule().fromSMILES(rxnString.split('_')[0].split('+')[0])
r2 = Molecule().fromSMILES(rxnString.split('_')[0].split('+')[1])
p1 = Molecule().fromSMILES(rxnString.split('_')[1].split('+')[0])
p2 = Molecule().fromSMILES(rxnString.split('_')[1].split('+')[1])
rSpecies1 = Species(molecule=[r1])
rSpecies2 = Species(molecule=[r2])
pSpecies1 = Species(molecule=[p1])
pSpecies2 = Species(molecule=[p2])
rSpecies1.generateResonanceIsomers()
rSpecies2.generateResonanceIsomers()
pSpecies1.generateResonanceIsomers()
pSpecies2.generateResonanceIsomers()

testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
reactionList = []
for moleculeA in rSpecies1.molecule:
    for moleculeB in rSpecies2.molecule:
        tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=['H_Abstraction'])
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
fileStore='/gss_gpfs_scratch/slakman.b/QMfiles',
scratchDirectory='/gss_gpfs_scratch/slakman.b/QMscratch'
)

qmReac = GaussianTSM062X(reaction, qmSettings, database)
qmReac.generateTSGeometryDirectGuess()
