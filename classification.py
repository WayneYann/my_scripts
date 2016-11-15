
# coding: utf-8

# In[1]:

import os

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, KineticsDatabase

from rmg_helper_functions import find_reaction_folder, load_mechanism, get_rmg_reaction

import pandas as pd


# In[2]:

database, reaction_list, species_dict = load_mechanism(
        'sih4_mech/chem_edge_annotated.inp', 'sih4_mech/species_edge_dictionary.txt')


# In[3]:

print reaction_list[5]
rxn = get_rmg_reaction(database, species_dict, 'Silylene_Insertion', reaction_list[5][1])
print rxn.reactants, rxn.products
rxn


# In[4]:

def get_features(reaction):
    """
    Given a rmgpy.reaction.Reaction() 'reaction', calculate the value of 
    the features for this reaction and return as a dictionary
    
    Currently specific to Silylene Insertion reactions
    """
    features = ['reaction', 'H2', 'Sirad_H', 'Sirad_Si_H', 'Sid_H', 'Sid_Si_H', 'sil_H2', 'sil_rad', 'sil_d', 'Success Rate']
    feature_dict = dict.fromkeys(features[:-1], False)
    feature_dict.update(dict.fromkeys(['reaction'], 
                                      '+'.join([r.molecule[0].toSMILES() for r in reaction.reactants]) 
                                      + '_' + '+'.join([p.molecule[0].toSMILES() for p in reaction.products])))
    
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
    
    for r in reaction.reactants:
        try:
            sil_h = r.molecule[0].getLabeledAtom('*3')
        except:
            pass
    neighbor_atoms = [neighbor[0] for neighbor in sil_h.bonds.items()]
    if len(neighbor_atoms) == 2 and all([neighbor.symbol == 'H' for neighbor in neighbor_atoms]):
        feature_dict.update(dict.fromkeys(['sil_H2'], True))
    else:
        if sil_h.radicalElectrons > 0:
            feature_dict.update(dict.fromkeys(['sil_rad'], True))
        for atom, bond in sil_h.bonds.items():
            if bond.order == 'D':
                feature_dict.update(dict.fromkeys(['sil_d'], True))

    return feature_dict


# In[5]:

print get_features(rxn)


# In[6]:

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
        if os.path.exists(os.path.join(scratch_dir, 'QMfiles/Reactions', family, 'false_positives', my_folder)):
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
                print "Not sure about reaction {0}, has contents {1}".format(my_folder, 
                            ' , '.join(os.listdir(os.path.join(scratch_dir, 'QMfiles/Reactions', 
                                                               family, my_folder, 'm062x'))))

    return success


# In[7]:

print get_success(rxn)


# In[8]:

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
        if rxn_family.strip() == family:
            reaction = get_rmg_reaction(rmg_database, species_dict, family, reaction_line)

            if reaction is not None:
                features_for_reaction = get_features(reaction)
                features_for_reaction['Success Code'] = get_success(reaction)
                data.append(features_for_reaction)
            else:
                print "No RMG reactions found for reaction {0}".format(reaction_line.split()[0])

    return pd.DataFrame().from_dict(data)


# In[9]:

family_data = get_family_data('Silylene_Insertion')


# In[10]:

print family_data.head()


# In[11]:

family_data = family_data[family_data['Success Code'].notnull()]
print family_data.shape


# In[118]:

import numpy as np
import matplotlib.pyplot as plt
import seaborn
get_ipython().magic(u'matplotlib inline')


# In[13]:

family_data.groupby('Success Code').size().plot(kind='bar')
plt.xlabel('Success Code')
plt.ylabel('Frequency')
plt.show()


# In[14]:

average_H2_value = family_data.groupby('Success Code')['H2', 'Sirad_H', 'sil_H2'].mean()
average_H2_value.plot.bar()


# In[29]:

from sklearn.cross_validation import ShuffleSplit, cross_val_score, train_test_split, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix


# In[58]:

y = family_data['Success Code']
X = family_data[[column for column in family_data.columns if column not in ['Success Code', 'reaction']]]


# In[59]:

log_reg = LogisticRegression(class_weight='balanced')#)#LogisticRegression(multi_class='multinomial', solver='lbfgs')
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=9)
log_reg.fit(X_train, y_train)


# In[48]:

log_reg.score(X_test, y_test)


# In[49]:

print log_reg.coef_
print log_reg.classes_
print X.columns


# In[50]:

zip(y_test, log_reg.predict(X_test))


# Try doing StratifiedKFold for cross validation- keep training and test set with proportional number of class labels 

# In[63]:

cv_skf = StratifiedKFold(y, random_state=9)
sum_scores = 0.0
for train_index, test_index in cv_skf:
    X_train, X_test = X.iloc[train_index,:], X.iloc[test_index,:]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    log_reg.fit(X_train, y_train)
    sum_scores += log_reg.score(X_test, y_test)
    print log_reg.score(X_test, y_test)
    print zip(y_test, log_reg.predict(X_test))
    print '\n'
print "Avg score: {0}".format(round(sum_scores/len(cv_skf), 3))


# In[64]:

predictions = log_reg.predict(X)
performance = pd.DataFrame(data=zip(y, predictions), columns = ['Actual', 'Prediction'])


# In[108]:

cm = confusion_matrix(y, predictions, labels=['S', 'FP', 'F1', 'F2'])
cm_df = pd.DataFrame(data=cm, columns=['S', 'FP', 'F1', 'F2'])
print cm_df


# In[109]:

percents = cm_df.iloc[:].transpose().apply(lambda x: x/x.sum()).transpose()


# In[120]:

print percents
percents.plot.bar()
plt.xticks(np.arange(4), ('S', 'FP', 'F1', 'F2') )


# #### Do a grid search to figure out best parameters for model (ie size of test set etc.)
# #### Include more features, remove unnecessary (Sid_Si_H is always False for this set)

# Remove features which may not be independent
# 
# For example, if H2 is true, then all the Si_H ones are false

# In[60]:

X = family_data[[column for column in family_data.columns if column not in ['Success Code', 'reaction', 'H2', 'sil_H2']]]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state = 9)
log_reg.fit(X_train, y_train)
log_reg.score(X_test, y_test)

