import pandas as pd
import numpy as np
import os
from damagerepair_counting import counting, prune

if not os.path.exists('../datasets/'):
    os.makedirs('../datasets/')

##### read data
data = pd.read_csv('../../step_data/Exercise.csv')

# fix age units
data['Age (months)'] = data['Age (weeks)']/4.34524


data['mouse numeric'] = data['Mouse Code'].astype(int)
data = data.sort_values(by=['mouse numeric','Age (months)'])
data = data.drop('mouse numeric', axis=1)

# fix death ages
data['status'] = np.nan
data['death age'] = np.nan
for label, group in data.groupby('Mouse Code'):
    if group.shape == group.dropna(how='all',subset=group.columns.values[1:]).shape:
        data.loc[data['Mouse Code'] == label, 'death.age'] = group['Age (months)'].max() - group['Age (months)'].min()
        data.loc[data['Mouse Code'] == label, 'status'] = 0
    else:
        data.loc[data['Mouse Code'] == label, 'death.age'] = group['Age (months)'].max() - group['Age (months)'].min()
        data.loc[data['Mouse Code'] == label, 'status'] = 1
        
data = data.dropna(how = 'all', subset = ['Age (weeks)'])

#remove missing rows
data = data.dropna(subset = data.columns.values[1:], how='all')#, 
data = data.dropna(subset = data.columns.values[11:], how='all')#, 

individuals = []

# remove short measurement times
for label, group in data.groupby(['Mouse Code']):
    if np.diff(group['Age (months)']).min() < 0.1:
        diff = np.diff(group['Age (months)'])
        ages = group['Age (months)'].values
        ids = [i  for i in range(group.shape[0]) if i != group.shape[0]-2]
        individuals.append(group.iloc[ids])
    else:
        individuals.append(group)
data = pd.concat(individuals, axis=0)
data = data.loc[data['Age (months)'] > 0]


data['exer'] = data['Mouse Group (1=sedentary, 2= exercise)']
deficits = list(data.columns.values[12:43])

# impute
for label, group in data.groupby(['Mouse Code']):
    data.loc[data['Mouse Code'] == label, deficits] = data.loc[data['Mouse Code'] == label, deficits].ffill()
    data.loc[data['Mouse Code'] == label, deficits] = data.loc[data['Mouse Code'] == label, deficits].bfill()

data, pruned_repair, pruned_damage = prune(data, deficits, 'Mouse Code', 'Age (months)')
#data_prune.to_csv('mouse_2_pruned.csv')


import matplotlib.pyplot as plt


fig, ax = plt.subplots(2,1,figsize=(5, 6))
ax=ax.flatten()

ax[0].title.set_text('i) Mouse dataset 2')

pruned_repair, counts = np.unique(pruned_repair, return_counts=True)
index = np.argsort(counts)[::-1]
ax[0].bar(np.arange(pruned_repair.shape[0]), counts[index], tick_label = pruned_repair[index])
ax[0].set_xticklabels(pruned_repair[index], rotation='vertical', fontsize=7)

output_repair = pd.DataFrame({'deficit': pruned_repair[index], 'pruned count': counts[index]})

pruned_damage, counts = np.unique(pruned_damage, return_counts=True)
index = np.argsort(counts)[::-1]
ax[1].bar(np.arange(pruned_damage.shape[0]), counts[index], tick_label = pruned_damage[index])
ax[1].set_xticklabels(pruned_damage[index], rotation='vertical', fontsize=7)

output_damage = pd.DataFrame({'deficit': pruned_damage[index], 'pruned count': counts[index]})

ax[0].set_xlabel('')
ax[0].set_ylabel('Isolated repair pruned')

ax[1].set_xlabel('')
ax[1].set_ylabel('Isolated damage pruned')

plt.tight_layout()
plt.savefig('Mouse2_pruned.pdf')

output_repair.to_csv('../figure_data/figure5_supplement4/mouse2_pruned_repair.csv')
output_damage.to_csv('../figure_data/figure5_supplement4/mouse2_pruned_damage.csv')


# create binary deficits from fractional
new_deficits = []
for d in deficits:
    
    for i in range(4):
        new_deficits.append(d + str(i))
        data[d+str(i)] = 0
    
    for j in range(data[d].shape[0]):
        i=0
        while int(data.iloc[j].loc[d] * 4) > i:
            data[d+str(i)].iloc[j] = 1
            i += 1
deficits = new_deficits

data = counting(data, new_deficits, 'Mouse Code', 'Age (months)')


data['mouse'] = data['Mouse Code']
data['n'] = data['deficit.count']
data['repair'] = data['repair.count']
data['damage'] = data['damage.count']
data['N'] = data['total.deficits']
data['f'] = data['FI']
data['sex'] = data['Sex']
data['exercise'] = data['exer'].apply(lambda x: 'yes' if x == 2 else 'no')

data = data[['mouse', 'sex', 'exercise', 'baseline.age', 'time', 
               'f','repair', 'damage', 'n','N','delta.t', 'status', 'death.age'] + \
              ['d%d'%i for i in range(len(deficits))] +\
              ['d%d.r'%i for i in range(len(deficits))] +\
              ['d%d.d'%i for i in range(len(deficits))]]


data['mouse numeric'] = data['mouse'].astype(int)
data['real age'] = data['baseline.age'] + data['time']
data = data.sort_values(by=['mouse numeric','real age'])
data = data.drop('mouse numeric', axis=1)
data = data.drop('real age', axis=1)

data.to_csv('../datasets/exercise_data_pruned.csv', index=False)


surv_data = []
for label, group in data.groupby('mouse'):
     surv_data.append(group[['mouse','sex','exercise','baseline.age', 'death.age', 'status']].iloc[0])
surv_data = pd.DataFrame(surv_data)

surv_data.to_csv('../datasets/exercise_surv_data_pruned.csv', index=False)
