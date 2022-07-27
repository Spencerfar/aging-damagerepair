import pandas as pd
import numpy as np
import argparse
import os
from damagerepair_counting import counting, prune

if not os.path.exists('../datasets/'):
    os.makedirs('../datasets/')

parser = argparse.ArgumentParser('Clean Schultz')
parser.add_argument('--folder', type=str, default = '../../step_data/')
args = parser.parse_args()

##### read data
data = pd.read_csv(args.folder + 'Schultz.csv')

# remove single times
single_times = []
for label, group in data.groupby('Mouse ID'):
    if len(group) <= 2:
        single_times.append(label)
data = data.loc[~data['Mouse ID'].isin(single_times)]

deficits = data.columns.values[1:30]

# fixing age units
data['age_months'] = data['age_days']/30.4167


data['Age at death (months)'] = np.nan
for label, group in data.groupby('Mouse ID'):
    
    death_age = group['age_months'].max() + group['Lifespan from point of FI test (days)'].min()/30.4167
    data.loc[data['Mouse ID'] == label, 'Age at death (months)'] = death_age


data, pruned_repair, pruned_damage = prune(data, list(deficits), 'Mouse ID', 'age_months')
#data_prune.to_csv('mouse_2_pruned.csv')


import matplotlib.pyplot as plt


fig, ax = plt.subplots(2,1,figsize=(5, 6))
ax=ax.flatten()

ax[0].title.set_text('j) Mouse dataset 3')

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
plt.savefig('Mouse3_pruned.pdf')

output_repair.to_csv('../figure_data/figure5_supplement4/mouse3_pruned_repair.csv')
output_damage.to_csv('../figure_data/figure5_supplement4/mouse3_pruned_damage.csv')


    
# create binary deficits from fractional
new_deficits = []
for d in deficits:

    new_deficit_columns = [d+str(i) for i in range(4)]
    new_cols = pd.DataFrame(np.zeros((data.shape[0], 4)), columns=new_deficit_columns, index=data.index)
    
    data = pd.concat((data, new_cols), axis=1)
    new_deficits = new_deficits + new_deficit_columns
    
    for j in range(data[d].shape[0]):
        i=0
        while int(data.iloc[j].loc[d] * 4) > i:
            data.loc[data.index[j],d+str(i)] = 1
            i += 1
deficits=new_deficits


data = counting(data, new_deficits, 'Mouse ID', 'age_months')

data['mouse'] = data['Mouse ID']
data['n'] = data['deficit.count']
data['repair'] = data['repair.count']
data['damage'] = data['damage.count']
data['N'] = data['total.deficits']
data['f'] = data['FI']
data['sex'] = 'Male'
data['age'] = data['time'] + data['baseline.age']
data['status'] = 1
data['death.age'] = data['Age at death (months)']


# remove short dt
data = data.loc[(data['delta.t']<3) | np.isnan(data['delta.t'])]

# recompute times
for label, group in data.groupby(['mouse']):
    time = group['age'] - group['age'].iloc[0]
    data.loc[data['mouse'] == label, 'time'] = time
    data.loc[data['mouse'] == label, 'baseline.age'] = group['age'].iloc[0]


data = data[['mouse', 'sex', 'baseline.age', 'time', 
               'f', 'repair', 'damage', 'n','N','delta.t', 'death.age',
               'status'] + \
              ['d%d'%i for i in range(len(deficits))] +\
              ['d%d.r'%i for i in range(len(deficits))] +\
              ['d%d.d'%i for i in range(len(deficits))]]


data.to_csv('../datasets/schultz_data_pruned.csv', index=False)

surv_data = []
for label, group in data.groupby('mouse'):
     surv_data.append(group[['mouse','sex','baseline.age', 'death.age', 'status']].iloc[0])
surv_data = pd.DataFrame(surv_data)

surv_data.to_csv('../datasets/schultz_surv_data_pruned.csv', index=False)
