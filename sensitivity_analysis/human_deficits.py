import pandas as pd
import numpy as np
import os
from pyreadstat import read_sav

from damagerepair_counting import counting, prune

if not os.path.exists('../datasets/'):
    os.makedirs('../datasets/')

# read data
#data = pd.read_csv('../datasets/human_data.csv')
data = pd.read_csv('../clean_data/extracted_elsa.csv')
folder = "~/Downloads/UKDA-5050-spss/spss/spss24/"
fin = read_sav(folder + "wave_5_financial_derived_variables.sav", 
               usecols = ['idauniq', 'totinc_bu_s', 'totinc_bu_t', 
                          'nettotw_bu_s', 'nettotw_bu_t', 'nettotw_bu_ni2'])[0]

# remove imputed data
fin = fin[fin['nettotw_bu_t'] <= 2.0]
fin = fin[fin['nettotw_bu_ni2'] <= 0.0]
fin = fin[fin['nettotw_bu_t'] != -995.0]
fin = fin[fin['nettotw_bu_t'] != -1.0]
fin = fin[fin['nettotw_bu_t'] != -999.0]
fin = fin[fin['nettotw_bu_t'] != -998.0]

data['wealth'] = np.nan
for label, group in data.groupby('id'):
    if label in fin['idauniq'].values:
        
        data.loc[data['id']==label, 'wealth'] = \
            fin.loc[fin['idauniq']==label,'nettotw_bu_s'].item()

N_adl = 10
N_iadl = 13
N = N_adl + N_iadl
deficits_n = ['d'+str(i) for i in range(N_adl + N_iadl)]
deficits_r = ['d'+str(i)+'.r' for i in range(N_adl + N_iadl)]
deficits_d = ['d'+str(i)+'.d' for i in range(N_adl + N_iadl)]


data[deficits_n] = data[deficits_n].replace(-1, np.nan)
data = data.dropna(subset=deficits_n)

data['n'] = data[deficits_n].sum(1)
data['f'] = data['n']/N

"""
data = data.groupby('id').filter(lambda x: x.shape[0] > 6)
#data = counting(data, deficits_n, 'id', 'age')

for label, group in data.groupby(['id']):

    #time = group[time_column] - group[time_column].iloc[0]

    data.loc[data['id'] == label, 'delta.t'] = np.append(np.diff(group['age']), np.nan)
    
data = data.groupby('id').filter(lambda x: np.all((x['delta.t'] <= 4) | np.isnan(x['delta.t'])) )

data = data.dropna(how='any',subset=['wealth'])

#data = data.replace(-1000.0, np.nan)
"""

# prune
data, pruned_repair, pruned_damage = prune(data, deficits_n, 'id', 'age')

cols = ['id', 'wave', 'age', 'sex', 'n', 'wealth'] + deficits_n
data = data[cols]
data = data.groupby('id').filter(lambda x: x.shape[0] > 6)
deficits = deficits_n


data = counting(data, deficits_n, 'id', 'age')

data['n'] = data['deficit.count']
data['repair'] = data['repair.count']
data['damage'] = data['damage.count']
data['N'] = 23.0
data['f'] = data['FI']
data['age'] = data['time'] + data['baseline.age']

data = data[['id', 'wave', 'sex', 'age', 'baseline.age', 'time', 
             'repair', 'damage', 'n', 'N', 'delta.t', 'wealth'] + \
              ['d%d'%i for i in range(len(deficits))]]


data = data.groupby('id').filter(lambda x: np.all((x['delta.t'] <= 4) | np.isnan(x['delta.t'])) )

data = data.dropna(how='any',subset=['wealth'])

data = data.replace(-1000.0, np.nan)

#deficits = ['d'+str(i)  for i in range(23)]


##### repair
results_status = {}
results_damage = {}
results_repair = {}
results_sex = {}
results_wealth = {}
results_baseline = {}
results_id = {}
results_diff_left = {}
results_diff_right = {}

for d in deficits:
    
    status = []
    damage_age = []
    damage_age_prev = []
    repair_age = []
    repair_age_prev = []
    sex = []
    wealth = []
    baseline = []
    index = []
    diff_left = []
    diff_right = []
    
    if True:
        for label, group in data.groupby(['id']):
            start = np.where(((group[d].iloc[:-1] == 0).values & \
                   (group[d].iloc[1:] == 1).values).astype(int) > 0)[0]
            
            for s in start:
                end = np.where(((group[d].iloc[:-1] == 1).values & \
                     (group[d].iloc[1:] == 0).values).astype(int) > 0)[0]
                end = end[end>s]
                
                dam_age = np.nan
                dam_age_prev = np.nan
                rep_age = np.nan
                rep_age_prev = np.nan
                rep = np.nan
                
                if len(end) > 0 and s+1==end[0]:
                    dam_age = group['age'].values[s+1]
                    rep_age = group['age'].values[s+2]
                    
                    dam_age_prev = group['age'].values[s]
                    rep_age_prev = group['age'].values[s+1]
                    rep = 1 - group[d].values[s+2]
                    
                elif len(end) > 0 and s+1<end[0]:
                    dam_age = group['age'].values[s+1]
                    rep_age = group['age'].values[end[0]+1]
                    
                    dam_age_prev = group['age'].values[s]
                    rep_age_prev = group['age'].values[end[0]]
                    
                    rep = 1 - group[d].values[end[0]+1]
                    
                    
                elif s + 1 < len(group[d])-1:
                    dam_age = group['age'].values[s+1]
                    rep_age = group['age'].values[-1]
                    
                    dam_age_prev = group['age'].values[s]
                    rep_age_prev = group['age'].values[-2]
                    rep = 1 - group[d].values[-1]
                    
                else:
                    continue
                
                if not np.isnan(rep):
                    status.append(rep)
                    damage_age.append(dam_age)
                    repair_age.append(rep_age)
                    damage_age_prev.append(dam_age_prev)
                    repair_age_prev.append(rep_age_prev)
                    sex.append(group['sex'].values[0])
                    wealth.append(group['wealth'].values[0])
                    baseline.append(group['baseline.age'].values[0])
                    index.append(label)
                    
                    if rep > 0:
                        diff_left.append(rep_age_prev - dam_age)
                        diff_right.append(rep_age - dam_age_prev)
                    else:
                        diff_left.append(rep_age_prev - dam_age)
                        diff_right.append(np.inf)
                    
                else:
                    print('nan repair')
                    
    results_status[d] = np.array(status)
    results_damage[d] = np.array(damage_age)
    results_repair[d] = np.array(repair_age)
    results_diff_left[d] = np.array(diff_left)
    results_diff_right[d] = np.array(diff_right)
    
    results_sex[d] = np.array(sex)
    results_wealth[d] = np.array(wealth)
    results_baseline[d] = np.array(baseline)
    results_id[d] = np.array(index)

all_df = []

for i,d in enumerate(deficits):
    
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], results_wealth[d],
                       results_damage[d], results_repair[d], 
                       results_status[d], results_baseline[d],
                                results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['id', 'sex', 'wealth', 'damage.age', 'repair.age',
                                 'status', 'baseline.age','diff.left','diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/human_repair_deficits_pruned.csv')



##### damage
results_status = {}
results_damage = {}
results_repair = {}
results_sex = {}
results_treatment = {}
results_id = {}
results_diff_left = {}
results_diff_right = {}

for d in deficits:
    status = []
    damage_age = []
    damage_age_prev = []
    repair_age = []
    repair_age_prev = []
    sex = []
    baseline = []
    wealth = []
    index = []
    diff_left = []
    diff_right = []
    
    if True:
        for label, group in data.groupby(['id']):
            start = np.where(((group[d].iloc[:-1] == 1).values & \
                   (group[d].iloc[1:] == 0).values).astype(int) > 0)[0]
            
            for s in start:
                end = np.where(((group[d].iloc[:-1] == 0).values & \
                     (group[d].iloc[1:] == 1).values).astype(int) > 0)[0]
                end = end[end>s]
                
                dam_age = np.nan
                dam_age_prev = np.nan
                rep_age = np.nan
                rep_age_prev = np.nan
                rep = np.nan
                
                if len(end) > 0 and s+1==end[0]:
                    rep_age = group['age'].values[s+1]
                    dam_age = group['age'].values[s+2]
                    
                    rep_age_prev = group['age'].values[s]
                    dam_age_prev = group['age'].values[s+1]
                    dam = group[d].values[s+2]
                    
                elif len(end) > 0 and s+1<end[0]:
                    rep_age = group['age'].values[s+1]
                    dam_age = group['age'].values[end[0]+1]
                    rep_age_prev = group['age'].values[s]
                    dam_age_prev = group['age'].values[end[0]]
                    dam = group[d].values[end[0]+1]
                    
                elif s + 1 < len(group[d])-1:
                    rep_age = group['age'].values[s+1]
                    dam_age = group['age'].values[-1]
                    
                    rep_age_prev = group['age'].values[s]
                    dam_age_prev = group['age'].values[-2]
                    dam = group[d].values[-1]
                    
                else:
                    continue
                
                if not np.isnan(dam):
                    status.append(dam)
                    damage_age.append(dam_age)
                    repair_age.append(rep_age)
                    damage_age_prev.append(dam_age_prev)
                    repair_age_prev.append(rep_age_prev)
                    sex.append(group['sex'].values[0])
                    wealth.append(group['wealth'].values[0])
                    baseline.append(group['baseline.age'].values[0])
                    index.append(label)
                    
                    if dam > 0:
                        diff_left.append(dam_age_prev - rep_age)
                        diff_right.append(dam_age - rep_age_prev)
                    else:
                        diff_left.append(dam_age_prev - rep_age)
                        diff_right.append(np.inf)
                    
    results_status[d] = np.array(status)
    results_damage[d] = np.array(damage_age)
    results_repair[d] = np.array(repair_age)
    results_diff_left[d] = np.array(diff_left)
    results_diff_right[d] = np.array(diff_right)
    
    results_sex[d] = np.array(sex)
    results_wealth[d] = np.array(wealth)
    results_baseline[d] = np.array(baseline)
    results_id[d] = np.array(index)

all_df = []

for i,d in enumerate(deficits):
    print(results_id[d].shape, results_sex[d].shape, results_wealth[d].shape, results_baseline[d].shape)
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], results_wealth[d],
                       results_damage[d], results_repair[d], 
                       results_status[d], results_baseline[d],
                               results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['id', 'sex', 'wealth', 'damage.age', 'repair.age',
                                 'status', 'baseline.age', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/human_damage_deficits_pruned.csv')
