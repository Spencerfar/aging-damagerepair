import pandas as pd
import numpy as np
import os

if not os.path.exists('../datasets/'):
    os.makedirs('../datasets/')

##### read data
data = pd.read_csv('../../step_data/Schultz.csv')

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

deficits=new_deficits

##### repair
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
    treatment = []
    index = []
    diff_left = []
    diff_right = []
    
    if True:
        for label, group in data.groupby(['Mouse ID']):
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
                    dam_age = group['age_months'].values[s+1]
                    rep_age = group['age_months'].values[s+2]
                    rep = 1 - group[d].values[s+2]
                    
                    dam_age_prev = group['age_months'].values[s]
                    rep_age_prev = group['age_months'].values[s+1]
                    
                    
                elif len(end) > 0 and s+1<end[0]:
                    dam_age = group['age_months'].values[s+1]
                    rep_age = group['age_months'].values[end[0]+1]
                    rep = 1 - group[d].values[end[0]+1]
                    
                    dam_age_prev = group['age_months'].values[s]
                    rep_age_prev = group['age_months'].values[end[0]]
                    
                elif s + 1 < len(group[d])-1:
                    dam_age = group['age_months'].values[s+1]
                    rep_age = group['age_months'].values[-1]
                    rep = 1 - group[d].values[-1]
                    
                    dam_age_prev = group['age_months'].values[s]
                    rep_age_prev = group['age_months'].values[-2]
                    
                else:
                    continue
                
                if not np.isnan(rep):
                    status.append(rep)
                    damage_age.append(dam_age)
                    repair_age.append(rep_age)
                    sex.append('Male')
                    index.append(label)
                    
                    if rep > 0:
                        diff_left.append(rep_age_prev - dam_age)
                        diff_right.append(rep_age - dam_age_prev)
                    else:
                        diff_left.append(rep_age_prev - dam_age)
                        diff_right.append(np.inf)
                    
                    
    results_status[d] = np.array(status)
    results_damage[d] = np.array(damage_age)
    results_repair[d] = np.array(repair_age)
    
    results_diff_left[d] = np.array(diff_left)
    results_diff_right[d] = np.array(diff_right)
    
    results_sex[d] = np.array(sex)
    results_treatment[d] = np.array(treatment)
    results_id[d] = np.array(index)

all_df = []
for i,d in enumerate(deficits):
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], 
                       results_damage[d], results_repair[d], results_status[d],
                               results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['mouse', 'sex', 'damage.age', 'repair.age',
                                 'status', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/schultz_repair_deficits.csv')



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
    treatment = []
    index = []
    diff_left = []
    diff_right = []
    
    if True:
        for label, group in data.groupby(['Mouse ID']):
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
                dam = np.nan
                
                if len(end) > 0 and s+1==end[0]:
                    rep_age = group['age_months'].values[s+1]
                    dam_age = group['age_months'].values[s+2]
                    dam = group[d].values[s+2]
                    
                    rep_age_prev = group['age_months'].values[s]
                    dam_age_prev = group['age_months'].values[s+1]
                    
                elif len(end) > 0 and s+1<end[0]:
                    rep_age = group['age_months'].values[s+1]
                    dam_age = group['age_months'].values[end[0]+1]
                    dam = group[d].values[end[0]+1]
                    
                    rep_age_prev = group['age_months'].values[s]
                    dam_age_prev = group['age_months'].values[end[0]]
                    
                elif s + 1 < len(group[d])-1:
                    rep_age = group['age_months'].values[s+1]
                    dam_age = group['age_months'].values[-1]
                    dam = group[d].values[-1]
                    
                    rep_age_prev = group['age_months'].values[s]
                    dam_age_prev = group['age_months'].values[-2]
                    
                else:
                    continue
                
                if not np.isnan(dam):
                    status.append(dam)
                    damage_age.append(dam_age)
                    repair_age.append(rep_age)
                    damage_age_prev.append(dam_age_prev)
                    repair_age_prev.append(rep_age_prev)
                    sex.append('Male')
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
    
    results_sex[d] = np.array(sex)
    results_treatment[d] = np.array(treatment)
    results_id[d] = np.array(index)
    results_diff_left[d] = np.array(diff_left)
    results_diff_right[d] = np.array(diff_right)

all_df = []

for i,d in enumerate(deficits):
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], 
                       results_damage[d], results_repair[d], results_status[d],
                               results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['mouse', 'sex', 'damage.age', 'repair.age',
                                 'status', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/schultz_damage_deficits.csv')
