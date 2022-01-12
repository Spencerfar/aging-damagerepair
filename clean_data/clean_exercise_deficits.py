import pandas as pd
import numpy as np

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

#remove missing rows
data = data.dropna(how = 'all', subset = ['Age (weeks)'])
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

##### repair
results_status = {}
results_damage = {}
results_repair = {}
results_sex = {}
results_exercise = {}
results_id = {}
results_baseline = {}
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
    exercise = []
    index = []
    diff_left = []
    diff_right = []
    
    for label, group in data.groupby(['Mouse Code']):
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
                dam_age = group['Age (months)'].values[s+1]
                rep_age = group['Age (months)'].values[s+2]
                rep = 1 - group[d].values[s+2]
                
                dam_age_prev = group['Age (months)'].values[s]
                rep_age_prev = group['Age (months)'].values[s+1]
                    
            elif len(end) > 0 and s+1<end[0]:
                dam_age = group['Age (months)'].values[s+1]
                rep_age = group['Age (months)'].values[end[0]+1]
                rep = 1 - group[d].values[end[0]+1]
                
                dam_age_prev = group['Age (months)'].values[s]
                rep_age_prev = group['Age (months)'].values[end[0]]
                    
            elif s + 1 < len(group[d])-1:
                dam_age = group['Age (months)'].values[s+1]
                rep_age = group['Age (months)'].values[-1]
                rep = 1 - group[d].values[-1]
                
                dam_age_prev = group['Age (months)'].values[s]
                rep_age_prev = group['Age (months)'].values[-2]
                
            else:
                continue
                
            if not np.isnan(rep):
                status.append(rep)
                damage_age.append(dam_age)
                repair_age.append(rep_age)
                damage_age_prev.append(dam_age_prev)
                repair_age_prev.append(rep_age_prev)
                sex.append(group['Sex'].values[0])
                baseline.append(group['Age (months)'].iloc[0])
                exercise.append('yes' if group['exer'].iloc[0] > 1 else 'no')
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
    
    results_baseline[d] = np.array(baseline)
    results_sex[d] = np.array(sex)
    results_exercise[d] = np.array(exercise)
    results_id[d] = np.array(index)

all_df = []
for i,d in enumerate(deficits):
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], results_exercise[d],results_baseline[d],
                       results_damage[d], results_repair[d], results_status[d],
                               results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['mouse', 'sex', 'exercise', 'baseline.age', 'damage.age', 'repair.age',
                                 'status', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/exercise_repair_deficits.csv')


##### damage
results_status = {}
results_damage = {}
results_repair = {}
results_sex = {}
results_exercise = {}
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
    exercise = []
    index = []
    diff_left = []
    diff_right = []
    
    for label, group in data.groupby(['Mouse Code']):
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
                rep_age = group['Age (months)'].values[s+1]
                dam_age = group['Age (months)'].values[s+2]
                dam = group[d].values[s+2]
                rep_age_prev = group['Age (months)'].values[s]
                dam_age_prev = group['Age (months)'].values[s+1]
                    
            elif len(end) > 0 and s+1<end[0]:
                rep_age = group['Age (months)'].values[s+1]
                dam_age = group['Age (months)'].values[end[0]+1]
                dam = group[d].values[end[0]+1]
                
                rep_age_prev = group['Age (months)'].values[s]
                dam_age_prev = group['Age (months)'].values[end[0]]
                    
            elif s + 1 < len(group[d])-1:
                rep_age = group['Age (months)'].values[s+1]
                dam_age = group['Age (months)'].values[-1]
                dam = group[d].values[-1]
                
                rep_age_prev = group['Age (months)'].values[s]
                dam_age_prev = group['Age (months)'].values[-2]
                
            else:
                continue
                
            if not np.isnan(dam):
                status.append(dam)
                damage_age.append(dam_age)
                repair_age.append(rep_age)
                damage_age_prev.append(dam_age_prev)
                repair_age_prev.append(rep_age_prev)
                sex.append(group['Sex'].values[0])
                exercise.append('yes' if group['exer'].iloc[0] > 1 else 'no')
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
    results_exercise[d] = np.array(exercise)
    results_id[d] = np.array(index)
    results_diff_left[d] = np.array(diff_left)
    results_diff_right[d] = np.array(diff_right)

all_df = []
for i,d in enumerate(deficits):
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], results_exercise[d],
                       results_damage[d], results_repair[d], results_status[d],
                               results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['mouse', 'sex', 'exercise', 'damage.age', 'repair.age',
                                 'status', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/exercise_damage_deficits.csv')
