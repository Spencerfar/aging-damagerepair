import pandas as pd
import numpy as np
import os

if not os.path.exists('../datasets/'):
    os.makedirs('../datasets/')

#####data
males = pd.read_csv('../../step_data/Enalapril_male.csv')
females = pd.read_csv('../../step_data/Enalapril_female.csv')

male_deaths = pd.read_csv('../../Male mice deaths.csv')
female_deaths = pd.read_csv('../../Female mice deaths.csv')

##### remove some mice
# unknown death age
females = females[females['Mouse Code'] != 'F19.1']
female_deaths = female_deaths[female_deaths['Mouse Code'] != 'F19.1']

# only 1 time point
females = females[females['Mouse Code'] != 'F5.2']
female_deaths = female_deaths[female_deaths['Mouse Code'] != 'F5.2']

# only 1 time point
females = females[females['Mouse Code'] != 'F8.2']
female_deaths = female_deaths[female_deaths['Mouse Code'] != 'F8.2']


male_deaths.loc['Age at death (days)']= pd.to_numeric(male_deaths['Age at death (days)'], downcast="float")
female_deaths.loc['Age at death (days)']= pd.to_numeric(female_deaths['Age at death (days)'], downcast="float")

# fix months units
females['Age (months) rounded'] = round(females['Age (weeks)']/4.34524).astype(int)
females['Age (months)'] = females['Age (weeks)']/4.34524
males['Age (months) rounded']=round(males['Age (weeks)']/4.34524).astype(int)
males['Age (months)'] = males['Age (weeks)']/4.34524


females = females.sort_values(by=['Mouse Code', 'Age (months)'])
female_deaths = female_deaths.sort_values(by=['Mouse Code'])
male_deaths['status'] = male_deaths['Died before 25 months FI (1=yes,0=no)']
female_deaths['status'] = female_deaths['Died before 23 months FI (1=yes,0=no)']

# fix death ages
for label, group in males.groupby('Mouse Code'):
    
    if male_deaths.loc[male_deaths['Mouse Code'] == label, 'status'].values == 0:
        male_deaths.loc[male_deaths['Mouse Code'] == label, 'Age at death (months)'] = \
            group['Age (months)'].max()

for label, group in males.groupby('Mouse Code'):
    
    if male_deaths.loc[male_deaths['Mouse Code'] == label, 'status'].values == 0:
        male_deaths.loc[male_deaths['Mouse Code'] == label, 'Age at death (months)'] = \
            group['Age (months)'].max()


for label, group in females.groupby('Mouse Code'):
    if female_deaths.loc[female_deaths['Mouse Code'] == label, 'status'].values == 0:
        female_deaths.loc[female_deaths['Mouse Code'] == label, 'Age at death (months)'] = \
            group['Age (months)'].max()

# fix ages for drugs
males = males.loc[males['Age (months)'] <= 25]
females = females.loc[females['Age (months)'] <= 23]
males['treat'] = males['Mouse Group (1=drug, 3=control, 2=off drug)']
females['treat'] = females['Mouse Group (1=drug, 3=control)']

for data in [males, females]:
    keys = data['Mouse Code'].unique()
    deficits = list(data.columns.values[12:43])

# impute
for data in (males, females):
    for label, group in data.groupby(['Mouse Code']):
        
        data.loc[data['Mouse Code'] == label, deficits] = data.loc[data['Mouse Code'] == label, deficits].ffill()
        data.loc[data['Mouse Code'] == label, deficits] = data.loc[data['Mouse Code'] == label, deficits].bfill()


# create binary deficits from fractional
new_male_deficits = []
for d in deficits:
    
    for i in range(4):
        new_male_deficits.append(d + str(i))
        males[d+str(i)] = 0
    
    for j in range(males[d].shape[0]):
        i=0
        while int(males.iloc[j].loc[d] * 4) > i:
            males[d+str(i)].iloc[j] = 1
            i += 1

new_female_deficits = []
for d in deficits:
    
    for i in range(4):
        new_female_deficits.append(d + str(i))
        females[d+str(i)] = 0
    
    for j in range(females[d].shape[0]):
        i=0
        while int(females.iloc[j].loc[d] * 4) > i:
            females[d+str(i)].iloc[j] = 1
            i += 1

deficits = new_female_deficits

##### repair
results_status = {}
results_damage = {}
results_repair = {}
results_sex = {}
results_treatment = {}
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
    baseline = []
    treatment = []
    index = []
    diff_left = []
    diff_right = []
    
    for (data, deaths) in [(males, male_deaths), (females, female_deaths)]:
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
                    
                    dam_age_prev = group['Age (months)'].values[s]
                    rep_age_prev = group['Age (months)'].values[s+1]
                    
                    rep = 1 - group[d].values[s+2]
                    
                    
                elif len(end) > 0 and s+1<end[0]:
                    dam_age = group['Age (months)'].values[s+1]
                    rep_age = group['Age (months)'].values[end[0]+1]
                    
                    dam_age_prev = group['Age (months)'].values[s]
                    rep_age_prev = group['Age (months)'].values[end[0]]
                    
                    rep = 1 - group[d].values[end[0]+1]
                    
                    
                elif s + 1 < len(group[d])-1:
                    dam_age = group['Age (months)'].values[s+1]
                    rep_age = group['Age (months)'].values[-1]
                    
                    dam_age_prev = group['Age (months)'].values[s]
                    rep_age_prev = group['Age (months)'].values[-2]
                    
                    rep = 1 - group[d].values[-1]
                    
                else:
                    continue
                
                if not np.isnan(rep):
                    status.append(rep)
                    damage_age.append(dam_age)
                    repair_age.append(rep_age)
                    damage_age_prev.append(dam_age_prev)
                    repair_age_prev.append(rep_age_prev)
                    baseline.append(group['Age (months)'].iloc[0])
                    sex.append(group['Sex'].values[0])
                    treatment.append('drug' if group['treat'].iloc[0] < 3 else 'control')
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
    results_baseline[d] = np.array(baseline)
    
    results_sex[d] = np.array(sex)
    results_treatment[d] = np.array(treatment)
    results_id[d] = np.array(index)

all_df = []

for i,d in enumerate(deficits):
    
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], results_treatment[d],results_baseline[d],
                       results_damage[d], results_repair[d], results_status[d],
                                results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['mouse', 'sex', 'treatment', 'baseline.age', 'damage.age', 'repair.age',
                                 'status', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/enalapril_repair_deficits.csv')





##### damage
results_status = {}
results_damage = {}
results_repair = {}
results_sex = {}
results_treatment = {}
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
    baseline = []
    treatment = []
    index = []
    diff_left = []
    diff_right = []
    
    for (data, deaths) in [(males, male_deaths), (females, female_deaths)]:
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
                rep = np.nan
                
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
                    baseline.append(group['Age (months)'].values[0])
                    sex.append(group['Sex'].values[0])
                    treatment.append('drug' if group['treat'].iloc[0] < 3 else 'control')
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
    results_baseline[d] = np.array(baseline)
    results_diff_left[d] = np.array(diff_left)
    results_diff_right[d] = np.array(diff_right)
    results_id[d] = np.array(index)

all_df = []

for i,d in enumerate(deficits):
    df = pd.DataFrame(np.array([results_id[d], results_sex[d], results_treatment[d],results_baseline[d],
                       results_damage[d], results_repair[d], results_status[d],
                               results_diff_left[d], results_diff_right[d]]).T,
                      columns = ['mouse', 'sex', 'treatment', 'baseline.age','damage.age', 'repair.age',
                                 'status', 'diff.left', 'diff.right'])
    df['deficit'] = d
    all_df.append(df)

all_df = pd.concat(all_df)
all_df.to_csv('../datasets/enalapril_damage_deficits.csv')
