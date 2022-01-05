import pandas as pd
import numpy as np
from damagerepair_counting.counting import counting

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

for label, group in females.groupby('Mouse Code'):
    if female_deaths.loc[female_deaths['Mouse Code'] == label, 'status'].values == 0:
        female_deaths.loc[female_deaths['Mouse Code'] == label, 'Age at death (months)'] = \
            group['Age (months)'].max()

# fix ages for drugs
males = males.loc[males['Age (months)'] <= 25]
females = females.loc[females['Age (months)'] <= 23]
males = males.loc[males['Age (months)'] <= 25]
females = females.loc[females['Age (months)'] <= 22]
males['treat'] = males['Mouse Group (1=drug, 3=control, 2=off drug)']
females['treat'] = females['Mouse Group (1=drug, 3=control)']

for data in [males, females]:
    keys = data['Mouse Code'].unique()
    deficits = list(data.columns.values[12:43])

print('Number missing before imputation:', data[deficits].isna().sum().sum()/(len(deficits)*data.shape[0]))
print('Percent missing before imputation:', data[deficits].isna().sum())
    
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

new_deficits = new_female_deficits


# fix death age units
for label, group in males.groupby(['Mouse Code']):
    if male_deaths.loc[male_deaths['Mouse Code'] == label, 'status'].item() == 1:
        
        male_deaths.loc[male_deaths['Mouse Code'] == label, 'Age at death (months)'] = \
            float(male_deaths.loc[male_deaths['Mouse Code'] == label, 'Age at death (days)'].item())*0.0328767 # convert to months

for label, group in females.groupby(['Mouse Code']):
    if female_deaths.loc[female_deaths['Mouse Code'] == label, 'status'].item() == 1:
        female_deaths.loc[female_deaths['Mouse Code'] == label, 'Age at death (months)'] = \
            float(female_deaths.loc[female_deaths['Mouse Code'] == label, 'Age at death (days)'].item())*0.0328767

male_data = counting(males, new_deficits, 'Mouse Code', 'Age (months)', male_deaths, 'Age at death (months)', 'status')
female_data = counting(females, new_deficits, 'Mouse Code', 'Age (months)', female_deaths, 'Age at death (months)', 'status')
data = pd.concat((male_data, female_data), axis=0)

data['mouse'] = data['Mouse Code']
data['n'] = data['deficit.count']
data['repair'] = data['repair.count']
data['damage'] = data['damage.count']
data['N'] = data['total.deficits']
data['f'] = data['FI']
data['death.age'] = data['Age at death (months)']
data['sex'] = data['Sex']
data['treatment'] = data['treat'].apply(lambda x: 'drug' if x < 3 else 'control')

deficits=new_deficits

data = data[['mouse', 'sex', 'treatment', 'baseline.age', 'time', 'f', 'repair', 'damage', 'n','N','delta.t','death.age', 'status'] +\
     ['d%d'%i for i in range(len(deficits))] +\
     ['d%d.r'%i for i in range(len(deficits))] +\
     ['d%d.d'%i for i in range(len(deficits))]]

for label, group in data.groupby('mouse'):
    if group['death.age'].max() - (group['time']+group['baseline.age']).max() < 0:
        data.loc[data['mouse']==label,'death.age'] = (group['time']+group['baseline.age']).max()
        
data.to_csv('../datasets/enalapril_data.csv', index=False)

surv_data = []
for label, group in data.groupby('mouse'):
     surv_data.append(group[['mouse','sex','treatment','baseline.age', 'death.age', 'status']].iloc[0])
surv_data = pd.DataFrame(surv_data)

surv_data.to_csv('../datasets/enalapril_surv_data.csv', index=False)
