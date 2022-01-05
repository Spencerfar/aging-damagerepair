import pandas as pd
import numpy as np
from damagerepair_counting.counting import counting

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


data.to_csv('../datasets/schultz_data.csv', index=False)

surv_data = []
for label, group in data.groupby('mouse'):
     surv_data.append(group[['mouse','sex','baseline.age', 'death.age', 'status']].iloc[0])
surv_data = pd.DataFrame(surv_data)

surv_data.to_csv('../datasets/schultz_surv_data.csv', index=False)
