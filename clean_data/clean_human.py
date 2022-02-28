import numpy as np
import pandas as pd
import os
from pyreadstat import read_sav
from damagerepair_counting import counting

if not os.path.exists('../datasets/'):
    os.makedirs('../datasets/')

data = pd.read_csv('extracted_elsa.csv')
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

cols = ['id', 'wave', 'age', 'sex', 'n', 'wealth'] + deficits_n
data = data[cols]
data = data.groupby('id').filter(lambda x: x.shape[0] > 6)
deficits = deficits_n


data = counting(data, deficits, 'id', 'age')

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

data.to_csv('../datasets/human_data.csv')
