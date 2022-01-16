from pyreadstat import read_sav
import numpy as np
import pandas as pd

#Read files. Column names vary wildly per file, so these need to be selected manually.
folder = "~/Downloads/UKDA-5050-spss/spss/spss24/"

files = ["wave_1_core_data_v3.sav", "wave_2_core_data_v4.sav", "wave_2_nurse_data_v2.sav",
        "wave_3_elsa_data_v4.sav", "wave_4_elsa_data_v3.sav", "wave_4_nurse_data.sav",
        "wave_5_elsa_data_v4.sav", 
         "wave_6_elsa_data_v2.sav", 
         "wave_6_elsa_nurse_data_v2.sav","wave_7_elsa_data.sav", 
         "wave_8_elsa_data_eul_v2.sav", "wave_8_elsa_nurse_data_eul_v1.sav",
        "wave_9_elsa_data_eul_v1.sav"]

cols = []

# wave1
cols.append(['idauniq','dhager','dhsex','dhdobyr'])

# wave2
cols.append(['idauniq','indsex','dhdobyr'])

# wave2 nurse
cols.append(['idauniq','confage'])

# wave3 
cols.append(['idauniq','dhager','dhsex','dhdobyr'])

# wave4
cols.append(['idauniq','dhager','dhsex'])

# wave4 nurse
cols.append(['idauniq','confage','dobyear'])

# wave5
cols.append(['idauniq','indager','dhsex','indobyr'])

# wave6
cols.append(['idauniq','indager','DhSex','Indobyr'])

# wave6 nurse
cols.append(['idauniq'])

# wave7
cols.append(['idauniq','indager','DhSex','Indobyr'])

# wave8
cols.append(['idauniq','indager','indobyr'])

# wave8 nurse
cols.append(['idauniq','indsex'])

# wave9
cols.append(['idauniq','indager','indsex','indobyr'])

keys = ['wave1', 'wave2', 'wave2_nurse',
       'wave3','wave4', 'wave4_nurse', 'wave5', 'wave6', 'wave6_nurse', 'wave7',
       'wave8','wave8_nurse', 'wave9']

waves = dict()
for i,f in enumerate(files):
    waves[keys[i]] = read_sav(folder+f, usecols=cols[i])[0]

waves_FI = dict()
for i in [1,2,3,4,6,7,8,9]:
    waves_FI['wave'+str(i)] = pd.read_csv('extracted_elsa_deficits_wave'+str(i)+'.csv')
    waves_FI['wave'+str(i)].fillna(-1.0, inplace=True)
    
eol2 = read_sav(folder+"elsa_eol_w2_archive_v1.sav", usecols=['idauniq','EiDateY'])[0]
eol2['idauniq'] = eol2['idauniq'].astype(int)
eol2 = eol2.set_index('idauniq')

eol3 = read_sav(folder+"elsa_eol_w3_archive_v1.sav",usecols=['idauniq','EiDateY'])[0]
eol3['idauniq'] = eol3['idauniq'].astype(int)
eol3 = eol3.set_index('idauniq')

eol4 = read_sav(folder+"elsa_eol_w4_archive_v1.sav",usecols=['idauniq','EiDateY'])[0]
eol4['idauniq'] = eol4['idauniq'].astype(int)
eol4 = eol4.set_index('idauniq')

eol6 = read_sav(folder+"elsa_endoflife_w6archive.sav",usecols=['idauniq','EiDateY'])[0]
eol6['idauniq'] = eol6['idauniq'].astype(int)
eol6 = eol6.set_index('idauniq')


for key in keys:
    waves[key]['idauniq'] = waves[key]['idauniq'].astype(int)
    waves[key] = waves[key].set_index('idauniq')

for key in waves_FI.keys():
    waves_FI[key]['idauniq'] = waves_FI[key]['idauniq'].astype(int)
    waves_FI[key] = waves_FI[key].set_index('idauniq')

combined_waves = dict()

combined_waves['wave1'] = pd.concat([waves['wave1'], 
                                     waves_FI['wave1']], 
                                    axis = 1, sort=False)

combined_waves['wave2'] = pd.concat([waves['wave2_nurse'], 
                                     waves['wave2'], 
                                     waves_FI['wave2'],
                                     eol2], axis=1, sort=False) #concat when index same but columns diff

combined_waves['wave3'] = pd.concat([waves['wave3'], 
                                     waves_FI['wave3'], 
                                     eol3], axis=1, sort=False)

combined_waves['wave4'] = pd.concat([waves['wave4_nurse'], 
                                     waves['wave4'], 
                                     waves_FI['wave4'], 
                                     eol4], axis=1, sort=False)

combined_waves['wave5'] = pd.concat([waves['wave5']],
                                    axis = 1, sort=False) #waves_FI['wave5']

combined_waves['wave6'] = pd.concat([waves['wave6_nurse'], 
                                     waves['wave6'], 
                                     waves_FI['wave6'], 
                                     eol6], axis=1, sort=False)

combined_waves['wave7'] = pd.concat([waves['wave7'], 
                                     waves_FI['wave7']],
                                    axis = 1, sort=False)

combined_waves['wave8'] = pd.concat([waves['wave8_nurse'], 
                                     waves['wave8'], 
                                     waves_FI['wave8']], 
                                    axis=1, sort=False)

combined_waves['wave8'] = pd.concat([waves['wave8_nurse'], 
                                     waves['wave8'], 
                                     waves_FI['wave8']], 
                                    axis=1, sort=False)

combined_waves['wave9'] = pd.concat([waves['wave9'],  
                                     waves_FI['wave9']], 
                                    axis=1, sort=False)

combined_waves['wave1'].rename(columns={'dhager':'age', 'dhsex':'sex','dhdobyr':'dob'}, inplace=True)

combined_waves['wave2'].rename(columns={'confage':'age','indsex':'sex','dhdobyr':'dob','EiDateY':'dod'}, inplace=True)

combined_waves['wave3'].rename(columns={'dhager':'age','dhsex':'sex','dhdobyr':'dob','EiDateY':'dod'}, inplace=True)

combined_waves['wave4'].rename(columns={'confage':'age', 'dhsex':'sex','dobyear':'dob','EiDateY':'dod'}, inplace=True)

combined_waves['wave5'].rename(columns={'indager':'age','dhsex':'sex','indobyr':'dob','EiDateY':'dod'}, inplace=True)

combined_waves['wave6'].rename(columns={'indager':'age', 'DhSex':'sex','Sex':'sex','Indobyr':'dob','EiDateY':'dod'}, inplace=True)

combined_waves['wave7'].rename(columns={'indager':'age','DhSex':'sex','Indobyr':'dob'}, inplace=True)

combined_waves['wave8'].rename(columns={'indager':'age', 'indsex':'sex','indobyr':'dob'}, inplace=True)

combined_waves['wave9'].rename(columns={'indager':'age', 'indsex':'sex','indobyr':'dob'}, inplace=True)

data = pd.concat(list(combined_waves.values()), keys=[1,2,3,4,5,6,7,8,9], sort=False)


data['dod'].fillna(-1,inplace=True)
data['dod'] = data['dod'].astype(int)

data['age'].replace(99.0,90,inplace=True) #collapsed at 90
data['age'].replace(90.0,-1.0,inplace=True) #collapsed at 90
data['age'].replace(-9.0,-1.0,inplace=True) #refusal
data['age'].replace(-8.0,-1.0,inplace=True) #don't know
data['age'].replace(-7.0,-1.0,inplace=True) #c
data['age'].fillna(-1,inplace=True)
data['age'] = data['age'].astype(int)

data['sex'].fillna(-1, inplace=True)
data['sex'].replace(1.0, 0.0, inplace=True)
data['sex'].replace(2.0, 1.0, inplace=True)

# impute sex if some are missing
data['sex'].replace(-1.0,np.nan,inplace=True)
data['sex'] = data.groupby('idauniq')['sex'].transform(lambda x: x.fillna(method='ffill'))




# Fill missing age values by using the previous age and the fact that waves are ~2 years apart.
for index, group in data.groupby('idauniq'):
    if np.any(np.diff(group['age']) <= 0) or np.any((group['age']) <= 0):
        waves = group.xs(index,level=1).index.values
        for w, wave in enumerate(waves):
            if w > 0:
                #print(wave, index, data.loc[(waves[w-1],index),'age'])
                if data.loc[(wave,index),'age'] <= data.loc[(waves[w-1],index),'age'] and data.loc[(wave,index),'age'] > 0 and data.loc[(waves[w-1],index),'age'] > 0:
                    data.loc[(wave,index),'age'] = data.loc[(waves[w-1],index),'age'] + 2*(waves[w] - waves[w-1])

                if data.loc[(waves[w-1],index),'age'] > 0 and data.loc[(wave,index),'age'] < 0:
                    data.loc[(wave,index),'age'] = data.loc[(waves[w-1],index),'age'] + 2*(waves[w] - waves[w-1])


for index, group in data.groupby('idauniq'):
    if np.any(np.diff(group['age']) <= 0) or np.any((group['age']) <= 0):
        waves = group.xs(index,level=1).index.values[::-1]
        for w, wave in enumerate(waves):
            if w > 0:
                
                if data.loc[(wave,index),'age'] >= data.loc[(waves[w-1],index),'age'] and data.loc[(wave,index),'age'] > 0 and data.loc[(waves[w-1],index),'age'] > 0:
                    data.loc[(wave,index),'age'] = data.loc[(waves[w-1],index),'age'] - 2*np.abs(waves[w] - waves[w-1])

                if data.loc[(waves[w-1],index),'age'] > 0 and data.loc[(wave,index),'age'] < 0:
                    data.loc[(wave,index),'age'] = data.loc[(waves[w-1],index),'age'] - 2*np.abs(waves[w] - waves[w-1])


for index, group in data.groupby('idauniq'):
    if np.any(group['age'] < 0):
        waves = group.xs(index,level=1).index.values
        for w, wave in enumerate(waves):
            if w == 0 and len(waves) > 1 and data.loc[(wave,index),'age'] < 0 and data.loc[(waves[w+1],index),'age'] > 0:
                data.loc[(wave, index), 'age'] = data.loc[(waves[w+1],index),'age'] - 2*np.abs(waves[w] - waves[w-1])
                
data['sex'] = data['sex'].fillna(-1)
deficits = ['d%d'%i for i in range(23)]
data[deficits] = data[deficits].fillna(-1.0)

# Remove all individuals who don't have at least 1 measurement, have sex missing, or have missing ages that couldn't be approximated through the year of interview.
indexes = [x[1] for x in data.index.values]
dropping = []
for i,index in enumerate(np.unique(indexes)):
    
    selected = data.xs(index,level=1)
    
    remove = True
    for w,wave in enumerate(selected.index.values):
        current = data.loc[(wave,index),deficits]
        
        if np.any(current.values > -1.0) and np.all(selected['age'].values > 0) \
        and np.all(selected['sex'].values >= 0):
            remove = False
            break
            
    if remove:
        dropping.append(index)
data.drop(dropping,level=1,inplace=True)

data = data.reset_index()
data.rename(columns={'level_0':'wave','idauniq':'id'}, inplace=True)

columns = ['id','wave','age','sex'] + deficits
data.fillna(-1.0, inplace=True)

transition_data = data[columns].sort_values(by=['id','wave'])
transition_data = transition_data.dropna(subset=['age'], how = 'any')

transition_data[columns].to_csv('extracted_elsa.csv',index=False)
