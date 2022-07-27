import numpy as np
import pandas as pd

def counting(data_raw, deficits, id_column, time_column, deaths = None, death_column = None, status_column = None):

    # check input
    assert isinstance(data_raw, pd.DataFrame), 'data_raw={} must be a dataframe'.format(data_raw)
    assert isinstance(id_column, str), 'id_column={} must be a string'.format(id_column)
    assert isinstance(time_column, str), 'time_column={} must be a string'.format(time_column)
    assert isinstance(deficits, list), 'deficits={} must be a list of strings'.format(deficits)
    for item in deficits:
        assert isinstance(item, str), 'deficits={} must be a list of strings'.format(deficits)
    if deaths is not None:
        assert isinstance(deaths, pd.DataFrame), 'deaths={} must be a dataframe'.format(deaths)
        assert isinstance(death_column, str), 'death_column={} must be a string if a deaths dataframe is passed'.format(death_column)
        assert isinstance(status_column, str), 'status_column={} must be a string if a deaths dataframe is passed'.format(status_column)
    
    data = data_raw.copy(deep=True)

    
    # add new columns to dataframe
    new_cols = ['time', 'baseline.age', 'repair.count', 'damage.count', 'delta.t', 'deficit.count', 'total.deficits']
    
    if deaths is not None:
        new_cols = new_cols + ['status', 'death.age']
    
    new_deficits = ['d%d'%i for i in range(len(deficits))]
    if any((True for x in new_deficits if x in new_cols)):
        new_cols = new_cols + new_deficits
    
    new_repair_deficits = ['d%d.r'%i for i in range(len(deficits))]
    new_cols = new_cols + new_repair_deficits
    
    new_damage_deficits = ['d%d.d'%i for i in range(len(deficits))]
    new_cols = new_cols + new_damage_deficits

    new_cols_df = pd.DataFrame(np.ones((data.shape[0], len(new_cols)))*np.nan, columns = new_cols, index=data.index)
    data = pd.concat((data, new_cols_df), axis=1)
    
    
    # go through each individual subject
    for label, group in data_raw.groupby([id_column]):
    
        time = group[time_column] - group[time_column].iloc[0]
        
        data.loc[data[id_column] == label, 'time'] = time
        data.loc[data[id_column] == label, 'baseline.age'] = group[time_column].iloc[0]

        if deaths is not None:
            
            data.loc[data[id_column] == label, status_column] = deaths.loc[deaths[id_column] == label, status_column].item()
            
            if deaths.loc[deaths[id_column] == label, status_column].item() == 1:
                data.loc[data[id_column] == label, death_column] = \
                    deaths.loc[deaths[id_column] == label, death_column].item()
            else:
                data.loc[data[id_column] == label, death_column] = \
                    deaths.loc[deaths[id_column] == label, death_column].item()
        
        N = (~group[deficits].isna()).sum(axis=1)
        n = group[deficits].sum(axis=1)
        
        data.loc[data[id_column] == label, 'deficit.count'] = n
        data.loc[data[id_column] == label, 'total.deficits'] = N
        data.loc[data[id_column] == label, 'FI'] = (n/N)
        
        deficit_values = group[deficits].values
        
        damaged_all = ((group[deficits].iloc[:-1] == 0).values & \
                   (group[deficits].iloc[1:] == 1).values)
        damaged = damaged_all.sum(axis=1)
        
        repaired_all = ((group[deficits].iloc[:-1] == 1).values & \
                    (group[deficits].iloc[1:] == 0).values)
        repaired = repaired_all.sum(axis=1)
        
        data.loc[data[id_column] == label, 'delta.t'] = np.append(np.diff(group[time_column]), np.nan)
        
        damaged = np.append(damaged, np.nan)
        repaired = np.append(repaired, np.nan)
        
        data.loc[data[id_column] == label, 'damage.count'] = damaged
        data.loc[data[id_column] == label, 'repair.count'] = repaired
        
        for di,d in enumerate(['d%d'%i for i in range(len(deficits))]):
            data.loc[data[id_column] == label, d] = deficit_values[:,di]
        
        for di,d in enumerate(['d%d.r'%i for i in range(len(deficits))]):
            data.loc[data[id_column] == label, d] = np.append(repaired_all[:,di],np.nan)
        
        for di,d in enumerate(['d%d.d'%i for i in range(len(deficits))]):
            data.loc[data[id_column] == label, d] = np.append(damaged_all[:,di],np.nan)
        
    return data
