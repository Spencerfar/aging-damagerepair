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
    
    data = data_raw.copy()
    
    data['time'] = np.nan
    data['baseline.age'] = np.nan
    data['repair.count'] = np.nan
    data['damage.count'] = np.nan
    data['delta.t'] = np.nan
    data['deficit.count'] = np.nan
    data['total.deficits'] = np.nan

    if deaths is not None:
        data['status'] = np.nan
        data['death.age'] = np.nan

    if deficits != ['d%d'%i for i in range(len(deficits))]:
        for d in ['d%d'%i for i in range(len(deficits))]:
            data[d] = np.nan
    
    for d in ['d%d.r'%i for i in range(len(deficits))]:
        data[d] = np.nan
    
    for d in ['d%d.d'%i for i in range(len(deficits))]:
        data[d] = np.nan

    for label, group in data.groupby([id_column]):
    
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
