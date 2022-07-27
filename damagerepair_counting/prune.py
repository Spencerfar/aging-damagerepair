import numpy as np
import pandas as pd

# prunes damage/repair events if they only occur for one time step
def prune(data_raw, deficits, id_column, time_column):

    # check input
    assert isinstance(data_raw, pd.DataFrame), 'data_raw={} must be a dataframe'.format(data_raw)
    assert isinstance(id_column, str), 'id_column={} must be a string'.format(id_column)
    assert isinstance(time_column, str), 'time_column={} must be a string'.format(time_column)
    assert isinstance(deficits, list), 'deficits={} must be a list of strings'.format(deficits)
    for item in deficits:
        assert isinstance(item, str), 'deficits={} must be a list of strings'.format(deficits)
    
    data = data_raw.copy()

    pruned_repair = []
    pruned_damage = []
    
    for label, group in data_raw.groupby([id_column]):

        if len(group) > 5:

            for d in deficits:
                
                updated_deficit = group[d].values
                for t_itr in range(2, len(group)-2):
                    
                    if group[d].iloc[t_itr-1] < group[d].iloc[t_itr] and group[d].iloc[t_itr] > group[d].iloc[t_itr+1] and group[d].iloc[t_itr-1] == group[d].iloc[t_itr+1] and group[d].iloc[t_itr-2] == group[d].iloc[t_itr+2] and group[d].iloc[t_itr-2] == group[d].iloc[t_itr-1] and group[d].iloc[t_itr+2] == group[d].iloc[t_itr+1]:
                        
                        updated_deficit[t_itr] = group[d].iloc[t_itr-1]
                        pruned_damage.append(d)
                        
                        
                    elif group[d].iloc[t_itr-1] > group[d].iloc[t_itr] and group[d].iloc[t_itr] < group[d].iloc[t_itr+1] and group[d].iloc[t_itr-1] == group[d].iloc[t_itr+1] and group[d].iloc[t_itr-2] == group[d].iloc[t_itr+2] and group[d].iloc[t_itr-2] == group[d].iloc[t_itr-1] and group[d].iloc[t_itr+2] == group[d].iloc[t_itr+1]:

                        updated_deficit[t_itr] = group[d].iloc[t_itr-1]
                        pruned_repair.append(d)
                data.loc[data[id_column] == label, d] = updated_deficit
    return data, pruned_repair, pruned_damage
