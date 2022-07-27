import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['axes.titlesize'] = 7

enalapril_repair = pd.read_csv('../datasets/enalapril_repair_deficits.csv', index_col=0)
exercise_repair = pd.read_csv('../datasets/exercise_repair_deficits.csv', index_col=0)
schultz_repair = pd.read_csv('../datasets/schultz_repair_deficits.csv', index_col=0)
elsa_repair = pd.read_csv('../datasets/human_damage_deficits.csv', index_col=0)

#fig, ax = plt.subplots(1, 3, figsize=(8,4))
#fig, ax = plt.subplots(1, 3, figsize=(8,2.5))
fig, ax = plt.subplots(1, 4, figsize=(8,2.5))
ax=ax.flatten()

ax[0].title.set_text('l) Mouse dataset 1')


deficits = enalapril_repair['deficit']
deficits = np.array([d[:-1] for d in deficits], dtype=str)

deficits = deficits[enalapril_repair['diff.right'] < np.inf]
deficit, counts = np.unique(deficits, return_counts=True)

deficits_full = pd.read_csv('../../step_data/Enalapril_male.csv').columns.values[12:43]

deficit = list(deficit[np.argsort(counts)][::-1])
counts = list(counts[np.argsort(counts)][::-1])
for d in deficits_full:
    if d not in deficit:
        deficit.append(d)
        counts.append(0)


ax[0].bar(np.arange(len(counts)), counts, tick_label = deficit)
ax[0].set_xticklabels(deficit, rotation='vertical', fontsize=4)

output_1 = pd.DataFrame({'deficit': deficit, 'repair count': counts})
output_1.to_csv('../figure_data/figure5_supplement4/mouse1_repair_count.csv')

ax[0].set_yticks([0,50,100,150])
ax[0].set_yticklabels([0,50,100, 150], fontsize=6)

ax[1].title.set_text('m) Mouse dataset 2')


deficits = exercise_repair['deficit']
deficits = np.array([d[:-1] for d in deficits], dtype=str)

deficits = deficits[exercise_repair['diff.right'] < np.inf]
deficit, counts = np.unique(deficits, return_counts=True)

deficits_full = pd.read_csv('../../step_data/Exercise.csv').columns.values[12:43]

deficit = list(deficit[np.argsort(counts)][::-1])
counts = list(counts[np.argsort(counts)][::-1])
for d in deficits_full:
    if d not in deficit:
        deficit.append(d)
        counts.append(0)

ax[1].bar(np.arange(len(counts)), counts, tick_label = deficit)
ax[1].set_xticklabels(deficit, rotation='vertical', fontsize=4)

output_2 = pd.DataFrame({'deficit': deficit, 'repair count': counts})
output_2.to_csv('../figure_data/figure5_supplement4/mouse2_repair_count.csv')

ax[1].set_yticks([0,25,50,75])
ax[1].set_yticklabels([0,25,50,75], fontsize=6)


ax[2].title.set_text('n) Mouse dataset 3')


deficits = schultz_repair['deficit']
deficits = np.array([d[:-1] for d in deficits], dtype=str)

deficits = deficits[schultz_repair['diff.right'] < np.inf]
deficit, counts = np.unique(deficits, return_counts=True)

deficits_full = pd.read_csv('../../step_data/Schultz.csv').columns.values[1:30]

deficit = list(deficit[np.argsort(counts)][::-1])
counts = list(counts[np.argsort(counts)][::-1])
for d in deficits_full:
    if d not in deficit:
        deficit.append(d)
        counts.append(0)


ax[2].bar(np.arange(len(counts)), counts, tick_label = deficit)
ax[2].set_xticklabels(deficit, rotation='vertical', fontsize=4)


output_3 = pd.DataFrame({'deficit': deficit, 'repair count': counts})
output_3.to_csv('../figure_data/figure5_supplement4/mouse3_repair_count.csv')

ax[2].set_yticks([0,25,50,75])
ax[2].set_yticklabels([0,25,50,75], fontsize=6)

#ax[0].set_xlabel('Deficit')
ax[0].set_ylabel('Repair counts', fontsize=7)






ax[3].title.set_text('o) Human dataset')


deficits = elsa_repair['deficit']


#deficits = np.array([d[:-1] for d in deficits], dtype=str)


deficits = deficits[elsa_repair['diff.right'] < np.inf]
deficit, counts = np.unique(deficits, return_counts=True)

deficits_full = pd.read_csv('../clean_data/extracted_elsa.csv').columns.values[4:]


human_deficit_names = np.array(['Walk 100 yards', 'Sit ~2 hours', 'Up from chair',
                       'Several flights of stairs',
                       'One flight of stairs',
                       'Stooping/kneeling/crouching',
                       'Reaching/extending arms',
                       'Pulling/pushing large objects',
                       'Lifting/carrying over 10 lbs',
                       'Picking up a 5p coin',
                       'Dressing' ,
                       'Walking across a room' ,
                       'Bathing/showering' ,
                       'Eating/cutting food'  ,
                       'Getting in/out of bed' ,
                       'Using the toilet'  ,
                       'Using a map in a strange place'  ,
                       'Preparing a hot meal' ,
                       'Shopping for groceries' ,
                       'Making telephone calls' ,
                       'Taking medications'  ,
                       'Work around the house/garden' ,
                       'Managing money'])


deficit = list(deficit[np.argsort(counts)][::-1])
counts = list(counts[np.argsort(counts)][::-1])
for d in deficits_full:
    if d not in deficit:
        deficit.append(d)
        counts.append(0)

deficit = [int(d[1:]) for d in deficit]

ax[3].bar(np.arange(len(counts)), counts, tick_label = human_deficit_names[deficit])
ax[3].set_xticklabels(human_deficit_names[deficit], rotation='vertical', fontsize=4)
ax[3].set_yticks([0,250,500,750,1000])
ax[3].set_yticklabels([0,250,500,750,1000], fontsize=6)


output_human = pd.DataFrame({'deficit': human_deficit_names[deficit], 'repair count': counts})
output_human.to_csv('../figure_data/figure5_supplement4/human_repair_count.csv')


plt.tight_layout()
#plt.savefig('Mouse_repair_counts.pdf')
plt.subplots_adjust(wspace=0.2)
plt.savefig('repair_counts_all.pdf')
