#%% Classification
import json
infile = open("TCGA Data/TCGA-KICH/files.2018-06-04.json", "r")
file_list = json.loads(infile.read())
infile.close()

infile = open("TCGA Data/TCGA-KICH/clinical.project-TCGA-KICH.2018-06-03.json", "r")
case_list = json.loads(infile.read())
infile.close()

id_map = {}

for file in file_list:
    if file['file_name'] not in id_map:
        if len(file['cases']) != 1:
            print("Error:",file['file_name'],"has more than one case id!")
            exit(-1)
        else:
            id_map[file['file_name']] = file['cases'][0]['case_id']
    else:
        print("Error:",file['file_name'],"not unique!")
        exit(-1)

filtered_case_list = []

for case in case_list:
    if case['diagnoses'][0]['days_to_last_follow_up'] != None:
        filtered_case_list.append(case['case_id'])


from pathlib import Path
import gzip

expr_map = {}

replicate_count = 0
path_list = Path("TCGA Data/TCGA-KICH").glob('**/*.FPKM.txt.gz')
for path in path_list:
    if id_map[str(path).split('\\')[-1]] not in expr_map and id_map[str(path).split('\\')[-1]] in filtered_case_list:
        expr_map[id_map[str(path).split('\\')[-1]]] = {}
    else:
        continue
    
    infile = gzip.open(str(path),"rt")
    for line in infile:
        gene_id = line.split('\t')[0]
        expr_level = float(line.split('\t')[1])
        if gene_id not in expr_map[id_map[str(path).split('\\')[-1]]]:
            expr_map[id_map[str(path).split('\\')[-1]]][gene_id] = []
        expr_map[id_map[str(path).split('\\')[-1]]][gene_id].append(expr_level)
    infile.close()

for case_id in expr_map:
    for gene_id in expr_map[case_id]:
        expr_map[case_id][gene_id] = sum(expr_map[case_id][gene_id])/len(expr_map[case_id][gene_id])

import pandas as pd
df_expr = pd.DataFrame.from_dict(expr_map)
print(df_expr.shape[0],df_expr.shape[1])
df_expr.to_csv("TCGA Data/TCGA-KICH/ExpressionTable.csv")

cluster_map = {}

median_array = df_expr.median(axis=1)
median_map = {}
for i in range(df_expr.shape[0]):
    median_map[df_expr.index[i]] = median_array[i]

for case_id in expr_map:
    cluster_map[case_id] = {}
    num0 = 0
    num1 = 0
    for gene_id in expr_map[case_id]:
        if abs(expr_map[case_id][gene_id]-median_array[gene_id]) < 0.0000000001:
            cluster_map[case_id][gene_id] = -1
        elif expr_map[case_id][gene_id] > median_array[gene_id]:
            cluster_map[case_id][gene_id] = 1
        else:
            cluster_map[case_id][gene_id] = 0

df_cluster = pd.DataFrame.from_dict(cluster_map)         
df_cluster.to_csv("TCGA Data/TCGA-KICH/ClusterTable.csv")

#%% Survival Analysis

df_cluster = pd.read_csv("TCGA Data/TCGA-KIRP/ClusterTable.csv", index_col=0)

survival_map = {'T':{}, 'E':{}, 'group': {}}

for case in case_list:
    if case['case_id'] not in df_cluster.columns:
        continue
    if case['diagnoses'][0]['days_to_last_follow_up'] != None:
        survival_map['T'][case['case_id']] = float(case['diagnoses'][0]['days_to_last_follow_up'])
    else:
        continue
    
    if case['diagnoses'][0]['vital_status'] == 'alive':
        survival_map['E'][case['case_id']] = 1
    elif case['diagnoses'][0]['vital_status'] == 'dead':
        survival_map['E'][case['case_id']] = 0
    else:
        print("Error: No clear status indicated!")
        exit(-1)

if len(survival_map['T']) != len(survival_map['E']):
    print("Error: Wrong dimemsion!")

print(len(survival_map['T']))


from lifelines.statistics import logrank_test

max_p = []

for gene_id in df_cluster.index:
    for case_id in df_cluster.columns:
        if case_id not in survival_map['T'] and case_id not in survival_map['E']:
            continue
        if df_cluster[case_id][gene_id] == 1:
            survival_map['group'][case_id] = "High Expression"
        if df_cluster[case_id][gene_id] == 0:
            survival_map['group'][case_id] = "Low Expression"
            
    df = pd.DataFrame.from_dict(survival_map)
    T = df['T']
    E = df['E']
    groups = df['group']
    ix = (groups == 'High Expression')
    
    if len(T[~ix]) != len(T[ix]):
        continue
    
    results = logrank_test(T[~ix], T[ix], event_observed_A=E[~ix], event_observed_B=E[~ix])
    max_p.append((gene_id, results.p_value))

sorted_max_p = sorted(max_p, key = lambda x:x[1])
print(sorted_max_p[0:6])    

#%% K-M Plot
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()

for gene_id in [gene for (gene, p_value) in sorted_max_p[0:6]]:
    for case_id in df_cluster.columns:
        if case_id not in survival_map['T'] and case_id not in survival_map['E']:
            continue
        if df_cluster[case_id][gene_id] == 1:
            survival_map['group'][case_id] = "High Expression"
        if df_cluster[case_id][gene_id] == 0:
            survival_map['group'][case_id] = "Low Expression"
            
    df = pd.DataFrame.from_dict(survival_map)
    T = df['T']
    E = df['E']
    groups = df['group']
    ix = (groups == 'High Expression')
    
    
    results = logrank_test(T[~ix], T[ix], event_observed_A=E[~ix], event_observed_B=E[~ix])
    
    kmf.fit(T[~ix], E[~ix], label='Low Expression')
    ax = kmf.plot()
    ax.set_title(gene_id)
    ax.text(0,0,"p-value:" + str(results.p_value))
    kmf.fit(T[ix], E[ix], label='High Expression')
    kmf.plot(ax=ax)

#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%

