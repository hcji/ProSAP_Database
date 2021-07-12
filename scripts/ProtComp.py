# -*- coding: utf-8 -*-
"""
Created on Mon May 24 15:51:17 2021

@author: hcji
"""

import json
import numpy as np
import pandas as pd

# CORUM
with open('raw/CORUM_Core.json') as js:
    data = json.load(js)

data = pd.DataFrame(data)
data = data[['ComplexID', 'ComplexName', 'Organism', 'subunits(UniProt IDs)', 'GO description']]
data.columns = ['ComplexID', 'ComplexName', 'Organism', 'Subunits_UniProt_IDs', 'GO description']

for i in range(data.shape[0]):
    s = data.iloc[i,3]
    s = s.replace(';', ',')
    data.iloc[i,3] = s

data.to_csv('separated/CORUM_Core.csv', index = False)


with open('raw/CORUM_All.json') as js:
    data = json.load(js)

data = pd.DataFrame(data)
data = data[['ComplexID', 'ComplexName', 'Organism', 'subunits(UniProt IDs)', 'GO description']]
data.columns = ['ComplexID', 'ComplexName', 'Organism', 'Subunits_UniProt_IDs', 'GO description']

for i in range(data.shape[0]):
    s = data.iloc[i,3]
    s = s.replace(';', ',')
    data.iloc[i,3] = s

data.to_csv('separated/CORUM_All.csv', index = False)


# Humap2
data = pd.read_csv('raw/humap2_complexes_20200809.txt')
data = data[['HuMAP2_ID', 'Uniprot_ACCs']]
data.columns = ['ComplexID', 'Subunits_UniProt_IDs']

for i in range(data.shape[0]):
    s = data.iloc[i,1]
    s = s.replace(' ', ',')
    data.iloc[i,1] = s

data.to_csv('separated/Humap2.csv', index = False)


# Havugimana
data = pd.read_csv('raw/Havugimana_2012.csv')
data = data[['Predicted complex', 'Co-complex membership_UniProtKB AC']]
data.columns = ['ComplexID', 'Subunits_UniProt_IDs']
data.to_csv('separated/Havugimana_2012.csv', index = False)


# Kristensen
data = pd.read_csv('raw/Kristensen_2012.csv')
res = []
for i in range(1, max(data['Complex no'])+1):
    idx = 'Complex_' + str(i)
    wh = data.loc[data.loc[:,'Complex no'] == i,:]
    su = str(wh.iloc[0,4])
    for j in range(1, wh.shape[0]):
        su += ',' + str(wh.iloc[j,4])
    res.append([idx, su])
res = pd.DataFrame(res)
res.columns = ['ComplexID', 'Subunits_UniProt_IDs']
res.to_csv('separated/Kristensen_2012.csv', index = False)


# Heusel
data = pd.read_csv('raw/Heusel_BioPlex_2019.csv')
data = pd.DataFrame(data['subunits_coeluting'].drop_duplicates())
data.columns = ['Subunits_UniProt_IDs']
for i in range(data.shape[0]):
    s = data.iloc[i, 0]
    data.iloc[i, 0] = s.replace(';', ',')
idx = ['Complex_{}'.format(i) for i in range(len(data))]
data['ComplexID'] = idx
data = data[['ComplexID', 'Subunits_UniProt_IDs']]
data.to_csv('separated/Heusel_BioPlex_2019.csv', index = False)

data = pd.read_csv('raw/Heusel_StringDB_2019.csv')
data = pd.DataFrame(data['subunits_coeluting'].drop_duplicates())
data.columns = ['Subunits_UniProt_IDs']
for i in range(data.shape[0]):
    s = data.iloc[i, 0]
    data.iloc[i, 0] = s.replace(';', ',')
idx = ['Complex_{}'.format(i) for i in range(len(data))]
data['ComplexID'] = idx
data = data[['ComplexID', 'Subunits_UniProt_IDs']]
data.to_csv('separated/Heusel_StringDB_2019.csv', index = False)

data = pd.read_csv('raw/Heusel_CSB_2019.csv')
data = pd.DataFrame(data['subunits_coeluting'].drop_duplicates())
data.columns = ['Subunits_UniProt_IDs']
for i in range(data.shape[0]):
    s = data.iloc[i, 0]
    data.iloc[i, 0] = s.replace(';', ',')
idx = ['Complex_{}'.format(i) for i in range(len(data))]
data['ComplexID'] = idx
data = data[['ComplexID', 'Subunits_UniProt_IDs']]
data.to_csv('separated/Heusel_CSB_2019.csv', index = False)


# PCprophet
data = pd.read_excel('raw/PCprophet.xlsx')
data = data[['ComplexID', 'Members']]
data = data.drop_duplicates(ignore_index=True)
data.columns = ['ComplexName', 'Subunits_UniProt_IDs']
data['ComplexID'] = ['Complex_{}'.format(i) for i in range(len(data))]
data = data[['ComplexID', 'ComplexName', 'Subunits_UniProt_IDs']]
data.to_csv('separated/PCprophet.csv', index = False)
