# -*- coding: utf-8 -*-
"""
Created on Mon May 24 13:48:47 2021

@author: hcji
"""


import numpy as np
import pandas as pd

from tqdm import tqdm

BioPlex_293 = pd.read_csv('db/raw/BioPlex_293T_Network_10K_Dec_2019.tsv', sep='\t')
BioPlex_HCT = pd.read_csv('db/raw/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv', sep='\t')

# BioPlex 293T
pA, pB = [], []
for i in tqdm(range(BioPlex_293.shape[0])):
    a = BioPlex_293.iloc[i, 2].split('-')[0]
    b = BioPlex_293.iloc[i, 3].split('-')[0]
    ch = 0
    if a in pA:
        w = np.where(np.array(pA) == a)[0]
        for j in w:
            if pB[j] == b:
                ch = 1
                continue
    if ch == 0:
        pA.append(a)
        pB.append(b)
    
pp = pd.DataFrame([pA, pB])
pp = pp.T
pp.columns = ['Protein A', 'Protein B']
pp.to_csv('db/BioPlex_293T.csv', index = False)


# bioPlex HCT116
pA, pB = [], []
for i in tqdm(range(BioPlex_HCT.shape[0])):
    a = BioPlex_HCT.iloc[i, 2].split('-')[0]
    b = BioPlex_HCT.iloc[i, 3].split('-')[0]
    ch = 0
    if a in pA:
        w = np.where(np.array(pA) == a)[0]
        for j in w:
            if pB[j] == b:
                ch = 1
                continue
    if ch == 0:
        pA.append(a)
        pB.append(b)
    
pp = pd.DataFrame([pA, pB])
pp = pp.T
pp.columns = ['Protein A', 'Protein B']
pp.to_csv('db/BioPlex_HCT116.csv', index = False)


## bioPlex Combine
BioPlex_all = pd.concat([BioPlex_HCT, BioPlex_293])

pA, pB = [], []
for i in tqdm(range(BioPlex_all.shape[0])):
    a = BioPlex_all.iloc[i, 2].split('-')[0]
    b = BioPlex_all.iloc[i, 3].split('-')[0]
    ch = 0
    if a in pA:
        w = np.where(np.array(pA) == a)[0]
        for j in w:
            if pB[j] == b:
                ch = 1
                continue
    if ch == 0:
        pA.append(a)
        pB.append(b)
    
pp = pd.DataFrame([pA, pB])
pp = pp.T
pp.columns = ['Protein A', 'Protein B']
pp.to_csv('db/BioPlex_all.csv', index = False)

