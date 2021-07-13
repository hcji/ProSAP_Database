# -*- coding: utf-8 -*-
"""
Created on Tue May 25 15:24:11 2021

@author: hcji
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.metrics import pairwise_distances

files = os.listdir('separated')
dbn, name, peps, subs = [], [], [], []
for f in files:
    n = f.split('.')[0]
    p = 'separated/{}'.format(f)
    d = pd.read_csv(p)
    for i, s in enumerate(d['Subunits_UniProt_IDs']):
        if 'ComplexName' in list(d.columns):
            n = d['ComplexName'][i]
        else:
            n = ''
        p = s.split(',')
        subs.append(p)
        name.append(n)
        dbn.append(f.split('.')[0])
        for pp in p:
            if pp not in peps:
                peps.append(pp)
    print('{} is finished'.format(f), '\n')


vecs = []
for s in subs:
    v = np.zeros(len(peps))
    for p in s:
        i = peps.index(p)
        v[i] = 1
    vecs.append(v)
vecs = np.array(vecs)

dist_matrix = pairwise_distances(vecs, metric='jaccard', n_jobs=-1)
sim_matrix = 1 - dist_matrix

used = []
new_name, new_subs, new_dbn = [], [], []
for i in tqdm(range(len(peps))):
    if i in used:
        continue
    pp, db = [], []
    wh = np.where(sim_matrix[i,:] > 0.7)[0]
    wh = wh[wh >= i]
    for j in wh:
        pp += subs[j]
        db += [dbn[j]]
        used.append(j)
    pp = list(set(pp))
    db = list(set(db))
    new_subs.append(pp)
    new_dbn.append(db)
    new_name.append(name[i])

n = np.array([len(i) for i in new_subs])
k = np.where(n >= 2)[0]

res = []
for i, j in enumerate(k):
    n = 'Complex_{}'.format(i)
    m = new_name[j]
    s = new_subs[j]
    d = new_dbn[j]
    p = s[0]
    b = d[0]
    if len(d) > 1:
        for dd in d[1:]:
            b += ', ' + dd
    for pp in s[1:]:
        p += ', ' + pp
    res.append([n, m, p, b]) 
res = pd.DataFrame(res)
res.columns = ['Complex_ID', 'Complex_Name', 'Subunits_UniProt_IDs', 'Database']
res.to_csv('Consensus_Complex.csv', index=False)


