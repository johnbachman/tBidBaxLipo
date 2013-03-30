from tbidbaxlipo.data.gridv1 import *
import pandas as pd
from tbidbaxlipo.util.numsort import sorted_copy as sort_numeric
import numpy as np

flattened = []

lipo_concs_str = sort_numeric(data_tbidmaj.keys())
tbid_concs_str = sort_numeric(data_tbidmaj[lipo_concs_str[0]].keys())
bax_concs_str = sort_numeric(data_tbidmaj[lipo_concs_str[0]][tbid_concs_str[0]].keys())

tuples = []

for lipo_conc_str in lipo_concs_str:
    for bax_conc_str in bax_concs_str:
        for tbid_conc_str in tbid_concs_str:
            tuples.append(tuple([float(lipo_conc_str), float(tbid_conc_str), float(bax_conc_str)]))
            flattened.append(data_tbidmaj[lipo_conc_str][tbid_conc_str][bax_conc_str])

np_flattened = np.array(flattened)
index = pd.MultiIndex.from_tuples(tuples, names=['Liposomes','tBid','Bax'])

df = pd.DataFrame(np_flattened, index=index, columns=time)
dt = df.T

"""
for lipo_conc_str in lipo_concs_str:
    for bax_conc_str in bax_concs_str:
        for tbid_conc_str in tbid_concs_str:
            for time_sec, value in zip(time, data_tbidmaj[lipo_conc_str][tbid_conc_str][bax_conc_str]):
                flattened.append([lipo_conc_str, bax_conc_str, tbid_conc_str, time_sec, value])

index = pd.MultiIndex.from_tuples([tuple(row[0:4]) for row in flattened],
                                  names=['lipos','bax','tbid','time'])
values = [row[4] for row in flattened]
s = pd.Series(values, index=index)

p = pd.Panel(data_tbidmaj)
p = p.reindex_axis(sort_numeric(p.items), axis=0)
p = p.reindex_axis(sort_numeric(p.major_axis), axis=1)
p = p.reindex_axis(sort_numeric(p.minor_axis), axis=2)
p.items.name = 'lipos'
p.major_axis.name = 'bax'
p.minor_axis.name = 'tbid'
df = p.to_frame()
"""
