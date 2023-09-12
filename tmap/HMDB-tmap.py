"""
tmap for the HMDB dataset.
Tmap 
1. LSH Forest indexing 
2. Build c-approximate k-NN
3. Calculate MST of the k-NN
4. Lay tree on Euclidian plane using OGDF modular C++ library
"""


#import packages
import os
import tmap as tm

import faerun as Faerun
import pandas as pd
import pickle

from faerun import host
from map4 import MAP4Calculator
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from timeit import default_timer as timer


## Creation of tmap
lf = tm.LSHForest(1024, 128)

directory = os.getcwd()
fragments = pd.read_pickle('HMDB_prep_tmap.pkl')
tmap_plot = os.path.join(directory, 'tmap')
tmap_coordinates = os.path.join(directory, 'TMAP_coordinates.dat')
tmap_topology = os.path.join(directory, 'TMAP_topology.dat')

print('load data', flush=True)
print(fragments.shape, flush=True)
print('\n')

# Encode SMILES into Map4
dim = 1024
MAP4 = MAP4Calculator(dimensions=dim)
ENC = tm.Minhash(dim)
"""
MAP4 performs much better than MHFP6 for mapping the Human Metabolome Database (HMDB)
MHFP6 fails to properly distinguish between related metabolites and the map consists 
of very large groups of molecules appearing as “grapes"
"""
start = timer()
fps = []
dsmiles = fragments['smiles_desalted'].copy()


# Chunks to ease computational time 

def chunk_df(df, chunksize=1000):
    chunks_ls = list()
    num_chunks = len(df) // chunksize + (1 if len(df) % chunksize else 0)
    for i in range(num_chunks):
        chunks_ls.append(df[i * chunksize: (i + 1) * chunksize])
    return chunks_ls

chunks2 = chunk_df(dsmiles, chunksize=1000)

x = 0   
for chunk in chunks2:
    x += 1
    print(f'calculating MAP4 for chunk {x}')
    map4fps = []
    for index, i in chunk.items():
       mol = Chem.MolFromSmiles(i) # transform smiles to mol 
       mol_map4 = MAP4.calculate(mol) 
       map4fps.append(mol_map4)
    fps.extend(map4fps)
    print('chunk done...')


print(f'Encoding smiles took: {(timer() - start) * 0.01666} minutes.')

start = timer()

#LSH Forest indexing
print('start LSH Forest indexing', flush=True)
lf.batch_add(fps)
lf.index()

#configuration of k-NN graph layout
cfg = tm.LayoutConfiguration()
cfg.k = 50                             #100 in tmap drugbank
cfg.kc = 50
cfg.mmm_repeats = 2
cfg.node_size = 1/70                   #decreasing resolves overlap in very crowded tree
cfg.sl_extra_scaling_steps = 10
cfg.sl_repeats = 2

#creation of k-NN graph coordinates and topology
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
x = list(x)
y = list(y)
s = list(s)
t = list(t)


#save the coordinates and topology in pickle files
print('saving coordinates and topology files', flush=True)
with open(tmap_coordinates, 'wb+') as f:
    pickle.dump((x, y, dsmiles), f, protocol=pickle.HIGHEST_PROTOCOL)
with open(tmap_topology, 'wb+') as g:
    pickle.dump((s, t), g, protocol=pickle.HIGHEST_PROTOCOL)

#configuration of TMAP visualization
f = Faerun.Faerun(view='front', coords=False, clear_color='#ffffff', alpha_blending=True)
f.add_scatter('tmap_1M', {'x': x,
                            'y': y,
                            'c': [
                            list(fragments.hac.values),
                            list(fragments.mw.values),
                            list(fragments.alogp.values),
                            list(fragments.hba.values),
                            list(fragments.hbd.values),
                            list(fragments.psa.values),
                            list(fragments.rotb.values),
                            list(fragments.arom.values),
                            list(fragments.alerts.values)],
                            'labels': fragments['smiles']},
             shader='smoothCircle',
             point_scale=2.0,
             max_point_size=20,
             interactive=True,
             categorical=[False, False, False, False, False, False, False, False, False],
             colormap=['viridis', 'viridis', 'viridis', 'viridis', 'viridis', 'viridis', 'viridis', 'viridis', 'viridis'],
             series_title=[
                 'Heavy atom count',
                 'Molecular weight',
                 'ALogP',
                 'Hydrogen bond acceptors',
                 'Hydrogen bond donors',
                 'Polar surface area',
                 'Rotatable bonds',
                 'Aromatic rings',
                 'Alerts'],
             has_legend=True,
             legend_title='REAL fragments')
f.add_tree('tmap_1M__tree', {'from': s, 'to': t}, point_helper='tmap_1M', color='#222222')

f.plot('tmap_1M', template='smiles')
print(f'Creation of TMAP took {(timer() - start) / 3600} hours.')

#save TMAP data
with open(tmap_plot, 'wb+') as handle:
    pickle.dump(f.create_python_data(), handle, protocol=pickle.HIGHEST_PROTOCOL)

#start a Faerun web server
#host(tmap_plot, title='TMAP of REAL fragments', label_type='default', theme='dark')



