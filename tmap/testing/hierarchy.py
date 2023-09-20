import os
import pandas as pd
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from pyvis.network import Network
from timeit import default_timer as timer

#get input and output files
directory = os.getcwd()
tmap_frag_topology = os.path.join(directory, 'TMAP_topology.dat')
tmap_frag_coordinates =  os.path.join(directory,'HMDB-smiles_coordinates.dat')
frag_prep_tmap = os.path.join(directory, 'HMDB-smiles_prep_tmap.pkl')
centroid_txt = os.path.join(directory, 'centroid.txt')
frag_neighbors_df = os.path.join(directory, 'neighbors.pkl')
frag_hierarchy_df = os.path.join(directory, 'hierarchy.pkl')
frags_plevel = os.path.join(directory, 'fragments_per_level.csv')
frag_hierarchy_plot = os.path.join(directory, 'hierarchy_plot.html')


directory = os.getcwd()
tmap_frag_topology = os.path.join(directory, 'HMDB-smiles_topology.dat')
tmap_frag_coordinates = os.path.join(directory, 'HMDB-smiles_coordinates.dat')
frag_prep_tmap = os.path.join(directory, 'HMDB-smiles_prep_tmap.pkl')



#get fragment and topology data
fragments = pd.read_pickle(frag_prep_tmap)

coordinates = list(pd.read_pickle(tmap_frag_coordinates))
coordinates = pd.DataFrame(coordinates, index=['x', 'y', 'dsmiles']).transpose()

topology = list(pd.read_pickle(tmap_frag_topology))
topology = pd.DataFrame(topology, index=['src', 'dest']).transpose()

with open(centroid_txt, 'r') as ctxt:
    centroid = int(ctxt.read())

print(f'shape fragments: {fragments.shape} \nshape coordinates: {coordinates.shape} \nshape topology: {topology.shape} \ncentroid: {centroid}', flush=True)
start = timer()

##get the hierarchy of the TMAP tree
#set network and solver
nt = Network(layout='Hierarchical')
nt.hrepulsion()

#get edge_data for visual
nt_data = topology.copy()
sources = nt_data['src']
destinations = nt_data['dest']
edge_data = zip(sources, destinations)

#add nodes and edges to network
for i in edge_data:
    src = i[0]
    dest = i[1]
    nt.add_node(src, src, title=str(src))
    nt.add_node(dest, dest, title=str(dest))
    nt.add_edge(src, dest)

#get neighbors per node dataframe and save file
neighbor_map = nt.get_adj_list()
neighbors = pd.DataFrame.from_dict(neighbor_map, orient='index').transpose()
neighbors.to_pickle(frag_neighbors_df)
print('neighbors per node df saved as neighbors.pkl', flush=True)
print(f'creating neighbors per node df took: {(timer() - start) / 3600} hours.', flush=True)

#def to create df containing parent and first children
def hierarchy_tmap(centroid, neighbors):
    hierarchy = pd.DataFrame([centroid]).astype(int)
    child = neighbors[centroid].dropna().astype(int)
    neighbors = neighbors.drop(columns=centroid)
    hierarchy = pd.concat([hierarchy, child], ignore_index=True, axis=1)
    return hierarchy, neighbors

#def to add children columns to hierarchy dataframe
def add_children(hierarchy, neighbors):
    child_ls = []
    children = []
    for index, i in hierarchy.iloc[:, -1].dropna().items():
        child = neighbors[i].dropna().astype(int)
        child_ls.extend(child)
        [children.append(x) for x in child_ls if x not in children]
        neighbors = neighbors.drop(columns=i)
    for i in pd.Series(children, dtype='int'):
        if i in hierarchy.values:
            children.remove(i)
    hierarchy = pd.concat([hierarchy, pd.Series(children, dtype='float64')], ignore_index=True, axis=1)
    return hierarchy, neighbors

start = timer()

#create hierarchy dataframe and the reduced version of the neighbors dataframe
print('create df containing centroid and first children', flush=True)
hierarchy, neighbors = hierarchy_tmap(centroid, neighbors)

print('add children to hierarchy dataframe', flush=True)

#add children to hierarchy dataframe until the neighbors dataframe is empty (i.e. stop when all nodes are in the hierarchy dataframe)
while neighbors.empty == False and hierarchy.iloc[:, -1].isna().sum() != len(hierarchy):
    hierarchy, neighbors = add_children(hierarchy, neighbors)

#remove last column of hierarchy dataframe if it only contains NaN values
print('remove last column if NaN-only', flush=True)
nan_last_col = hierarchy.iloc[:,-1].isna().sum()
if nan_last_col == len(hierarchy):
    hierarchy = hierarchy.iloc[:,:-1]

print(f'creating hierarchy df of TMAP fragments tree took: {(timer() - start) / 3600} hours.', flush=True)
hierarchy.to_pickle(frag_hierarchy_df)
frag_plevel = hierarchy.count()
frag_plevel.to_csv(frags_plevel)
print('hierarchy df and number of values per level saved as hierarchy.pkl and fragments_per_level.csv', flush=True)

print('number of fragments per level', flush=True)
pd.set_option('display.max_rows', None)

print('hierarchy of the generated fragments tree:       (parent: column 0)', flush=True)

hierarchy.to_csv('Hierarchy.csv')
start = timer()

#get level data from hierarchy data
value_cols = pd.DataFrame(hierarchy.stack().values.astype(int))
levels_src = pd.DataFrame()

for column in hierarchy:
    data = hierarchy[column].dropna().astype(int)
    val_lev = pd.DataFrame(data)
    val_lev = val_lev.rename(columns={column: 'src'})
    val_lev['level_src'] = column
    levels_src = pd.concat([levels_src, val_lev]).reset_index(drop=True)
levels_dest = levels_src.copy().rename(columns={'src': 'dest', 'level_src': 'level_dest'})

#merge src and dest level data with topology data
nt2_data = topology.copy()
nt2_data = nt2_data.merge(levels_src, how='left', on='src')
nt2_data = nt2_data.merge(levels_dest, how='left', on='dest')

#create hierarchical layout network based on level
nt2 = Network(height='100%', width='100%', font_color='blue', layout='Hierarchical', heading='REAL fragments')
nt2.hrepulsion()

#get edge_data2 for visual
sources = nt2_data['src']
destinations = nt2_data['dest']
levels_src = nt2_data['level_src']
levels_dest = nt2_data['level_dest']
edge_data2 = zip(sources, destinations, levels_src, levels_dest)

#add nodes with level and edges to network
for i in edge_data2:
    src = i[0]
    dest = i[1]
    level_src = i[2]
    level_dest = i[3]
    nt2.add_node(src, src, title=str(src), level=level_src)
    nt2.add_node(dest, dest, title=str(dest), level=level_dest)
    nt2.add_edge(src, dest)

neighbor_map2 = nt2.get_adj_list()

#set size of nodes based on number of neighbors
for node in nt2.nodes:
    node['title'] += ' Neighbors: <br>' + '<br>'.join(str(neighbor_map2[node['id']]))
    node['value'] = len(neighbor_map2[node['id']])



