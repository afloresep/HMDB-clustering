"""
tmap for the  dataset.
Tmap 
1. LSH Forest indexing 
2. Build c-approximate k-NN
3. Calculate MST of the k-NN
4. Lay tree on Euclidian plane using OGDF modular C++ library
"""

import subprocess

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

import argparse
import shutil
from pyvis.network import Network





def create_tmap(pkl_file_name, 
                LSH_dim = 1024, 
                prefix_trees = 128, 
                cfg_k=50, 
                cfg_kc=50,
                cfg_mmm=2, 
                cfg_node=1/70,
                cfg_sl_steps=10, 
                cfg_sl_repeats=2, 
                map4_dim = 1024,
                ):
    
    lf = tm.LSHForest(LSH_dim, prefix_trees)

    directory = os.getcwd()
    fragments = pd.read_pickle(pkl_file_name+'_prep_tmap.pkl')
    tmap_coordinates = os.path.join(directory, pkl_file_name+'_coordinates.dat')
    tmap_plot = os.path.join(directory, pkl_file_name+'_tmap')    
    tmap_topology = os.path.join(directory, pkl_file_name+'_topology.dat')

    print('\n')
    print('Data loaded', flush=True)
    print(f'Data shape (molecules, features) = {fragments.shape}', flush=True)
    print('\n')

    MAP4 = MAP4Calculator(dimensions = map4_dim)
    ENC  = tm.Minhash(map4_dim)

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
    print('\n')
    print('starting LSH Forest indexing', flush=True)

    lf.batch_add(fps)
    lf.index()

    cfg = tm.LayoutConfiguration()
    cfg.k = cfg_k # Number of nearest neighbors to be calculated
    cfg.kc = cfg_kc # Scalar by which is multiplied (kc) 
    cfg.mmm_repeats = cfg_mmm # Number of repeats of the per-level layout algorithm 
    cfg.node_size = cfg_node # Size of the nodes
    cfg.sl_extra_scaling_steps = cfg_sl_steps # Number of repeats of scaling
    cfg.sl_repeats = cfg_sl_repeats # Number of repeats of scaling layout algorithm 

    #creation of k-NN graph coordinates and topology
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
    x = list(x)
    y = list(y)
    s = list(s)
    t = list(t)

    tmap_name = pkl_file_name+'_tmap'

    #save the coordinates and topology in pickle files
    print('saving coordinates and topology files', flush=True)
    with open(tmap_coordinates, 'wb+') as f:
        pickle.dump((x, y, dsmiles), f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(tmap_topology, 'wb+') as g:
        pickle.dump((s, t), g, protocol=pickle.HIGHEST_PROTOCOL)


    #configuration of TMAP visualization
    f = Faerun.Faerun(view='front', coords=False, clear_color='#ffffff', alpha_blending=True)
    f.add_scatter('tmap_', {'x': x,
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
    f.add_tree('tmap__tree', {'from': s, 'to': t}, point_helper='tmap_', color='#222222')

    f.plot('tmap_', template='smiles')
    print(f'Creation of TMAP took {(timer() - start) / 3600} hours.')

    #save TMAP data
    with open(tmap_plot, 'wb+') as handle:
        pickle.dump(f.create_python_data(), handle, protocol=pickle.HIGHEST_PROTOCOL)

    #start a Faerun web server
    #host(tmap_plot, title='TMAP of REAL fragments', label_type='default', theme='dark')


def hierarchy(pkl_file_name):
    

    directory = os.getcwd()
    fragments = pd.read_pickle(pkl_file_name+'_prep_tmap.pkl')

    original_csv = pd.read_csv(pkl_file_name+'.csv')
    tmap_plot = os.path.join(directory, pkl_file_name+'_tmap')    
    centroid_txt = os.path.join(directory, 'centroid.txt')
    frag_neighbors_df = os.path.join(directory, 'neighbors.pkl')
    frag_hierarchy_df = os.path.join(directory, 'hierarchy.pkl')
    frags_plevel = os.path.join(directory, 'fragments_per_level.csv')
    frag_hierarchy_plot = os.path.join(directory, 'hierarchy_plot.html')


    tmap_frag_topology = os.path.join(directory, pkl_file_name+'_topology.dat')
    

    tmap_frag_coordinates = os.path.join(directory, pkl_file_name+'_coordinates.dat')
    
    frag_prep_tmap = os.path.join(directory, pkl_file_name+'_prep_tmap.pkl')

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
    #print(frag_plevel)

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


    hierarchy = pd.read_csv('Hierarchy.csv')


    index_to_idnumber = {int(row[0]+2): row[1] for row in original_csv[['idnumber']].itertuples()}

    def replace_with_idnumber(index):
        return index_to_idnumber.get(index, index)

    # Apply the mapping to the entire first DataFrame
    first = hierarchy.applymap(replace_with_idnumber)
    first.to_csv('Hierarchy.csv')


def FindCentroid(file_name):
    #get input and output files
    directory = os.getcwd()
    tmap_frag_topology = os.path.join(directory, file_name+'_topology.dat')
    frag_prep_tmap = os.path.join(directory, file_name+'_prep_tmap.pkl')
    tmap_frag_coordinates = os.path.join(directory, file_name+'_coordinates.dat')
    centroid_file = os.path.join(directory, 'centroid.txt')

    #get fragment and topology data
    fragments = pd.read_pickle(frag_prep_tmap)
    coordinates = list(pd.read_pickle(tmap_frag_coordinates))
    coordinates = pd.DataFrame(coordinates, index=['x', 'y', 'dsmiles']).transpose()
    print(coordinates.shape, flush=True)
    topology = list(pd.read_pickle(tmap_frag_topology))
    topology = pd.DataFrame(topology, index=['src', 'dest']).transpose()
    print(topology.shape, flush=True)

    ##find the centroid of the TMAP tree
    class Tree:
        def __init__(self, n):
            self.size = n + 1
            self.cur_size = 0
            self.tree = [[] for _ in range(self.size)]
            self.iscentroid = [False] * self.size
            self.ctree = [[] for _ in range(self.size)]

        def dfs(self, src, visited, subtree):
            visited[src] = True
            subtree[src] = 1
            self.cur_size += 1
            for adj in self.tree[src]:
                if not visited[adj] and not self.iscentroid[adj]:
                    self.dfs(adj, visited, subtree)
                    subtree[src] += subtree[adj]

        def findCentroid(self, src, visited, subtree):
            iscentroid = True
            visited[src] = True
            heavy_node = 0
            for adj in self.tree[src]:
                if not visited[adj] and not self.iscentroid[adj]:
                    if subtree[adj] > self.cur_size//2:
                        iscentroid = False
                    if heavy_node == 0 or subtree[adj] > subtree[heavy_node]:
                        heavy_node = adj
            if iscentroid and self.cur_size - subtree[src] <= self.cur_size//2:
                return src
            else:
                return self.findCentroid(heavy_node, visited, subtree)

        def findCentroidUtil(self, src):
            visited = [False] * self.size
            subtree = [0] * self.size
            self.cur_size = 0
            self.dfs(src, visited, subtree)
            for i in range(self.size):
                visited[i] = False
            centroid = self.findCentroid(src, visited, subtree)
            self.iscentroid[centroid] = True
            return centroid

        def decomposeTree(self, root):
            centroid = self.findCentroidUtil(root)
            print('centroid: ', flush=True)
            print(centroid)
            return centroid

        def addEdge(self, src, dest):
            self.tree[src].append(dest)
            self.tree[dest].append(src)

    start = timer()


    tree = Tree(len(coordinates))
    for i in topology.itertuples():
        tree.addEdge(i.src, i.dest)
    centroid = tree.decomposeTree(1)

    print(f'Finding centroid took: {timer() - start}sec.')

    file = open(centroid_file, 'w')
    file.write(str(centroid))
    file.close()


def move(output_folder):

    # Move files starting with args.input_file into the new folder
    new_directory = 'TMAP-'+output_folder
    os.makedirs(new_directory, exist_ok=True)

    for filename in os.listdir():
        if filename.startswith(output_folder):
            shutil.move(filename, os.path.join(new_directory, filename))

    # Move tmap.html and tmap.js files into the new folder
    shutil.move('tmap_.html', os.path.join(new_directory, 'tmap_'+output_folder+'.html'))
    shutil.move('tmap_.js', os.path.join(new_directory, 'tmap_.js'))
    shutil.move('hierarchy.csv', os.path.join(new_directory, 'hierarchy.csv'))
    shutil.move('fragments_per_level.csv', os.path.join(new_directory, 'fragments_per_level.csv'))
    shutil.move('centroid.txt', os.path.join(new_directory, 'centroid.txt'))
    shutil.move('hierarchy.pkl', os.path.join(new_directory, 'hierarchy.pkl'))
    shutil.move('neighbors.pkl', os.path.join(new_directory, 'neighbors.pkl'))

    
def main():
    parser = argparse.ArgumentParser(description="Create tmap.")
    parser.add_argument("input_file", help="Input .pkl file name generated from data_prep")
    parser.add_argument("--LSH_dim", type=int, default=1024, help="LSH dimension (default: 1024)")
    parser.add_argument("--prefix_trees", type=int, default=128, help="Prefix trees (default: 128)")
    parser.add_argument("--cfg_k", type=int, default=50, help="k for layout configuration (default: 50)")
    parser.add_argument("--cfg_kc", type=int, default=50, help="kc for layout configuration (default: 50)")
    parser.add_argument("--cfg_mmm", type=int, default=2, help="mmm_repeats for layout configuration (default: 2)")
    parser.add_argument("--cfg_node", type=float, default=1/70, help="node_size for layout configuration (default: 1/70)")
    parser.add_argument("--cfg_sl_steps", type=int, default=10, help="sl_extra_scaling_steps for layout configuration (default: 10)")
    parser.add_argument("--cfg_sl_repeats", type=int, default=2, help="sl_repeats for layout configuration (default: 2)")
    parser.add_argument("--map4_dim", type=int, default=1024, help="MAP4 dimension (default: 1024)")
    
    args = parser.parse_args()

    
    subprocess.run(["python", "data_prep.py", args.input_file])

    create_tmap(pkl_file_name=args.input_file, 
        LSH_dim=args.LSH_dim,
        prefix_trees=args.prefix_trees,
        cfg_k=args.cfg_k,
        cfg_kc=args.cfg_kc,
        cfg_mmm=args.cfg_mmm,
        cfg_node=args.cfg_node,
        cfg_sl_steps=args.cfg_sl_steps,
        cfg_sl_repeats=args.cfg_sl_repeats,
        map4_dim=args.map4_dim
    )


    FindCentroid(args.input_file)

    #    subprocess.run(["python", 'centroid.py'])
    

    hierarchy(pkl_file_name=args.input_file)

    move(args.input_file)


if __name__ == "__main__":
    main()

