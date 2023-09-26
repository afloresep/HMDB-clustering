
import subprocess

#import packages
import os
import tmap as tm

import faerun as Faerun
import pandas as pd
import pickle
import networkx as nx

from faerun import host
from map4 import MAP4Calculator
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from timeit import default_timer as timer

import argparse
import shutil
from pyvis.network import Network
import pandas as pd


import networkx as nx
import matplotlib.pyplot as plt


def create_tmap(pkl_file_name, heads, 
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

    fragments = pd.read_csv(pkl_file_name+'_prep_tmap.csv')
    lf = tm.LSHForest(LSH_dim, prefix_trees)
    
    print('\n')
    print('Data loaded', flush=True)
    print(f'Data shape (molecules, features) = {fragments.shape}', flush=True)
    print('\n')

    MAP4 = MAP4Calculator(dimensions = map4_dim)

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

    ld = 0   
    for chunk in chunks2:
        ld += 1
        print(f'calculating MAP4 for chunk {ld}')
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

    x_values = list(x)
    y_values = list(y)
    t_id = []
    s_id = []
    ids = list(fragments['idnumber'])

    # Replace index for actual name of molecule
    for index in t:
        t_id.append(fragments['idnumber'][index])
    t_id.append(fragments['idnumber'][1])
    for index_2 in s:
        s_id.append(fragments['idnumber'][index_2])
    s_id.append(fragments['idnumber'][1])

    data = {
    'id': ids,
    'x': x_values,
    'y': y_values,
    's': s_id,  # Add a random value to make it 497 elements
    't': t_id,  # Add a random value to make it 497 elements
    # Create a DataFrame with name, x, y, s, t
    }

    df = pd.DataFrame(data) 

    # Save the DataFrame to a CSV file
    df.to_csv('graph_info.csv', index=False)

    print('Building graph...')
    # Create an empty graph
    G = nx.Graph()

    # Iterate through the DataFrame to add nodes
    for _, row in df.iterrows():
        node_id = row['id']
        x_coord = row['x']
        y_coord = row['y']
        G.add_node(node_id, x=x_coord, y=y_coord)

    # Iterate through the DataFrame to add edges
    for _, row in df.iterrows():
        source_node = row['id']
        target_nodes = [row['s'], row['t']]

        # Add edges between the source node and target nodes
        for target_node in target_nodes:
            G.add_edge(source_node, target_node)


    degree_centrality = nx.degree_centrality(G)

    # Sort nodes by degree centrality in descending order
    sorted_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)

    # Select the top n nodes as representatives
    n = heads  
    representative_nodes = [node for node, centrality in sorted_nodes[:n]]

    # save nodes along with their smiles 
    with open('representative_nodes.txt', 'w') as file:
        for node in representative_nodes:
            filtered_fragment = fragments[fragments['idnumber'] == node]
            if not filtered_fragment.empty:
                smiles_value = filtered_fragment.iloc[0]['smiles']
                file.write(f"{smiles_value},{node}\n")
    
    print('Saving nodes to .txt file')


    # Create a scatter plot
    plt.figure(figsize=(64, 48)) 

    # Loop through the rows of the DataFrame and plot each point
    for index, row in df.iterrows():
        x_coordinates = row['x']
        y_coordinates = row['y']
        point_id = row['id']
        
        # Determine the size of the point based on whether it's in the highlight list
        size = 50  # Default size for points not in the highlight list
        if point_id in representative_nodes:
            size = 150  # Larger size for points in the highlight list

        # Check if the point should be highlighted in red
        if point_id in representative_nodes:
            plt.scatter(x_coordinates, y_coordinates, label=point_id, c='red', s=size)
        else:
            plt.scatter(x_coordinates, y_coordinates, c='blue')

    # Add labels and legend
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()

    # Show the plot
    plt.title('Scatter Plot with Highlighted Points')
    plt.grid(True)
    plt.savefig('tmap_with_nodes.png', bbox_inches='tight')

    print('plot with nodes saved to png')


    def tmap_visualization(x,y,s,t):

        #configuration of TMAP visualization
        f = Faerun.Faerun(view='front', coords=False, clear_color='#ffffff', alpha_blending=True)
        f.add_scatter('tmap', {'x': x,
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
        
        f.add_tree('tmap_tree', {'from': s, 'to': t}, point_helper='tmap', color='#222222')

        f.plot('tmap', template='smiles')

    tmap_visualization(x,y,s,t)
    print('TMAP built and save to html')

def main():
    parser = argparse.ArgumentParser(description="Create tmap.")
    parser.add_argument("input_file", help="Input .pkl file name generated from data_prep")
    parser.add_argument("--heads", type=int, default=100, help="Number of molecules to be selected from the tree based on their centrality")
    parser.add_argument("--LSH_dim", type=int, default=1024, help="LSH dimension (default: 1024)")
    parser.add_argument("--prefix_trees", type=int, default=128, help="Prefix trees (default: 128)")
    parser.add_argument("--cfg_k", type=int, default=50, help="k for layout configuration (default: 50)")
    parser.add_argument("--cfg_kc", type=int, default=50, help="kc for layout configuration (default: 50)")
    parser.add_argument("--cfg_mmm", type=int, default=2, help="mmm_repeats for layout configuration (default: 2)")
    parser.add_argument("--cfg_node", type=float, default=1/70, help="node_size for layout configuration (default: 1/30)")
    parser.add_argument("--cfg_sl_steps", type=int, default=10, help="sl_extra_scaling_steps for layout configuration (default: 10)")
    parser.add_argument("--cfg_sl_repeats", type=int, default=2, help="sl_repeats for layout configuration (default: 2)")
    parser.add_argument("--map4_dim", type=int, default=1024, help="MAP4 dimension (default: 1024)")
    
    args = parser.parse_args()

    subprocess.run(["python", "data-preparation.py", args.input_file])

    create_tmap(pkl_file_name=args.input_file,
        heads=args.heads,
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


    #FindCentroid(args.input_file)

    #subprocess.run(["python", 'centroid.py'])
    

    #hierarchy(pkl_file_name=args.input_file)

    #move(args.input_file)


if __name__ == "__main__":
    main()
