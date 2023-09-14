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

import argparse
import shutil

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

    move(args.input_file)


if __name__ == "__main__":
    main()

