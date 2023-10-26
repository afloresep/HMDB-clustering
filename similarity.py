import os
import pickle
import pandas as pd
import numpy as np
import networkx as nx
from map4 import MAP4Calculator
from rdkit import Chem
from timeit import default_timer as timer
import argparse
import shutil
from pyvis.network import Network
from rdkit import Chem
from tmap import tm
from map4 import MAP4Calculator
from multiprocessing import Pool
import time
from mhfp.encoder import MHFPEncoder
import faerun as Faerun


def calculate_similarity(args):

    mol_reference, dataset, dimensions, similarity_threshold, id = args
    MAP4 = MAP4Calculator(dimensions=dimensions)
    enc = tm.Minhash(dimensions)  # Create a new Minhash object in each worker
    similarity_dict = {}
    map4_reference = MAP4.calculate(mol_reference)
    
    for id2, smiles2 in zip(dataset['idnumber'], dataset['smiles']):
        mol_b = Chem.MolFromSmiles(smiles2)
        map4_b = MAP4.calculate(mol_b)
        distance = enc.get_distance(map4_reference, map4_b)
        if distance < similarity_threshold:
            similarity_dict[id2] = (smiles2, distance)
    
    return id, similarity_dict

def main():
    dim = 1024
    dataset = pd.read_csv('HMDB-smiles_prep_tmap.csv')
    representative_path = 'output/representative_nodes.txt'
    
    with open(representative_path, 'r') as file:
        similarity_threshold = 0.6
        pool = Pool()  # Create a pool of worker processes
        similarity_dicts = []

        start_time = time.time()


        for line in file:
            smiles, id = line.strip().split(",")
            mol_reference = Chem.MolFromSmiles(smiles)
            print(f'Searching similar molecules for {id}')
            args = (mol_reference, dataset, dim, similarity_threshold, id)
            similarity_dicts.append(pool.apply(calculate_similarity, (args,)))

        pool.close()
        pool.join()
        end_time = time.time()  # Record the end time

        elapsed_time = end_time - start_time
        print(f'Total time: {elapsed_time:.2f} seconds')

    similarity_dict = {id: sim_dict for id, sim_dict in similarity_dicts}

    df = pd.DataFrame(similarity_dict.items(), columns=['id', 'similarity'])
    df['smiles'], df['distance'] = zip(*df['similarity'])

    excel_file = 'output_optimized.xlsx'
    df.to_excel(excel_file, index=False)

if __name__ == "__main__":
    main()
