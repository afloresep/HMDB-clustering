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

import tmap as tm
from map4 import MAP4Calculator

from multiprocessing import Pool
import time

from mhfp.encoder import MHFPEncoder
import faerun as Faerun


def calculate_similarity(args):

    mol_reference, dataset, dimensions, similarity_threshold, id = args
    MAP4 = MAP4Calculator(dimensions=dimensions)
    enc = tm.Minhash(dimensions)  # Create a new Minhash object in each worker
    similar_molecules = []
    
    map4_reference = MAP4.calculate(mol_reference)
    
    for id2, smiles2 in zip(dataset['idnumber'], dataset['smiles']):
        mol_2 = Chem.MolFromSmiles(smiles2)
        map4_2 = MAP4.calculate(mol_2)
        distance = enc.get_distance(map4_reference, map4_2)
        if distance < similarity_threshold:
            similar_molecules.append((id2, distance))
    
    return similar_molecules

def main():
    dim = 1024
    
    dataset = pd.read_csv('HMDB-smiles_prep_tmap.csv')
    representative_path = 'output/representative_nodes.txt'
    
    with open(representative_path, 'r') as file:
        similarity_threshold = 0.5
        pool = Pool()  # Create a pool of worker processes
        similarity_dicts = []

        start_time = time.time()

        similarity_dictionary = {}

        for line in file:
            smiles, id = line.strip().split(",")
            mol_reference = Chem.MolFromSmiles(smiles)
            print(f'Searching similar molecules for {id}')
            args = (mol_reference, dataset, dim, similarity_threshold, id)
            similarity_dictionary[id] = (pool.apply(calculate_similarity, (args,)))
            similarity_dicts.append(pool.apply(calculate_similarity, (args,)))

        pool.close()
        pool.join()
        end_time = time.time()  # Record the end time

        elapsed_time = (end_time - start_time)/60
        print(f'Total time: {elapsed_time:.2f} minutes')

    df = pd.DataFrame.from_dict(similarity_dictionary, orient='index')
    print(similarity_dictionary)
    # Save the DataFrame to an Excel file
    excel_file = 'similarity_search.xlsx'
    df.to_excel(excel_file, index=False)

if __name__ == "__main__":
    main()
