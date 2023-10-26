import pandas as pd
import subprocess

#import packages
import os
import tmap as tm

import faerun as Faerun
import numpy as np
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


from rdkit import Chem
import tmap as tm
from map4 import MAP4Calculator


dim = 1024

MAP4 = MAP4Calculator(dimensions=dim)
ENC = tm.Minhash(dim)

dataset = pd.read_csv('HMDB-smiles_prep_tmap.csv')
representative_path = 'output/representative_nodes.txt'

data = dataset['smiles'] 
idnumber = dataset['idnumber']
similarity_dict = {}


with open(representative_path, 'r') as file:
    # Read each line in the file
    for line in file:
        smiles, id = line.strip().split(",")
        mol_reference = Chem.MolFromSmiles(smiles)
        map4_reference = MAP4.calculate(mol_reference)
        similar_molecules = []
        print(f'Searching similar molecules for {id}') 
        for smiles2, id2 in zip(data,idnumber):
            mol_b = Chem.MolFromSmiles(smiles2)
            map4_b = MAP4.calculate(mol_b)
            if ENC.get_distance(map4_reference, map4_b) < 0.7:
                similar_molecules.append((id2, smiles2, ENC.get_distance(map4_reference, map4_b)))
                print(f'{id2} was added with a distance of {ENC.get_distance(map4_reference, map4_b)}')
            similarity_dict[id] = similar_molecules
        print(f'for {id}: {len(similar_molecules)} molecules found')


# Create a DataFrame from the dictionary
df = pd.DataFrame.from_dict(similar_molecules, orient='index')

# Transpose the DataFrame for a more readable structure
df = df.transpose()

# Save the DataFrame to an Excel file
excel_file = 'output.xlsx'
df.to_excel(excel_file, index=False)
