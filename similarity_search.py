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

data = dataset['smiles'] 
idnumber = dataset['idnumber']
dict = {}

smiles_a = 'NC(=O)Cc1c[nH]c2ccccc12'
mol_a = Chem.MolFromSmiles(smiles_a)
map4_a = MAP4.calculate(mol_a)

for smiles, id in zip(data,idnumber):
    mol = Chem.MolFromSmiles(smiles)
    map4_b = MAP4.calculate(mol)
    if ENC.get_distance(map4_a, map4_b) < 0.7:
        print(f'{id} distance = {ENC.get_distance(map4_a, map4_b)}')
