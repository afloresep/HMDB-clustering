import numpy as np 
import pandas as pd 
from rdkit import Chem
from rdkit.Chem import QED
import os
import re



#get input and output files
directory = os.getcwd()
frag_prep_tmap = os.path.join(directory, 'HMDB_prep_tmap.pkl')
frag_hierarchy_df = os.path.join(directory, 'HMDB-hierarchy.pkl')

data = pd.read_csv('HMDB-smiles.csv')

def desalt(data):
    salts = ['\\.I', 'I\\.', '\\.Cl', 'Cl\\.', 'Br\\.']
    data = re.sub('|'.join(salts), '', data)
    return data


def calculate_parameters(df):
    # Create empty lists for each parameter
    smiles_desalted = []
    hac_list = []
    mw_list = []
    alogp_list = []
    hba_list = []
    hbd_list = []
    psa_list = []
    rotb_list = []
    arom_list = []
    alerts_list = []

    # Loop through each SMILES code in the DataFrame
    for smiles in df['smiles']:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Handle invalid SMILES here if needed
            continue

        p = QED.properties(mol)

        # Append parameter values to respective lists
        smiles_desalted.append(desalt(smiles))
        hac_list.append(Chem.rdchem.Mol.GetNumHeavyAtoms(mol))
        mw_list.append(p.MW)
        alogp_list.append(p.ALOGP)
        hba_list.append(p.HBA)
        hbd_list.append(p.HBD)
        psa_list.append(p.PSA)
        rotb_list.append(p.ROTB)
        arom_list.append(p.AROM)
        alerts_list.append(p.ALERTS)

    # Add the new parameter columns to the DataFrame
    df['smiles_desalted'] = smiles_desalted
    df['hac'] = hac_list
    df['mw'] = mw_list
    df['alogp'] = alogp_list
    df['hba'] = hba_list
    df['hbd'] = hbd_list
    df['psa'] = psa_list
    df['rotb'] = rotb_list
    df['arom'] = arom_list
    df['alerts'] = alerts_list

    return df[['smiles', 'smiles_desalted', 'hac', 'mw', 'alogp', 'hba', 'hbd', 'psa',
       'rotb', 'arom', 'alerts']]

fragments = calculate_parameters(data)
fragments = fragments.drop_duplicates(subset=['smiles'])
print(fragments.shape, flush=True)


fragments.to_pickle(frag_prep_tmap)
fragments.to_csv('HMDB_prep_tmap.csv', index=False) # Also save to csv