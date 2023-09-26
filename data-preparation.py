import os
import re
import argparse

import csv
import numpy as np 
import pandas as pd 


from rdkit import Chem
from rdkit.Chem import QED

def to_cxsmiles(input_file):
    """
    Script to change from csv to cxsmiles in desired format
    """
    # Create a 'data' subdirectory if it doesn't exist
    if not os.path.exists("data"):
        os.makedirs("data")

    # Define the fieldnames for the output CSV
    fieldnames = ['smiles', 'idnumber', 'Type']

    # Open the input and output files
    with open(input_file, 'r') as csv_file, open(os.path.join("data", input_file+'.cxsmiles'), 'w', newline='') as output_csv:
        # Create a CSV reader and writer
        reader = csv.DictReader(csv_file)
        writer = csv.DictWriter(output_csv, fieldnames=fieldnames, delimiter='\t')

        # Write the header
        writer.writeheader()

        # Iterate through the rows in the input CSV and write to the output CSV
        for row in reader:
            writer.writerow(row)

    print(f"File 'data/{input_file}.cxsmiles' has been created")

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

    return df


def main():

    parser = argparse.ArgumentParser(description="Convert a CSV file to cxsmiles format.")
    parser.add_argument("input_file", help="Input CSV file name")
    args = parser.parse_args()

    data = pd.read_csv(args.input_file+'.csv')

    fragments = calculate_parameters(data)
    fragments = fragments.drop_duplicates(subset=['smiles'])

    print(fragments.shape, flush=True)
    fragments.to_csv(args.input_file+'_prep_tmap.csv', index=False) # Also save to csv


if __name__ == "__main__":
    main()
