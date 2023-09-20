import csv


"""
Script to change from csv to cxsmiles in desired format 
"""

# Input and output file names
input_file = 'HMDB-smiles.csv'
output_file = 'HMDB-smiles.cxsmiles'

# Open the input and output files
with open(input_file, 'r') as csv_file, open(output_file, 'w', newline='') as output_csv:
    # Create a CSV reader and writer
    reader = csv.DictReader(csv_file)
    fieldnames = ['smiles', 'idnumber', 'Type']
    writer = csv.DictWriter(output_csv, fieldnames=fieldnames, delimiter='\t')
    
    # Write the header
    writer.writeheader()
    
    # Iterate through the rows in the input CSV and write to the output CSV
    for row in reader:
        writer.writerow(row)

print(f"File '{output_file}' has been created")
