"""
Reducen cxsmiles to a more feasable file
Doing this locally so cannot compile 17M fragments
"""

import pandas as pd 
from csv import reader
# Specify the path to your CXSMILES file
cxsmiles_file = "fragments_17M.cxsmiles"
shorter_cxsmiles_file = "shorter.cxsmiles"


# Number of lines to extract for testing
num_lines_to_extract = 100  # Adjust as needed

# Initialize a list to store the selected lines
test_lines = []

# Open the CXSMILES file and read lines
with open(cxsmiles_file, 'r') as file:
    for i, line in enumerate(file):
        # Append the current line to the test_lines list
        test_lines.append(line)

        # Stop reading after the desired number of lines
        if i + 1 >= num_lines_to_extract:
            break
print(test_lines)


# Save the selected lines as a new CXSMILES file
with open(shorter_cxsmiles_file, 'w') as new_file:
    new_file.writelines(test_lines)

print(f"Saved {num_lines_to_extract} lines to {shorter_cxsmiles_file}")

#load fragments dataset
with open(shorter_cxsmiles_file) as h:
    h_read = reader(h, delimiter='\t')
    frag = pd.DataFrame(h_read)


print(frag)
frag = frag.rename(columns=frag.iloc[0]).iloc[1:]
frag = frag[['smiles', 'idnumber']]
print(frag.shape, flush=True)


fragments = frag.copy()

print(fragments)