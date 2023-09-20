# TMAP

**Objective:**
This repository aims to visualize and cluster a highly diverse dataset obtained from the Human Metabolite Data Base. The input dataset comprises approximately 4,000 molecules. To use the script with a different input dataset than the one provided, please adhere to the following format in a CSV file:

```
idnumber,smiles
your_molecule,molecule_smiles
```

For example:

```
idnumber,smiles
HMDB0000058,Nc1ncnc2c1ncn2C1OC2COP(=O)(O)OC2C1O
```

**Execution:**
The visualization and clustering are achieved by calling a single script, `create_tmap.py`, which, in turn, utilizes other scripts to produce the desired outcome. In essence, this script follows several steps:

### 1. Data Preparation

By invoking `data_preparation.py`, the script prepares the data. This script takes SMILES and id_name as input and calculates various parameters:

- Desalted Smiles
- Heavy atom count (hac)
- Molecular weight (mw)
- Estimated lipophilicity (alogp)
- Number of hydrogen bond acceptors (hba)
- Number of hydrogen bond donors (hbd)
- Polar surface area (psa)
- Number of rotatable bonds (rotb)
- Number of aromatic rings (arom)
- Number of structural alerts (alerts)

Additionally, it generates other files in pickle (pkl) format that will be used in subsequent steps.

### 2. Create TMAP

TMAP, which stands for Topological Mapping, is a data visualization method designed to address the challenge of representing large datasets with millions of data points and high dimensionality in a two-dimensional format. It aims to visualize data while preserving both global and local features, providing sufficient detail for human inspection and interpretation.

![Reference](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-0416-x)

In essence, TMAP accomplishes these four steps:

1. LSH Forest indexing
2. Build c-approximate k-NN
3. Calculate MST of the k-NN
4. Lay the tree on the Euclidean plane using the OGDF modular C++ library

### 3. Find the Centroid

TMAP generates a tree structure, making it impossible to directly identify cluster centers. To overcome this limitation, the tree's topology is leveraged to create a DataFrame representing the hierarchy of the tree. This enables the selection of fragments that represent specific branches or portions.

To find the centroid of the tree, a node is identified where each branch contains fewer nodes than half the total number of nodes/fragments (n) in the tree. This process involves using a Python script (centroid.py), similar to the one used in OpenGenus IQ17 for centroid decomposition. The ID of the centroid node is saved to a text file (centroid.txt).

### 4. Hierarchy

The script identifies the children of the centroid node and adds them to a new column in the DataFrame. It also repeats this process for children of children while avoiding the addition of nodes that are already present in the DataFrame. The final DataFrame is converted to CSV for visualization purposes.
