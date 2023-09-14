## Table of Contents
- [Overview](#overview)
- [Dataset](#dataset)
- [Code](#code)
- [Results](#results)
- [Standardization](#standardization)
- [Contributing](#contributing)
- [License](#license)

# Overview

# Dataset

Data used can be found in `/datasets/MR1.smi`. Sourced from HMDB.
Map4 fingerprints for these molecules are in `/datasets/MR1-map4-clean.xlsx` or `/datasets/MR1-map4-clean.csv`

# Code

### Original Approach (`original-aproach.ipynb`)

This method unfortunately does not yield the intended results. The reason being that we initially create molecular fingerprints based on similarity between all possible pairs of molecules. We then perform clustering on these pairs, effectively clustering molecules based on the similarity of their pair fingerprints. This approach doesn't capture the molecule's inherent similarity.

### Cluster Fingerprints (`cluster-fp.ipynb`)

Here, we modified the approach by moving away from the MFP_matrix and opting for K-means clustering of the fingerprints. However, this change didn't align with our expectations. The challenge lies in the fingerprint type; the morgan fingerprint we used wasn't optimal for HMDB's diverse dataset containing phospholipids, drug-like molecules, etc. Consequently, despite the code's correctness, the fingerprint failed to provide the needed information for effective clustering.

Addressing the fingerprint, we implemented 'map4' fingerprints from [reymon-group](https://github.com/reymond-group/map4/blob/master/README.md). We tested map4 with data from the Human Metabolome Database. Each of the ~4000 molecules' map4 fingerprints were exported and cleaned into a `MR1_map4-clean.csv` CSV file. Additionally, a smaller dataset of ~100 molecules was created for faster computational times.

We debated the pros and cons of standardizing our data before running PCA and clustering algorithms. This discussion is available in `test.ipynb` or at the end of this file.

Initially, we conducted a test using map4 fingerprints to observe their impact on clustering, detailed in `test.ipynb`. These results displayed higher sensitivity than the previous approach. While not ideal, this outcome is reasonable due to the dataset's chemical diversity. The outcome led us to believe that map4 could be valuable for our clustering problem, prompting us to take a more comprehensive approach in our code.

### Find Optimal K (`test-findK.ipynb`)

This notebook was created to determine the optimal K value for our dataset.


### Recursive K-means for std data (`std-recursive-kmeans.ipynb`)

This code is for testing what would happen if we did a new K-means clustering on highly populated clusters after doing a first one. 
If this goes correctly, I will try to implement a recursive method for clustering clusters with more than 20% of total samples.

Also, plotting was improved

## Standardize or Not to Standardize, That Is the Question

Typically, data is standardized before performing a PCA. 
However, in the case of map4 fingerprints, I find that standarization might not be necessary due to being messured in the same scale (despite variations spanning three orders of magnitud). 
One would use standarization when contrasting disparate units like the height of a tree and the girth of the same tree. 

Altough HMDB database contains diverse lipids, phospholipids... which could explain the difference in magnitud for some fingerprints? So maybe standarization is due?

When it comes to using PCA, data standarization is usually the go-to step. However, the situation is a bit different with map4 fingerprints. Despite having variations that span three orders of magnitude (1e05 / 1e08) these fingerprints are actually measured on a consistent scale (map4). This makes me wonder if standarization is really necessary...
Statistics is not my area of expertise so i could be wrong!
We're comparing apples to apples so to speak... For instance, you'd standardize when comparing something like the height of a tree and the tree's girth – they're completely different units.


# Results 

So, for what we can see in the `test.ipynb` notebook, standarization does not affect how the molecules are clustered... 




# TMAP 

### `data_preparation.py` 
prepares data. Input should be a csv with the following format: 

```
idnumber,smiles
HMDB0000058,Nc1ncnc2c1ncn2C1OC2COP(=O)(O)OC2C1O
```

Script to preparate data. It takes SMILES and id_name and calculates different parameters.
    *	Heavy atom count (hac)
    *	Molecular weight (mw)
    *	Estimated lipophilicity (alogp)
    *	Number of hydrogen bond acceptors (hba)
    *	Number of hydrogen bond donors (hbd)
    *	Polar surface area (psa)
    *	Number of rotatable bonds (rotb)
    *	Number of aromatic rings (arom)
    *	Number of structural alerts (alerts)

    
### `HMDB-tmap.py`

This script is designed to create a TMAP visualization based on input data in a .pkl file format created by `data_preparation.py`. 
The TMAP visualization involves various configurations and options that can be customized using command-line arguments.

## Usage

### Running the Script

To run the script, open your terminal and execute the following command:

```
python HMDB-tmap.py input_file [options]
```

- `input_file`: This is the input data file in .pkl format that you want to use for creating the TMAP visualization.

### Command-Line Options

The script supports several command-line options that allow you to customize the TMAP creation process:

- `--LSH_dim`: LSH dimension (default: 1024)
- `--prefix_trees`: Prefix trees (default: 128)
- `--cfg_k`: k for layout configuration (default: 50)
- `--cfg_kc`: kc for layout configuration (default: 50)
- `--cfg_mmm`: mmm_repeats for layout configuration (default: 2)
- `--cfg_node`: node_size for layout configuration (default: 1/70)
- `--cfg_sl_steps`: sl_extra_scaling_steps for layout configuration (default: 10)
- `--cfg_sl_repeats`: sl_repeats for layout configuration (default: 2)
- `--map4_dim`: MAP4 dimension (default: 1024)

Example:

```
python HMDB-tmap.py HMDB-smiles --map4_dim 512
```

This command runs the script using `HMDB-smiles.pkl` as input data and modifies the `map4_dim` option to 512, overriding the default value.

### Output

The script will create a TMAP visualization and save it as `tmap.html` and `tmap.js` files in the same directory as the script. It will also move these files and any files starting with the input file name into a new folder with the same name as the input file.

### Images

![TMAP generated from ~4000 molecules from HMDB](https://github.com/afloresep/HMDB-clustering/blob/master/tmap/Screenshot%202023-09-12%20at%2017.59.40.png)




