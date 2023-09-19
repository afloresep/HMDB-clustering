## Table of Contents
- [Overview](#overview)
- [Dataset](#dataset)
- [Code](#code)
- [Results](#results)
- [Standardization](#standardization)
- [Contributing](#contributing)
- [License](#license)

# Overview
![TMAP generated from ~4000 molecules from HMDB](https://github.com/afloresep/HMDB-clustering/blob/master/tmap/Screenshot%202023-09-12%20at%2017.59.40.png)
The goal of this repository is to visually explore the structure of extensive and varied datasets through the utilization of TMAP and MAP4 fingerprint techniques.


# Dataset

Data used can be found in `/datasets/MR1.smi`. Sourced from HMDB.
Map4 fingerprints for these molecules are in `/datasets/MR1-map4-clean.xlsx` or `/datasets/MR1-map4-clean.csv`


# TMAP 

### `data_preparation.py` 
Prepares data. Input should be a csv with the following format: 

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

    
### `create_tmap.py`

This script is designed to create a TMAP visualization based on input data in a .pkl file format created by `data_preparation.py`. 
The TMAP visualization involves various configurations and options that can be customized using command-line arguments.

### `centroid.py`

Script designed to find centroid using topology files generated. 

## Usage


## Setting up the Conda Environment

To ensure that you have the required dependencies to run this project, you can create a Conda environment based on the provided `tmap_environment.yml` file. Follow these steps:

1. **Clone the Repository**: First, clone this repository to your local machine:

   ```
   git clone https://github.com/yourusername/yourproject.git
   ```

2. **Navigate to the Project Directory**: Change your current directory to the project's root directory:

   ```
   cd yourproject
   ```

3. **Create a New Conda Environment**: Create a new Conda environment using the provided YAML file. You can choose to use the same environment name (`tmap`) or specify a different one:

   ```
   conda env create --name tmap --file tmap_environment.yml
   ```

4. **Activate the Environment**: Activate the newly created environment:

   ```
   conda activate tmap
   ```

Now, you have set up the Conda environment with all the required packages. You can run the project within this environment.

To deactivate the environment when you're done, simply run:

```
conda deactivate
```

### Running the Script

To run the script, open your terminal and execute the following command:

```
python create_tmap.py input_file [options]
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
python create_tmap.py HMDB-smiles --map4_dim 512
```

This command runs the script using `HMDB-smiles.pkl` as input data and modifies the `map4_dim` option to 512, overriding the default value.

### Output

The script will create a TMAP visualization and save it as `tmap.html` and `tmap.js` files in the same directory as the script. It will also move these files and any files starting with the input file name into a new folder with the same name as the input file.


### Images

![TMAP generated from ~4000 molecules from HMDB](https://github.com/afloresep/HMDB-clustering/blob/master/tmap/Screenshot%202023-09-12%20at%2017.59.40.png)


