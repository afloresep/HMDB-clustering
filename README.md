# TMAP

![TMAP generated from ~4000 molecules from HMDB](https://github.com/afloresep/HMDB-clustering/blob/master/tmap/Screenshot%202023-09-12%20at%2017.59.40.png)

This repository aims to visualize and a highly diverse dataset obtained from the Human Metabolite Data Base. Then, based on the graph, select as many representative as the user wants.
The input dataset comprises approximately 4,000 molecules. To use the script with a different input dataset than the one provided, please adhere to the following format in a CSV file:

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
The visualization and selection of molecules are achieved by calling a single script, `main.py`, which, in turn, utilizes other scripts to produce the desired outcome. In essence, this script follows several steps:

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

All these will be used when visualizing our molecules in TMAP. Depending on their value on each parameter each node color will variate. 

### 2. Create TMAP

TMAP, which stands for Topological Mapping, is a data visualization method designed to address the challenge of representing large datasets with millions of data points and high dimensionality in a two-dimensional format. It aims to visualize data while preserving both global and local features, providing sufficient detail for human inspection and interpretation.

[Reference](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-0416-x)

In essence, TMAP accomplishes these four step:

1. [LSH Forest indexing](#lsh-forest-indexing)
2. [Build c-approximate k-NN](#c-approximate-knn)
3. [Calculate MST of the k-NN](#calculate-mst)
4. [Lay the tree on the Euclidean plane using the OGDF modular C++ library](#ogdf-library)

You can read more about this at the bottom of the page. 

### 3. Select n representative molecules

TMAP generates a network structure, making it impossible to directly identify centers.To overcome this, the code employs the concept of degree centrality, a fundamental measure in graph theory. Degree centrality assesses the importance of a node within a network based on the number of edges (connections) it has. In the context of your large dataset visualization, it helps identify nodes that play pivotal roles or serve as key junctions in the data structure.

Degree centrality is a measure of how well-connected a node is within the graph. A higher degree centrality indicates that a particular node is connected to more neighbors, making it a potentially significant node in terms of information flow or connectivity
The number of nodes selected is defined by the variable --heads (default= 100)


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
   conda create env --name tmap --file tmap_environment.yml
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

To run the script, open your terminal, go to tmap folder and execute the following command:

```
python create_tmap.py HMDB-smiles 
```

### Command-Line Options

The script supports several command-line options that allow you to customize the TMAP creation process:

- `--LSH_dim`: LSH dimension (default: 1024)
- `--heads`: Number of representative nodes selected from the graph. (default:100)
- `--prefix_trees`: Prefix trees (default: 128): More trees can lead to better accuracy but may also increase computation time and memory usage.
- `--cfg_k`: k for layout configuration (default: 50).  Number of nearest neighbors or connections considered when arranging data points in the layout.
- `--cfg_kc`: kc for layout configuration (default: 50)
- `--cfg_mmm`: mmm_repeats for layout configuration (default: 2)
- `--cfg_node`: node_size for layout configuration (default: 1/70). Size of individual nodes or data points in the layout. Influence how spaced out or concentrated the data points appear in the visualization
- `--cfg_sl_steps`: sl_extra_scaling_steps for layout configuration (default: 10)
- `--cfg_sl_repeats`: sl_repeats for layout configuration (default: 2)
- `--map4_dim`: MAP4 dimension (default: 1024)

Example:

```
python main.py HMDB-smiles --heads 200
```

### Output

The script will create a TMAP visualization and save it as `.html` and `.js` files in the same directory as the script. 



### Images

![TMAP generated from ~4000 molecules from HMDB](https://github.com/afloresep/HMDB-clustering/blob/master/tmap/Screenshot%202023-09-12%20at%2017.59.40.png)




### LSH Forest indexing
<!-- Anchor point for LSH Forest indexing section -->
<a name="lsh-forest-indexing"></a>
Here's an improved version of the text with some clarifications and readability enhancements:

LSH, short for Locality-Sensitive Hashing, is a technique designed to approximate the similarity between pairs of high-dimensional data points. This method proves particularly valuable for efficiently identifying similar items within large datasets.

**How LSH Works: Hashing Data Points**

At the core of LSH is the process of "hashing" data points. Hashing involves applying a specialized function to a data point, transforming it into a fixed-size hash code. This hash function takes a data point as input and produces a unique output that represents that specific data point.

The hash funciton is designed in a way that ensures similar data points are mapped to the same or nearby hash codes with a high degree of probability (usually described as falling 'in the same cube'). In other words, if two data points share similarity, their resulting hash codes should be very close to each other within the hash code space. This strategic design enables efficient searches for similar items.

LSH Forest Indexing takes LSH a step further by introducing a forest of multiple hash tables instead of just one. When querying a data point, LSH Forest performs an approximate nearest neighbor search by simultaneously examining multiple hash tables. This approach significantly improves accuracy when compared to traditional LSH.


### Build c-approximate k-NN
<!-- Anchor point for Build c-approximate k-NN section -->
<a name="c-approximate-knn"></a>
Given the LSH Forest index, this step involves finding c-approximate nearest neighbors for each data point. Instead of calculating exact nearest neighbors, which can be computationally expensive for large datasets, approximate methods are used to find an approximate set of k-NNs for each data point.
The 𝑐–𝑘-NNG construction phase takes two arguments, namely 𝑘, the number of nearest-neighbors to be searched for, and 𝑘𝑐, the factor used by the augmented query algorithm. The variant of the query algorithm increases the time complexity of a single query from 𝑂(log𝑛) to 𝑂(𝑘⋅𝑘𝑐+log𝑛), resulting in an overall time complexity of 𝑂(𝑛(𝑘⋅𝑘𝑐+log𝑛)), where practically 𝑘⋅𝑘𝑐>log𝑛, for the 𝑐–𝑘-NNG construction.

### Calculate MST of the k-NN
<!-- Anchor point for Calculate MST of the k-NN section -->
<a name="calculate-mst"></a>
The MST, or Minimum Spanning Tree, is a graph that connects all data points in a way that minimizes the total edge weight. In this step, the algorithm constructs an MST based on the c-approximate k-NNs obtained in the previous step. The MST helps capture the global structure and relationships among data points.

The algorith used for MST is Kruskal’s algorithm. A greedy algorithm beacuse it makes locally optimal choices at each step (i.e., selecting the minimum-weight edge) with the hope of finding a globally optimal solution (the MST)

![image](https://github.com/afloresep/HMDB-clustering/assets/41540492/23a86368-52ca-416c-8f34-ac6f0152ab2b)


### Lay the tree on the Euclidean plane using the OGDF modular C++ library
<!-- Anchor point for OGDF section -->
<a name="ogdf-library"></a>
To maintain compactness in the visualization and because the Minimum Spanning Tree (MST) is inherently unrooted, a graph layout algorithm, rather than a traditional tree layout, is used. For the purpose of rendering MSTs with a significant number of vertices, often in the millions, a layout algorithm employing a spring-electrical model with a multilevel multipole-based force approximation is employed. This algorithm is made available through the Open Graph Drawing Framework (OGDF), which is a modular C++ library.
