## Data 

Data can be found in `MR1-map4-clean.xlsx` and `MR1-map4-clean.csv` datasets.



## Code

`original-aproach.ipynb` - This method does not work as intended. Reason being that we first create molecular fingerprints based on similiraty between all pairs possible for all the molecules and then do clustering on them. That means that we are clustering molecules based on data from for how similar are the pair fingerprint instead of how similar is the molecule itself. 

`cluster-fp.ipynb` - Here we removed the MFP_matrix that was calculating similiraty between fingerprints pairs and went for a K-means clustering of the fingerprints. The result did not match our expectations. The reason is that the morgan fingerprint that we were using was not ideal for a diverse dataset such as HMDB (with phospholipids, drug-like molecules etc.). Therefore, even tough the code was fine, the fingerprint could not add as much information as we needed for clustering molecules. Hence, we had very different molecules grouped in the same cluster.

Being the fingerprint the main issue and not the code itself, we implemented 'map4' fingerprints from [reymon-group](https://github.com/reymond-group/map4/blob/master/README.md). Map4 was tested with data from the Human Metabolome data base. 
The map4 fingerprint for each of the ~4000 molecules was exported and cleaned into a csv file named `MR1_map4-clean.csv`. Also a small dataset with only ~100 molecules was created for faster computational times. 

We first discussed the possibility of standarizing our data before running the PCA and clustering algorithms. The discussion for this can be seen in `test.ipynb` or in the bottom of this file. 

[Here](#standardize-or-not-to-standardize-that-is-the-question)


First we did a test for seeing how the map4 fp helped us clustering our molecules that can be seen in `test.ipynb`. The results were much more sensitive than the first ones. A great part of our dataset was grouped into the same cluster while way smaller groups formed in different clusters. Again, this is not ideal (only one centroid for half of our datasets and several other centroids for very small data portions) but it makes sense when dealing with such a chemically diverse group of metabolits. This shows that map4 should be useful for our clustering problems so we decided to go for a more in-depth approach for our code. 

`test-findK.ipynb` was created to calculate the best K for our dataset. 



### Standardize or not to standardize, that is the question {#standardize-or-not-to-standardize-that-is-the-question}


Typically, data is standardized before performing a PCA. 
However, in the case of map4 fingerprints, I find that standarization might not be necessary due to being messured in the same scale (despite variations spanning three orders of magnitud). 
One would use standarization when contrasting disparate units like the height of a tree and the girth of the same tree. 

Altough HMDB database contains diverse lipids, phospholipids... which could explain the difference in magnitud for some fingerprints? So maybe standarization is due?

Let's try doing it and see what happens...

When it comes to using PCA, data standarization is usually the go-to step. However, the situation is a bit different with map4 fingerprints. Despite having variations that span three orders of magnitude (1e05 / 1e08) these fingerprints are actually measured on a consistent scale (map4). This makes me wonder if standarization is really necessary...
Statistics is not my area of expertise so i could be wrong!

We're comparing apples to apples so to speak... For instance, you'd standardize when comparing something like the height of a tree and the tree's girth – they're completely different units.

However... the HMDB database contains all sorts of stuff, lipids, phospholipids... This variety might be the reason behind the magnitude differences we're seeing in some fingerprints. Could standarization help with that??

We'll give standarization a shot and see if it changes things. Maybe it'll clear things up, maybe not. Either way, we'll get a better sense of whether standarization plays a role or not

[Jump to the Question](#standardize-or-not-to-standardize-that-is-the-question)
