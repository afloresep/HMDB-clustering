## Table of Contents
- [Overview](#overview)
- [Dataset](#dataset)
- [Code](#code)
- [Results](#results)
- [Standardization](#standardization)
- [Contributing](#contributing)
- [License](#license)

## Overview

## Dataset

Data can be found in `MR1-map4-clean.xlsx` and `MR1-map4-clean.csv` datasets. Sourced from HMDB.

## Code

### Original Approach (`original-aproach.ipynb`)

This method unfortunately does not yield the intended results. The reason being that we initially create molecular fingerprints based on similarity between all possible pairs of molecules. We then perform clustering on these pairs, effectively clustering molecules based on the similarity of their pair fingerprints. This approach doesn't capture the molecule's inherent similarity.

### Cluster Fingerprints (`cluster-fp.ipynb`)

Here, we modified the approach by moving away from the MFP_matrix and opting for K-means clustering of the fingerprints. However, this change didn't align with our expectations. The challenge lies in the fingerprint type; the morgan fingerprint we used wasn't optimal for HMDB's diverse dataset containing phospholipids, drug-like molecules, etc. Consequently, despite the code's correctness, the fingerprint failed to provide the needed information for effective clustering.

Addressing the fingerprint, we implemented 'map4' fingerprints from [reymon-group](https://github.com/reymond-group/map4/blob/master/README.md). We tested map4 with data from the Human Metabolome Database. Each of the ~4000 molecules' map4 fingerprints were exported and cleaned into a `MR1_map4-clean.csv` CSV file. Additionally, a smaller dataset of ~100 molecules was created for faster computational times.

We debated the pros and cons of standardizing our data before running PCA and clustering algorithms. This discussion is available in `test.ipynb` or at the end of this file.

Initially, we conducted a test using map4 fingerprints to observe their impact on clustering, detailed in `test.ipynb`. These results displayed higher sensitivity than the previous approach. While not ideal, this outcome is reasonable due to the dataset's chemical diversity. The outcome led us to believe that map4 could be valuable for our clustering problem, prompting us to take a more comprehensive approach in our code.

### Find Optimal K (`test-findK.ipynb`)

This notebook was created to determine the optimal K value for our dataset.

## Standardization

### Standardize or Not to Standardize, That Is the Question {#standardize-or-not-to-standardize-that-is-the-question}

When it comes to PCA, standardizing data is common practice. However, the scenario changes with map4 fingerprints. Their inherent scale (despite variations spanning three orders of magnitude) challenges the necessity of standardization. Standardization is typically employed to bridge dissimilar units, like comparing the height and girth of a tree.

Yet, the HMDB database comprises diverse lipids and phospholipids, potentially explaining magnitude differences. Could standardization address this?

Let's put it to the test and see what unfolds. While PCA often demands standardization, the unique attributes of map4 fingerprints might yield surprising results. While I'm not an expert in statistics, this exploration is enlightening.

In a way, we're comparing apples to apples. Standardization would be akin to contrasting the height and girth of a tree – fundamentally different units.

However, HMDB's wide-ranging content includes various molecules, such as lipids and phospholipids. Could standardization offer a solution to these magnitude variations?

We're giving standardization a shot to observe potential impacts. This experiment might clarify things, or it might not. Regardless, it's a step towards uncovering whether standardization influences our results.

[Jump to the Question](#standardize-or-not-to-standardize-that-is-the-question)
