## Goals
- What type?
	- Low-dimensional representation
	- With different input feature types
- Interpretation
	- Input feature variance per latent dimension (HIGH)
		- Prior knowledge of combinations of features (e.g. pathways) (MEDIUM)
	- Correlation between metadata variable and latent dimension (HIGH)
	- Metadata prediction based on multivariate latent space (LOW)
	- Clustering (MEDIUM)
		- Which clusters are variable in which latent dimensions
		- Differential features between clusters
- QC of latent spaces
	- Distribution check of the latent dimensions (MEDIUM)
		- Is normalization needed?
	- Robustness (HIGH/MEDIUM)
		- For changes in input features
		- For changes in latent space
	- Correlation (HIGH)
		- Between latent dimensions
	- Clustering metrics (LOW)
		- Sillhouette scores
		- Jaccard Index
	- Coranking (LOW)
- Latent space comparisons (Too big of a problem)
- Visualization
	- Static 
	- Interactive
- Misc ideas
	- Pseudotime trajectory for images
## Design
- Python
- PyTorch
- Github
	- Codename: TBA
### Input
- Input training data
	- csv.gz matrix
	- Scanpy
- Model
	- PyTorch model 
- Latent space
	- observations x features
	- numpy matrix
- Metadata
	- observations x metadata variables
	- dataframe
### Output
- Example notebook
### Documentation
- Numpy docstrings

## Example data
- 2 folders:
	- Subcellular expression patterns
	- Multimodal prostate cancer dataset

## Work plan
- David: Boilerplate code:
	- Data loading 
	- Github structure
- Nacho: 
	- Input feature variance per latent dimension 
- Gabriele: 
	- Distribution check of the latent dimensions
- Ceyhun: 
	- Which clusters are variable in which latent dimensions
	- Differential features between clusters
