import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances, cosine_similarity
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.multivariate.pca import PCA
from sklearn import manifold
from sklearn.decomposition import PCA as skPCA

def calculate_aa_composition(sequence):
    """
    Calculate amino acid composition percentages for a protein sequence
    """
    AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
           'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    
    aa_counts = [sequence.count(aa) for aa in AAs]
    total = sum(aa_counts)
    aa_percents = pd.Series(aa_counts)/total * 100
    return aa_percents

def engineer_features(df):
    """
    Standardize features by scaling each column
    """
    df2 = df.copy()
    
    for col in df.columns:
        col_i = df2[col]
        col_i_mean = col_i.mean()
        col_i_std = col_i.std()
        df2[col] = (col_i - col_i_mean)/col_i_std
    
    return df2.dropna(axis=1)

def eigen_scaling(pca, scaling=0):
    """
    Scale eigenvalues for PCA visualization
    """
    const = ((pca.scores.shape[0]-1) * pca.eigenvals.sum()/ pca.scores.shape[0])**0.25
    
    if scaling == 0:
        scores = pca.scores
        loadings = pca.loadings
    elif scaling in [1, 2, 3]:
        scaling_fac = (pca.eigenvals / pca.eigenvals.sum())**(0.5 if scaling in [1, 2] else 0.25)
        scaling_fac.index = pca.scores.columns
        
        scores = pca.scores * (scaling_fac if scaling == 1 else 1) * const
        loadings = pca.loadings * (scaling_fac if scaling in [2, 3] else 1) * const
    
    return scores, loadings

def create_biplot(pca, x="comp_0", y="comp_1", scaling=2, color=None):
    """
    Create PCA biplot with scores and loadings
    """
    scores, loadings = eigen_scaling(pca, scaling=scaling)
    
    scores_copy = scores.copy()
    if color is not None:
        scores_copy["group"] = color
    
    sns.relplot(x=x, y=y, data=scores_copy, hue="group" if color else None,
                palette="muted", alpha=1)
    
    # Plot loading vectors
    for i in range(loadings.shape[0]):
        plt.arrow(0, 0, loadings.iloc[i, 0], loadings.iloc[i, 1],
                 color='black', alpha=0.7)
        plt.text(loadings.iloc[i, 0]*1.05, loadings.iloc[i, 1]*1.05,
                loadings.index[i])

def calculate_distance_matrices(df):
    """
    Calculate Euclidean and cosine distance matrices
    """
    euclidean_mat = pd.DataFrame(
        euclidean_distances(df),
        index=df.index,
        columns=df.index
    )
    
    cosine_mat = pd.DataFrame(
        1 - cosine_similarity(df),
        index=df.index,
        columns=df.index
    )
    
    return euclidean_mat, cosine_mat

def perform_nmds(distance_matrix):
    """
    Perform Non-metric Multidimensional Scaling
    """
    seed = np.random.RandomState(seed=3)
    nmds = manifold.MDS(
        n_components=2,
        metric=False,
        max_iter=3000,
        eps=1e-12,
        dissimilarity="precomputed",
        random_state=seed
    )
    
    positions = nmds.fit_transform(distance_matrix)
    clf = skPCA(n_components=2)
    positions = clf.fit_transform(positions)
    
    return pd.DataFrame(positions, columns=["dim1", "dim2"])

def plot_heatmap(matrix, title="Heatmap"):
    """
    Create clustered heatmap
    """
    sns.clustermap(matrix, 
                  cmap="RdBu",
                  annot=True,
                  fmt=".2f")
    plt.title(title)

if __name__ == "__main__":
    # Example workflow
    protein_sequences = {
        "PAK3": "MSDGLDNEEKPPAPPLRMNSN...",
        "CAPN6": "MGPPLKLFKNQKYQELKQEC...",
        "CHRDL1": "MRKKWILEDFHFIFFGVLC..."
    }
    
    # Calculate compositions
    compositions = pd.DataFrame({name: calculate_aa_composition(seq) 
                               for name, seq in protein_sequences.items()}).T
    
    # Engineer features
    engineered_features = engineer_features(compositions)
    
    # Perform PCA
    pca = PCA(data=engineered_features, method="svd", 
              standardize=True, normalize=True, ncomp=2)
    
    # Calculate distance matrices
    euclidean_mat, cosine_mat = calculate_distance_matrices(engineered_features)
    
    # Generate visualizations
    create_biplot(pca, color=engineered_features.index)
    plot_heatmap(euclidean_mat, "Euclidean Distance Matrix")
    plot_heatmap(cosine_mat, "Cosine Similarity Matrix")
    
    # Perform NMDS
    nmds_positions = perform_nmds(euclidean_mat)
