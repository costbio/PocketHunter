import pandas as pd
import matplotlib.pyplot as plt
import os
import tqdm
import math
import glob
import subprocess
import csv
import seaborn as sns
import scipy.cluster.hierarchy as sch
import numpy as np
import nglview as nv
from IPython.display import display

from scipy.spatial.distance import cdist

from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import hamming
from collections import Counter
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

def optimized_dbscan(df, cols_2cluster, logger):
    """
    Performs optimized DBSCAN on a dataframe to maximize the silhouette coefficient.

    Args:
        df: Pandas DataFrame containing the samples, with features in different columns.
        initial_eps_range: Tuple defining the initial search range for epsilon.
        initial_min_samples_range: Tuple defining the initial search range for min_samples.

    Returns:
        A tuple containing:
            - The best DBSCAN model found.
            - The maximum silhouette coefficient achieved.
            - A dictionary containing the search history (eps, min_samples, silhouette_score).  
            - The scaled dataframe used for clustering
    """

    best_silhouette = -1  # Initialize with a value lower than any possible silhouette score
    best_model = None

    num_cols = len(cols_2cluster)
    logger.info(f'Number of pocket-forming residues: {num_cols}')
    num_frames = len(np.unique(df['Frame'].values))

    # eps_range is defined based on a range of maximum similarity between pockets in the same cluster.
    # accordingly, two pockets may shared up to 1-10 pocket-forming amino acids to be similar.
    # if the difference is higher, then they won't be clustered together.
    # Optimum range is found based on silhouette coefficient determined below.
    eps_range = np.arange(start=1/num_cols, stop=10/num_cols, step=0.005)
    logger.info(f'Number of frames: {num_frames}')

    # min_samples is the minimum number of samples (frames) each cluster must include.
    # this is to filter out irrelevant pockets that might have been observed by chance.
    if num_frames > 100:
        min_samples = math.ceil(num_frames*0.005)
    else:
        min_samples = math.ceil(num_frames*0.02)

    max_samples = math.ceil(num_frames*0.2)

    min_sample_range = np.arange(min_samples, max_samples, max_samples*0.1)

    df.pop('Frame')

    # More efficient search strategy (you can customize this):
    optim_count = 1

    total_optim_count = len(eps_range) * len(min_sample_range)

    for eps in eps_range:  
        for min_sample in min_sample_range:
            min_sample = int(min_sample)

            try:
                model = DBSCAN(eps=eps, min_samples=min_sample, n_jobs=40, metric='hamming')
                labels = model.fit_predict(df[cols_2cluster])

                # Silhouette score requires at least 2 clusters
                n_clusters = len(set(labels)) - (1 if -1 in labels else 0) # account for noise points
                logger.info(f'Optimizing round {optim_count}/{total_optim_count} eps={eps} min_sample={min_sample} num_clusters={n_clusters}')
                if n_clusters >= 2:
                    silhouette = silhouette_score(df, labels)
                    #logger.info(f'Silhouette score: {silhouette}')
                    
                    if silhouette > best_silhouette:
                        best_silhouette = silhouette
                        best_model = model
                        logger.info(best_silhouette)
            except Exception as e: # Catch potential errors (e.g., all points assigned to noise)
                logger.info(f"Error with eps={eps}, min_samples={min_samples}: {e}")

            optim_count += 1

    return best_model, best_silhouette, df

# Function to get medoids
def get_cluster_medoids(data, labels, unique_residues):
    """
    Calculates medoids for each cluster and returns the corresponding rows from the original DataFrame.
    """
    medoid_indices = []
    for cluster_label in set(labels):
        if cluster_label!= -1:
            cluster_points = data[labels == cluster_label][unique_residues]
            # Removed.to_numpy() - hamming can handle Series directly
            distances = np.array([[hamming(point1, point2) for point2 in cluster_points.values] for point1 in cluster_points.values])
            medoid_index = np.argmin(distances.sum(axis=0))
            medoid_indices.append(cluster_points.index[medoid_index])
    return data.loc[medoid_indices]

def hierarchical_clustering(data, method='ward', threshold=1.0):
    """
    Perform hierarchical clustering on the given data.

    Args:
        data: DataFrame containing the data to cluster.
        method: Linkage method to use (default: 'ward').
        threshold: Distance threshold for forming clusters.

    Returns:
        A list of cluster labels.
    """
    linkage_matrix = linkage(data, method=method)
    cluster_labels = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return cluster_labels

def cluster_pockets(infile, outfolder, method, depth, min_prob, config, hierarchical=False):

    logger = config['logger']

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    
    data = pd.read_csv(infile)

    # Filter by min_prob
    data = data[data['probability'] >= min_prob]

    all_residues = []
    for index, row in data.iterrows():
        residues_in_row = row['residues'].split(' ')
        all_residues.extend(residues_in_row)

    # Remove '' from all_residues
    all_residues = [residue for residue in all_residues if residue != '']

    unique_residues = sorted(list(set(all_residues)))  # Get unique residues and sort them for consistency

    # Function to create binary vector
    def sequence_to_binary(sequence, residues_list):
        binary_vector = [1 if res in sequence else 0 for res in residues_list]
        return pd.Series(binary_vector, index=residues_list)

    # Apply the function and create the binary DataFrame
    logger.info('Binarizing the pocket data frame.')
    binary_df = data['residues'].apply(lambda res_str: sequence_to_binary(res_str, unique_residues))
    #binary_df.to_csv(os.path.join(outfolder,'binary_df.csv')) # For debugging.

    # Concatenate with the original DataFrame (optional)
    data = pd.concat([data, binary_df], axis=1)
    
    if method == 'dbscan':
        logger.info('Clustering using the DBSCAN algorithm.')

        # Merge Frame and pocket_index columns to be use it as index
        data['Frame_pocket_index'] = data['Frame'].astype(str) + '_' + data['pocket_index'].astype(str)

        # Set Frame column as the index.
        data.set_index('Frame_pocket_index', inplace=True)

        data_2cluster = data[['Frame']+list(unique_residues)].copy()
        #data_2cluster.to_csv(os.path.join(outfolder,'data_2cluster.csv')) # For debugging.

        best_model, best_silhouette, __ = optimized_dbscan(data_2cluster, unique_residues, logger)

        logger.info(f"Best Silhouette Score: {best_silhouette}")

        # Get the labels from the best model:
        if best_model:
            labels = best_model.labels_
            data['cluster'] = labels # Add cluster labels back to the original dataframe
            data.to_csv(os.path.join(outfolder,'pockets_clustered.csv'))
            cluster_members = {lbl: data[data['cluster'] == lbl] for lbl in np.unique(labels) if lbl != -1}

            if hierarchical:
                logger.info('Performing hierarchical clustering within DBSCAN clusters.')
                for cluster_label, cluster_data in cluster_members.items():
                    logger.info(f'Processing DBSCAN cluster {cluster_label} with {len(cluster_data)} members.')
                    hierarchical_labels = hierarchical_clustering(
                        cluster_data[unique_residues], method='ward', threshold=1.0
                    )
                    cluster_data['hierarchical_cluster'] = hierarchical_labels
                    cluster_data.to_csv(
                        os.path.join(outfolder, f'cluster_{cluster_label}_hierarchical.csv')
                    )
                logger.info('Hierarchical clustering complete.')

            for lbl in np.unique(labels):
                if lbl == -1:
                    percentage_noise = (labels == -1).sum() / len(labels) * 100
                    logger.info(f"Percentage of noise pockets: {percentage_noise:.2f}%")
                else:
                    percentage_cluster = (labels == lbl).sum() / len(labels) * 100
                    logger.info(f"Percentage of cluster {lbl}: {percentage_cluster:.2f}%")
            logger.info('Identifying cluster representatives...')
            medioids = get_cluster_medoids(data, labels, unique_residues)
            logger.info('Identifying cluster representatives... Done.')
            logger.info('Saving cluster representatives to a CSV file...')
            medioids.to_csv(os.path.join(outfolder,'cluster_representatives.csv'))
            logger.info('Clustering complete.')
            
        else:
            logger.info("No suitable DBSCAN model found (no valid clustering achieved).")

