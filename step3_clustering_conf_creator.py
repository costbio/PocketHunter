import pandas as pd
import matplotlib.pyplot as plt
import os
import tqdm
import glob
import subprocess
import csv
import seaborn as sns
import scipy.cluster.hierarchy as sch
import numpy as np
import nglview as nv
from IPython.display import display

from scipy.spatial.distance import cdist

from scipy.cluster import hierarchy
from sklearn.preprocessing import MultiLabelBinarizer



def cluster_pockets(p2rank_output_folder,pdb_location, depth, plot=True):

    csv_file_dir = f"{p2rank_output_folder}/pockets.csv"

    # Step 1: Read the CSV file
    data = pd.read_csv(csv_file_dir)

    residues = data['residues'].str.split()
    mlb = MultiLabelBinarizer()
    residues_encoded = mlb.fit_transform(residues)
    residue_list = mlb.classes_
    for i in range(0,len(residue_list)):
        res = residue_list[i]
        data[res] = residues_encoded[:,i]
    data['frame_pocket'] = ['Frame:{} Index:{}'.format(f, p) for f, p in zip(map(str, data['Frame']), map(str, data['pocket_index']))]
    data.index = data['frame_pocket']
    dataindex = data['frame_pocket']
    data_2cluster = data.drop(columns=['File name','Frame','pocket_index','probability','residues','frame_pocket'])
    linkage = hierarchy.linkage(data_2cluster, method='ward')
    full_heatmap = sns.clustermap(data_2cluster, method='ward', cmap='viridis', row_cluster=True, col_cluster=False)
    if plot:
        plt.show()

    # Get the clusters at the chosen depth
    # Depth is the value to choose the depth you want to explore
    cutree = hierarchy.cut_tree(linkage, n_clusters=depth)
    # Function to calculate Euclidean distance
    def euclidean_distance(p, q):
        return np.sqrt(np.sum((p - q)**2))
    # Get representatives from each cluster at the chosen depth
    representatives = []

    list_rep_pockets = list()

    for cluster_id in range(depth):
        indices = [i for i, c in enumerate(cutree) if c[0] == cluster_id]
        if indices:
            cluster_points = residues_encoded[indices]
            #print(cluster_points)
            centroid = np.mean(cluster_points, axis=0)  # Calculate centroid/mean of the cluster
            # Find euclidean distance of each cluster point to the centroid
            distances = [euclidean_distance(centroid, point) for point in cluster_points]
            # Find the point closest to the centroid
            representative_index = np.argmin(distances)

            representative = cluster_points[representative_index]
            # Find the row in the original data frame that has this representative_index
            representative_pockets = data.iloc[indices[representative_index]]

            list_rep_pockets.append(representative_pockets)

            # File the representative pocket and the cluster it belongs to
            representatives.append((cluster_id, representative_pockets))

            #print(representative_pockets[4], type(representative_pockets))

    #For heatmap, we need a list of amino acids. Let's extract.

    filenames = list(data["File name"]) # List of the files in dataframe
    aminoacids = {} #New dict to keep new data

    for file in filenames:
        path = os.path.join(p2rank_output_folder,f"{file}.pdb_residues.csv")
        print(f"Adding for {path}...")
        residue_file = pd.read_csv(path) # current conformation's resnum=resid data is read
        residue_file.columns = residue_file.columns.str.strip()
        for _, row in residue_file.iterrows():
            aminoacids[row['residue_label']] = row['residue_name']
    
    aminoacids = dict(sorted(aminoacids.items()))

    residue_names = list(aminoacids.values())

    # After dropping the unnecessary columns
    data_df_for_aminoacids = data.drop(columns=['residues', 'frame_pocket', 'probability', 'File name', 'Frame', 'pocket_index'])

    # Extract the list of remaining columns (residue numbers)
    residue_columns = data_df_for_aminoacids.columns.tolist()

    # Print the list of columns
    #print(residue_columns)

    # Convert the target amino acids to a set for faster lookup
    target_set = set([int(x.split('_')[1]) for x in residue_columns])

    aminoacid_list = {key: value for key, value in aminoacids.items() if key in target_set}

    aminoacid_list = list(aminoacid_list.values())

    aminoacid_list = [item.strip() for item in aminoacid_list]

    amino_acid_map = {
    'MET': 'M', 'ASP': 'D', 'PRO': 'P', 'GLU': 'E', 'LEU': 'L', 'HIS': 'H',
    'CYS': 'C', 'ARG': 'R', 'TYR': 'Y', 'GLY': 'G', 'ALA': 'A', 'PHE': 'F',
    'GLN': 'Q', 'VAL': 'V', 'LYS': 'K', 'SER': 'S', 'ILE': 'I', 'THR': 'T',
    'ASN': 'N', 'TRP': 'W'}

    one_letter_amino_acids = [amino_acid_map[aa] for aa in aminoacid_list]

    ylabel = []
    for i in range(len(one_letter_amino_acids)):
        x = f"{residue_columns[i]}_{one_letter_amino_acids[i]}"
        ylabel.append(x)

    #print(new_dict)




    df_rep_pockets = pd.DataFrame(list_rep_pockets)

    df_rep_pockets['frame_pocket'] = ['Frame:{} Index:{}'.format(f, p) for f, p in zip(map(str, df_rep_pockets['Frame']), map(str, df_rep_pockets['pocket_index']))]
    df_rep_pockets.index = df_rep_pockets['frame_pocket']
    repdataindex = df_rep_pockets['frame_pocket']
    data_2cluster = df_rep_pockets.drop(columns=['File name','Frame','pocket_index','probability','residues','frame_pocket'])
    linkage = hierarchy.linkage(data_2cluster, method='ward')
    rep_heatmap = sns.clustermap(data_2cluster, method='ward', cmap='viridis', row_cluster=True, col_cluster=False)
    if plot:
        plt.show()


    for rep in representatives:
        filename = rep[1]['File name']
        file_path = os.path.join(pdb_location, f"{filename}.pdb")
        print(f"Cluster {filename}:")
        print(f"File name: {file_path}")
        print(rep[1]['residues'])

        if not os.path.isfile(file_path):
            print(f"Warning: File does not exist - {pdb_location}/{file_path}.pdb")
            continue

        if plot:
            try:
                view = nv.show_structure_file(file_path)
                view.add_cartoon()
                resnums = [el[2:] for el in rep[1]['residues'].split(' ') if el.startswith('P_')]
                resnums_nglview_sel = ' or '.join(resnums)
                view.add_surface(f':P and ({resnums_nglview_sel})')
                display(view)
            except Exception as e:
                print(f"Error loading file {file_path}: {e}")


    return df_rep_pockets, full_heatmap, rep_heatmap, dataindex, repdataindex, ylabel


def main(csv_file_path, output_file_path, pdb_location, clustering_depth, plot=True):
    """
    Main function to handle PDB file processing and clustering analysis.
    
    Parameters:
        csv_file_path (str): Path to the input CSV file.
        output_file_path (str): Path to save the output file.
        clustering_depth (int): Number of clusters for hierarchical clustering.
        pdb_location (str): Path to the pdb files' folder.
        plot (bool): Whether to plot the heatmaps.
    """
    try:
        # Ensure the required file exists
        if not os.path.exists(csv_file_path):
            raise FileNotFoundError(f"CSV file not found: {csv_file_path}")

        print("Starting clustering analysis...")
        
        # Perform clustering analysis on P2Rank output
        df_rep_pockets, full_heatmap, rep_heatmap, dataindex, repdataindex, ylabel = cluster_pockets(csv_file_path, pdb_location, clustering_depth, plot)

        # Save heatmaps as PNG files
        full_heatmap_path = os.path.splitext(output_file_path)[0] + "_full_heatmap.png"
        rep_heatmap_path = os.path.splitext(output_file_path)[0] + "_rep_heatmap.png"
        
        full_heatmap.savefig(full_heatmap_path)
        rep_heatmap.savefig(rep_heatmap_path)
        
        print(f"Full heatmap saved to {full_heatmap_path}")
        print(f"Representative heatmap saved to {rep_heatmap_path}")
    
    except Exception as e:
        print(f"An error occurred: {e}")

