import argparse
import os
import uuid
from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
import shutil
import glob
import subprocess
import csv
import re
import tarfile
import math
from collections import Counter

# External libraries
import mdtraj as md
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import cdist, hamming
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

# ==========================================
# Helpers / Shared Functions
# ==========================================

def setup_logging(log_file):
    """
    Configure and return a logger that writes to both console and file.
    Implements rotating file handler and timestamped log file names.
    """
    logger = logging.getLogger("pockethunter")
    logger.setLevel(logging.DEBUG)

    # Configure file handler with rotating capabilities
    file_handler = RotatingFileHandler(
        log_file,
        maxBytes=10*1024*1024,  # 10MB
        backupCount=5,
        encoding="utf-8"
    )
    file_handler.setLevel(logging.DEBUG)

    # Configure console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Format handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

def get_config(output_dir, overwrite):
    """Generate configuration based on the specified output directory."""

    os.makedirs(output_dir, exist_ok=overwrite)

    # Timestamp for job ID
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_id = f"{timestamp}_{str(uuid.uuid4())[:4]}"

    logger = setup_logging(os.path.join(output_dir,f"log_{job_id}.log"))

    config = {'logger': logger}

    return config

# ==========================================
# Functions from process_predicted.py
# ==========================================

def extract_frame_number(filename):
    """Extracts the number before '.pdb' and after '_' from a filename."""
    match = re.search(r'_(\d+)\.pdb', filename)

    if match:
        return match.group(1)
    else:
        return None

def load_pdb_list(pdb_list_file):
    with open(pdb_list_file, 'r') as f:
        # Assuming each line in pdb_list.ds is a PDB file path
        pdb_files = [line.strip() for line in f.readlines()]
    return pdb_files

def merge_to_csv(outfolder, config):
    logger = config['logger']

    # Define the output CSV filename
    csv_filename = os.path.join(outfolder, "pockets.csv")

    # Ensure the directory exists
    os.makedirs(os.path.dirname(csv_filename), exist_ok=True)

    logger.info('Compiling pocket predictions into summary csv files.')
    # Open CSV file for writing
    with open(csv_filename, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)

        # Write header
        csvwriter.writerow(["File name", 'Frame', "pocket_index", "probability", "residues"])

        # Get a list of all prediction files in the output folder.
        p2rank_output_folder = os.path.join(outfolder, 'p2rank_output')
        prediction_list = glob.glob(p2rank_output_folder+'/*.pdb_predictions.csv')

        # Iterate through PDB files
        for predictions_file in prediction_list:

            # Skip if predictions file doesn't exist
            if not os.path.exists(predictions_file):
                logger.info(f"Predictions file {predictions_file} not found. Skipping.")
                continue

            # Extract the frame number from the PDB file name
            frame_number = extract_frame_number(os.path.basename(predictions_file))
            if not frame_number:
                logger.info(f"Frame number not found for {predictions_file}. Skipping.")
                continue

            # Open predictions CSV file
            with open(predictions_file, 'r') as pred_file:
                reader = csv.DictReader(pred_file)
                for row in reader:
                    pocket_index = row['  rank']  # Pocket index
                    full_name = os.path.basename(predictions_file)[:-4]  # Get the PDB file name as the File name
                    probability = row[' probability']
                    residues = row[' residue_ids']

                    # Write data to the CSV file
                    csvwriter.writerow([full_name, frame_number, pocket_index, probability, residues])
    return csv_filename

# ==========================================
# Functions from extract_predict.py
# ==========================================

def xtc_to_pdb_conversion(xtc_file, topology, stride, outfolder, overwrite, config):
    """Extract frames from XTC file to a series of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting XTC to PDB extraction for {xtc_file}")
        traj = md.load_xtc(xtc_file, topology, stride=stride)

        total_frames = len(traj)
        logger.info(f"Total frames to extract: {total_frames}")

        # Calculate progress reporting intervals (every 10%)
        progress_interval = max(1, total_frames // 10)

        # Ensure output directory exists
        os.makedirs(outfolder, exist_ok=True)
        pdb_files = []
        for i, frame in enumerate(traj):
            real_frame_number = (i+1)*stride
            pdb_file = os.path.join(outfolder, f"{os.path.splitext(os.path.basename(xtc_file))[0]}_{real_frame_number}.pdb")
            frame.save_pdb(pdb_file)
            pdb_files.append(pdb_file)
            # Log progress every 10%
            if (i + 1) % progress_interval == 0 or (i + 1) == total_frames:
                progress_percent = ((i + 1) / total_frames) * 100
                logger.info(f"Extraction progress: {i + 1}/{total_frames} frames ({progress_percent:.1f}%)")

        logger.info(f"Extraction completed successfully. Saved {len(pdb_files)} PDB files.")
        return pdb_files
    except Exception as e:
        logger.error(f"Error during XTC to PDB extraction: {str(e)}")
        raise

def write_pdb_list(pdb_dir, output_file, config):
    try:
        """Write a list of PDB file paths relative to pdb_dir to a file."""
        logger = config['logger']
        logger.info(f"Starting PDB list creation for directory {pdb_dir}")

        pdb_list = glob.glob(os.path.join(pdb_dir, '*.pdb'))
        if not pdb_list:
            logger.warning(f"No PDB files found in {pdb_dir}")
            return None

        with open(output_file, 'w') as f:
            for pdb_file in pdb_list:
                rel_path = os.path.basename(pdb_file)
                f.write(rel_path + '\n')
                logger.debug(f"Added {rel_path} to PDB list")

        logger.info(f"PDB list written to {output_file} containing {len(pdb_list)} files")
        return output_file

    except Exception as e:
        logger.error(f"Error during PDB list creation: {str(e)}")
        raise

def compress_folder(infolder, output_filename, config):
    logger = config['logger']
    logger.info(f"Compressing folder: {infolder} to {output_filename}")
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(infolder, arcname=os.path.basename(infolder))

    # Delete infolder
    logger.info(f"Deleting folder: {infolder}")
    shutil.rmtree(infolder)
    logger.info(f"Compression completed successfully.")

def extract_folder(archive_path, destination_path, config):
    """
    Extracts a folder from a tar.gz archive to a specified destination.

    Args:
        archive_path (str): The path to the tar.gz archive.
        destination_path (str): The path to the destination directory.
        config (dict): A configuration dictionary containing a logger.
    """
    logger = config['logger']
    logger.info(f"Extracting folder from: {archive_path} to {destination_path}")

    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            # Get the top-level folder name from the archive.
            top_level_folder = tar.getnames()[0].split('/')[0]

            # Extract the entire archive to the destination.
            tar.extractall(path=destination_path)

            # Move the extracted top-level folder to the destination.
            extracted_folder_path = os.path.join(destination_path, top_level_folder)
            final_destination = os.path.join(destination_path, os.path.basename(top_level_folder))

            if extracted_folder_path != final_destination: # only move if they are different.
              shutil.move(extracted_folder_path, final_destination)

        logger.info(f"Extraction completed successfully to: {final_destination}")

    except FileNotFoundError:
        logger.error(f"Error: Archive not found at {archive_path}")
    except tarfile.ReadError:
        logger.error(f"Error: Invalid tar.gz archive at {archive_path}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

def run_p2rank(pdb_list_file, output_dir, numthreads, novis, compress, overwrite, config):
    """Run P2Rank with the specified list of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting P2Rank with {numthreads} threads")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        p2rank_output_dir = os.path.join(output_dir, 'p2rank_output')
        os.makedirs(p2rank_output_dir, exist_ok=overwrite)

        command = f'{dir_path}/tools/p2rank/prank predict {pdb_list_file} -o {p2rank_output_dir} -threads {numthreads}'

        if novis == True:
            command += ' -visualizations 0'
        logger.info(f"Executing command: {command}")

        # Redirect stdout and stderr to DEVNULL to suppress output

        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        logger.info(f"P2Rank completed successfully. Output written to {output_dir}")

        # Compiling results into a CSV
        logger.info("Compiling pocket data into a CSV")
        merge_to_csv(output_dir, config)

        if compress:
            compress_folder(p2rank_output_dir, f"{p2rank_output_dir}.tar.gz", config)
            return f"{p2rank_output_dir}.tar.gz"
        else:
            return output_dir

    except Exception as e:
        logger.error(f"Error during P2Rank execution: {str(e)}")
        raise

# ==========================================
# Functions from cluster.py
# ==========================================

def optimized_dbscan(df, cols_2cluster, logger):
    """
    Performs optimized DBSCAN on a dataframe to maximize the silhouette coefficient.
    """
    best_silhouette = -1  # Initialize with a value lower than any possible silhouette score
    best_model = None

    num_cols = len(cols_2cluster)
    logger.info(f'Number of pocket-forming residues: {num_cols}')
    num_frames = len(np.unique(df['Frame'].values))

    # eps_range is defined based on a range of maximum similarity between pockets in the same cluster.
    eps_range = np.arange(start=1/num_cols, stop=10/num_cols, step=0.005)
    logger.info(f'Number of frames: {num_frames}')

    # min_samples is the minimum number of samples (frames) each cluster must include.
    if num_frames > 100:
        min_samples = math.ceil(num_frames*0.005)
    else:
        min_samples = math.ceil(num_frames*0.02)

    max_samples = math.ceil(num_frames*0.2)

    min_sample_range = np.arange(min_samples, max_samples, max_samples*0.1)

    df.pop('Frame')

    # More efficient search strategy
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
    """
    linkage_matrix = linkage(data, method=method)
    cluster_labels = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return cluster_labels

def perform_clustering(infile, outfolder, method, depth, min_prob, config, hierarchical=False):

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

                    cluster_data = cluster_data.copy() #Create a copy to avoid SettingWithCopyWarning (which we don't want to see)
                    cluster_data.loc[:, 'hierarchical_cluster'] = hierarchical_labels
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

# ==========================================
# Functions from plot.py
# ==========================================

def generate_clustermap(infolder, config=None):

    OUTPUT_POCKETS_CSV = os.path.join(infolder, "pocket_clusters/pockets_clustered.csv")

    df_pockets_clustered = pd.read_csv(OUTPUT_POCKETS_CSV)

    # Exclude cluster -1 from the analysis (that is the "noise" cluster)
    df_pockets_clustered = df_pockets_clustered[df_pockets_clustered['cluster'] != -1]

    # --- Parse residues for each frame+cluster ---
    def parse_residues(res_str):
        if pd.isnull(res_str):
            return set()
        return set([r for r in str(res_str).split() if r])

    df_pockets_clustered['residue_set'] = df_pockets_clustered['residues'].apply(parse_residues)

    # Get all unique residues in this simulation
    all_residues = sorted(set.union(*df_pockets_clustered['residue_set']) if len(df_pockets_clustered) > 0 else set())

    # Sort by cluster, then by frame
    df_pockets_clustered.sort_values(['cluster', 'Frame'], inplace=True)

    unique_clusters = df_pockets_clustered['cluster'].unique()
    n_clusters = len(unique_clusters)

    # Calculate number of rows (frames) per cluster
    cluster_sizes = []
    cluster_residues = {}  # Store residues for each cluster
    for clust in unique_clusters:
        cluster_data = df_pockets_clustered[df_pockets_clustered['cluster'] == clust]
        cluster_sizes.append(len(cluster_data))
        # Get unique residues for this cluster
        cluster_residues[clust] = sorted(set.union(*cluster_data['residue_set']) if len(cluster_data) > 0 else set())

    cluster_sizes = np.array(cluster_sizes)
    # Use ratios for height allocation, but cap total height at 14
    min_height = 0.5  # Minimum height per cluster
    total_height = 14
    # Calculate raw ratios
    raw_ratios = cluster_sizes / cluster_sizes.sum()
    # Initial heights
    heights = raw_ratios * total_height
    # Enforce minimum height
    heights = np.maximum(heights, min_height)
    # Renormalize if sum exceeds total_height due to min_height enforcement
    if heights.sum() > total_height:
        heights = heights * (total_height / heights.sum())

    fig, axes = plt.subplots(
        nrows=n_clusters, ncols=1,
        figsize=(18, total_height),
        sharex=True, sharey=False,
        gridspec_kw={'hspace': 0, 'wspace': 0, 'height_ratios': heights}
    )
    if n_clusters == 1:
        axes = [axes]
    cmap_list = sns.color_palette("husl", n_clusters)

    # Store cluster boundaries and data for interaction
    cluster_boundaries = {}
    current_y = 0

    for idx, clust in enumerate(unique_clusters):
        ax = axes[idx]
        clust_df = df_pockets_clustered[df_pockets_clustered['cluster'] == clust]
        index = clust_df['Frame']
        onehot = pd.DataFrame(
            [
                [1 if res in row else 0 for res in all_residues]
                for row in clust_df['residue_set']
            ],
            index=index,
            columns=all_residues
        )

        # Store cluster boundary information
        cluster_boundaries[clust] = {
            'ax': ax,
            'y_start': current_y,
            'y_end': current_y + len(clust_df),
            'residues': cluster_residues[clust],
            'onehot': onehot
        }
        current_y += len(clust_df)

        sns.heatmap(
            onehot,
            cmap=sns.light_palette(cmap_list[idx], as_cmap=True),
            cbar=False,
            ax=ax,
            yticklabels=False,
            xticklabels=False,  # We'll handle x-labels interactively
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("")  # No axis titles

        # Add cluster label to the left of each subplot
        ax.text(
            -0.02, 0.5, f"Cluster {clust}",
            va='center', ha='right',
            fontsize=14, fontweight='bold',
            transform=ax.transAxes,
            rotation=0
        )

        ax.margins(0)

    # Interactive functionality
    class InteractiveClusterMap:
        def __init__(self, fig, axes, cluster_boundaries, all_residues):
            self.fig = fig
            self.axes = axes
            self.cluster_boundaries = cluster_boundaries
            self.all_residues = all_residues
            self.tooltip = None

            # Connect events
            self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)
            self.fig.canvas.mpl_connect('axes_leave_event', self.on_leave)

        def on_hover(self, event):
            if event.inaxes is None:
                self.clear_tooltip()
                return

            # Find which cluster we're hovering over
            for cluster_id, info in self.cluster_boundaries.items():
                if event.inaxes == info['ax']:
                    # Get the heatmap data for this cluster
                    onehot = info['onehot']

                    # Convert mouse coordinates to data coordinates
                    # event.xdata and event.ydata are in data coordinates
                    if event.xdata is not None and event.ydata is not None:
                        x_idx = int(event.xdata)
                        y_idx = int(event.ydata)

                        # Check bounds
                        if 0 <= x_idx < len(self.all_residues) and 0 <= y_idx < len(onehot):
                            residue = self.all_residues[x_idx]

                            # Check if this residue is present in this position (value = 1)
                            if onehot.iloc[y_idx, x_idx] == 1:
                                frame_name = onehot.index[y_idx]
                                self.show_tooltip(event, residue, cluster_id, frame_name)
                            else:
                                self.clear_tooltip()
                        else:
                            self.clear_tooltip()
                    break
            else:
                self.clear_tooltip()

        def on_leave(self, event):
            self.clear_tooltip()

        def show_tooltip(self, event, residue, cluster_id, frame_name):
            # Clear existing tooltip
            self.clear_tooltip()

            # Create tooltip text
            tooltip_text = f"Residue: {residue}\nCluster: {cluster_id}\nFrame: {frame_name}"

            # Create tooltip annotation
            self.tooltip = event.inaxes.annotate(
                tooltip_text,
                xy=(event.xdata, event.ydata),
                xytext=(20, 20),  # Offset from mouse position
                textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.8),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
                fontsize=10,
                zorder=1000
            )

            self.fig.canvas.draw_idle()

        def clear_tooltip(self):
            if self.tooltip:
                self.tooltip.remove()
                self.tooltip = None
                self.fig.canvas.draw_idle()

    # Create interactive handler
    interactive_handler = InteractiveClusterMap(fig, axes, cluster_boundaries, all_residues)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.suptitle(f"Residue Composition per Cluster (Hover over heatmap points to see residue details)", y=1.01)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()

# ==========================================
# Main Execution / Subcommand Wrappers
# ==========================================

def full_pipeline(args, config):
    try:
        """Run the full pipeline."""
        logger = config['logger']
        logger.info("Starting the full pocket detection pipeline.")
        outfolder = args.outfolder

        # Modify args.outfolder temporarily to enable extract_to_pdb to write to a subfolder inside main outfolder here.
        args.outfolder = os.path.join(outfolder,'pdbs')

        # Extract frames from XTC to PDB files
        extract_to_pdb(args, config)

        # Add args.infolder needed by detect_pockets.
        args.infolder = args.outfolder

        # Modify args.outfolder temporarily again for detect_pockets to write a subfolder.
        args.outfolder = os.path.join(outfolder,'pockets')

        # Detect pockets on PDB files
        detect_pockets(args, config)

        # Modify args.infolder for cluster_pockets again.
        args.infolder = args.outfolder

        # Set args.infile to the pockets.csv file generated by detect_pockets
        args.infile = os.path.join(args.infolder, "pockets.csv")

        # Modify args.outfolder again
        args.outfolder = os.path.join(outfolder,'pocket_clusters')

        # Perform clustering, identify representatives.
        cluster_pockets(args, config)

        logger.info("Full pocket detection pipeline completed successfully.")

    except Exception as e:
        logger.error(f"Error in pocket detection process: {str(e)}")
        raise

def extract_to_pdb(args, config):
    """Extract frames from XTC to PDB files."""
    xtc_file = args.xtc
    topology_file = args.topology
    stride = args.stride
    outfolder = args.outfolder
    overwrite = args.overwrite
    logger = config['logger']
    logger.info(f"Extracting frames from XTC file: {xtc_file}")

    # Convert XTC to PDB
    pdb_files = xtc_to_pdb_conversion(xtc_file, topology_file,
        stride, outfolder, overwrite, config)

    # Write PDB list to a file
    logger.info("Writing PDB list file")
    pdb_list_file = os.path.join(outfolder, 'pdb_list.ds')
    write_pdb_list(outfolder, pdb_list_file, config)

def detect_pockets(args, config):
    """Detect pockets on PDB files in input folder using p2rank."""
    infolder = os.path.abspath(args.infolder)
    outfolder = args.outfolder
    numthreads = args.numthreads
    logger = config['logger']
    pdb_list_file = os.path.join(infolder,'pdb_list.ds')
    novis = args.novis
    compress = args.novis
    overwrite = args.overwrite

    # Run P2Rank on the PDB list
    run_p2rank(pdb_list_file, outfolder, numthreads, novis, compress, overwrite, config)

def cluster_pockets(args, config):
    """Identifies representative pockets via clustering."""
    infile = os.path.abspath(args.infile)
    outfolder = args.outfolder
    method = args.method
    depth = args.depth
    min_prob = args.min_prob

    perform_clustering(
        infile=args.infile,
        outfolder=args.outfolder,
        method=args.method,
        depth=args.depth,
        min_prob=args.min_prob,
        config=config,
        hierarchical=args.hierarchical
    )

def plot_clustermap(args):
    """Generate cluster map visualization from clustering results."""
    infolder = os.path.abspath(args.infolder)

    generate_clustermap(
        infolder=infolder
    )

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="CLI tool to identify representative small-molecular binding pockets "
         " in protein molecular simulation trajectories using p2rank and via a simple clustering analysis"
    )
    parser.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Full pipeline
    parser_full = subparsers.add_parser("full_pipeline", help="Full pipeline for pocket prediction.")
    parser_full.add_argument("--xtc", required=True, help="Input XTC trajectory file.")
    parser_full.add_argument("--topology", required=True, help="Topology PDB/GRO file.")
    parser_full.add_argument("--depth", type=int, default=4, help="Clustering depth.")
    parser_full.add_argument("--method", type=str,  default="dbscan", choices=["hierarchical","dbscan"],
                                help='Clustering method. Choose between hierarchical and dbscan.')
    parser_full.add_argument("--min_prob",type=float, default=0.5, help="Minimum ligand-binding probability for clustering (default: 0.5).")
    parser_full.add_argument('--novis', action='store_true', help='Do not keep visualizations generated by p2rank.')
    parser_full.add_argument('--compress', action='store_true', help='Compress p2rank output directory.')
    parser_full.add_argument("--numthreads", type=int, default=4, help="Number of parallel threads (default: 4)."),
    parser_full.add_argument("--stride", type=int, default=1, help="Stride applied for frame extraction."),
    parser_full.add_argument("--overwrite", action='store_true', help="Overwrite existing output directory.")
    parser_full.add_argument('--outfolder', type=str, default='.', help='Output folder.')
    parser_full.add_argument('--hierarchical', action='store_true', help='Perform hierarchical clustering within DBSCAN clusters.')


    # Extract frames from XTC
    parser_convert = subparsers.add_parser("extract_to_pdb", help="Extract frames from XTC to PDB files.")
    parser_convert.add_argument("--xtc", required=True, help="Input XTC trajectory file.")
    parser_convert.add_argument("--topology", required=True, help="Topology PDB/GRO file.")
    parser_convert.add_argument("--stride", type=int, default=1, help="Stride applied for frame extraction.")
    parser_convert.add_argument('--outfolder', type=str, help='Output folder.')
    parser_convert.add_argument('--overwrite', action='store_true', help='Overwrite existing output directory.')

    # Predict pockets
    parser_detect = subparsers.add_parser("detect_pockets", help="Detect pockets on PDB files in input folder using p2rank.")
    parser_detect.add_argument("--numthreads", type=int, default=4, help="Number of parallel threads (default: 4).")
    parser_detect.add_argument('--infolder', type=str, help='Input folder containing PDB files (should be generated by extract_to_pdb).')
    parser_detect.add_argument('--outfolder', type=str, help='Output folder.')
    parser_detect.add_argument('--novis', action='store_true', help='Do not keep visualizations generated by p2rank.')
    parser_detect.add_argument('--compress', action='store_true', help='Compress p2rank output directory.')
    parser_detect.add_argument('--overwrite', action='store_true', help='Overwrite existing output directory.')

    # Cluster pockets
    parser_cluster = subparsers.add_parser("cluster_pockets",help="Clusters identified pockets based on binary feature vectors "
                                           "denoting residue presence in the pocket.")
    parser_cluster.add_argument("--infile", type=str, help='Input CSV file produced by predict_pockets (pockets.csv).')
    parser_cluster.add_argument("--outfolder", type=str, help='Output folder.')
    parser_cluster.add_argument("--method", type=str, default="dbscan", choices=["hierarchical","dbscan"],
                                help='Clustering method. Choose between hierarchical and dbscan.')
    parser_cluster.add_argument("--depth", type=int, default=4,
                                help='Clustering depth parameter for identification of representatives pockets from '
                                'hierarchical clustering.')
    parser_cluster.add_argument("--min_prob",type=float, default=0.5, help="Minimum ligand-binding probability for clustering (default: 0.5).")
    parser_cluster.add_argument('--hierarchical', action='store_true', help='Perform hierarchical clustering within DBSCAN clusters.')
    parser_cluster.add_argument('--overwrite', action='store_true', help='Overwrite existing output directory.')

    # Plot cluster map
    parser_plot = subparsers.add_parser("plot_clustermap", help="Generate cluster map visualization from clustering results.")
    parser_plot.add_argument("--infolder", type=str, help='Input folder containing clustering results (output from cluster_pockets).')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return  # Exit without running anything

    # Load configuration with user-specified output directory
    if not args.command == "plot_clustermap":
        config = get_config(args.outfolder, args.overwrite)

    # Execute commands
    if args.command == "full_pipeline":
        full_pipeline(args, config)
    elif args.command == "extract_to_pdb":
        extract_to_pdb(args, config)
    elif args.command == "detect_pockets":
        detect_pockets(args, config)
    elif args.command == "cluster_pockets":
        cluster_pockets(args, config)
    elif args.command == "plot_clustermap":
        plot_clustermap(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
