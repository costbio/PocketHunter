# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

PocketHunter is a command-line tool for detecting and characterizing potential small molecule binding pockets in protein molecular simulation trajectories. It identifies both druggable and cryptic binding pockets from XTC trajectory files using p2rank and clustering analysis.

**Important**: This tool is under active development and may contain bugs.

## Setup and Installation

### Initial Setup
```bash
# Install Python dependencies
pip install -r requirements.txt

# Install p2rank (required external dependency)
./first_setup.sh
```

The `first_setup.sh` script downloads and extracts p2rank 2.5 into `tools/p2rank/` within the repository.

### Dependencies
- mdtraj: Trajectory file handling
- pandas, numpy: Data processing
- scikit-learn, scipy: Clustering algorithms
- matplotlib, seaborn: Visualization
- p2rank (external): Pocket prediction engine

## Commands

### Running the Full Pipeline
```bash
python pockethunter.py full_pipeline --xtc path/to/trajectory.xtc --topology path/to/topology.pdb --numthreads 4 --outfolder output_dir --min_prob 0.7 --stride 10
```

### Running Individual Steps

**Extract PDB frames from XTC:**
```bash
python pockethunter.py extract_to_pdb --xtc path/to/trajectory.xtc --topology path/to/topology.pdb --outfolder output_dir --stride 10
```

**Predict pockets (requires extract_to_pdb output):**
```bash
python pockethunter.py detect_pockets --infolder path/to/pdbs --outfolder output_dir --numthreads 4
```

**Cluster pockets (requires detect_pockets output):**
```bash
python pockethunter.py cluster_pockets --infile path/to/pockets.csv --outfolder output_dir --method DBSCAN --min_prob 0.7
```

**Visualize clustered pockets:**
```bash
python pockethunter.py plot_clustermap --infolder path/to/clustered_results
```

### Test Data
Test data is available in `testdata/`:
- `test_data.xtc`: Sample trajectory
- `top.gro`: Sample topology

## Architecture

### Pipeline Flow

The tool implements a sequential pipeline with three core stages:

1. **Extract (extract_predict.py)**:
   - `xtc_to_pdb()`: Converts XTC trajectory frames to individual PDB files using mdtraj
   - `write_pdb_list()`: Creates a dataset file listing all PDB files for p2rank
   - `run_p2rank()`: Executes p2rank binary on PDB files to predict binding pockets

2. **Process (process_predicted.py)**:
   - `merge_to_csv()`: Compiles p2rank predictions from individual CSV files into a single `pockets.csv` with columns: File name, Frame, pocket_index, probability, residues

3. **Cluster (cluster.py)**:
   - Converts residue lists to binary feature vectors (1 if residue in pocket, 0 otherwise)
   - `optimized_dbscan()`: Performs DBSCAN clustering with automatic hyperparameter optimization via silhouette coefficient
   - `get_cluster_medoids()`: Identifies representative pockets using medoid calculation
   - Optional hierarchical sub-clustering within DBSCAN clusters (`--hierarchical` flag)

4. **Visualize (plot.py)**:
   - `plot_clustermap()`: Creates interactive heatmaps showing residue composition per cluster

### Key Design Patterns

**Configuration Management**: A `config` dict containing a logger is passed through all pipeline stages. The logger is initialized in `pockethunter.py` with rotating file handlers.

**Folder Orchestration**: The `full_pipeline()` function dynamically modifies `args.outfolder` and `args.infolder` to create subfolders (`pdbs/`, `pockets/`, `pocket_clusters/`) and chain outputs between stages.

**p2rank Integration**: The tool shells out to an external p2rank binary at `tools/p2rank/prank`. Output is suppressed via `subprocess.DEVNULL`.

**Clustering Strategy**:
- Binary feature vectors represent pockets based on constituent residues
- DBSCAN with Hamming distance for primary clustering
- Silhouette coefficient optimization determines epsilon and min_samples
- Medoids serve as cluster representatives
- Noise points (cluster -1) are tracked separately

### File Structure
```
pockethunter.py          # Main CLI entry point, argument parsing, pipeline orchestration
extract_predict.py       # XTC→PDB conversion, p2rank execution
process_predicted.py     # Merge p2rank CSVs into unified format
cluster.py               # DBSCAN/hierarchical clustering, medoid identification
plot.py                  # Interactive cluster visualization
first_setup.sh           # p2rank installation script
```

### Output Structure
```
output_folder/
├── log_<timestamp>_<uuid>.log
├── pdbs/
│   ├── *.pdb
│   └── pdb_list.ds
├── pockets/
│   ├── pockets.csv                    # Merged predictions
│   └── p2rank_output/                 # Raw p2rank results
└── pocket_clusters/
    ├── pockets_clustered.csv          # All pockets with cluster labels
    ├── cluster_representatives.csv    # Medoid pockets only
    └── cluster_<N>_hierarchical.csv   # Optional hierarchical results
```

## Common Issues

**p2rank not found**: Ensure `first_setup.sh` completed successfully. The binary should exist at `tools/p2rank/prank`.

**No clusters found**: Try lowering `--min_prob` or adjusting the stride to include more frames. DBSCAN requires sufficient data density.

**Memory issues with large trajectories**: Use a larger `--stride` value to reduce the number of extracted frames.
