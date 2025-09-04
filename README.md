# PocketHunter: A command-line tool for pocket detection in biomolecular simulation trajectories

**IMPORTANT:** This tool is currently under development, and may contain serious bugs. Use it at your own risk.

PocketHunter is a command-line tool for detection and characterization of potential small molecule binding pockets in protein molecular simulation trajectories. 

PocketHunter can be used to:

* characterize druggable small-molecule binding pockets using p2rank from protein molecular simulation trajectories in XTC format. 
* identify potential cryptic binding pockets.


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [License](#license)

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/costbio/PocketHunter.git
    cd PocketHunter
    ```

2. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

The tools performs the analysis in several sequential steps. 

* Extract PDBs
* Predict pockets
* Cluster pockets


For flexibility, the full pipeline can be run all at once, or by running the steps individually, then feeding results from one step to the next one.

### Full pipeline

The whole pipeline can be run with a single command as follows:

```bash
python pockethunter.py full_pipeline --xtc path/to/your.xtc --topology path/to/your_topology.pdb --numthreads 4 --outfolder path/to/output_dir --min_prob 0.7 --stride 10 
```

### Extract PDBs 

```bash
python pockethunter.py extract_pdbs --xtc path/to/your.xtc --topology path/to/your_topology.pdb --outfolder path/to/output_dir --stride 10
```
### Predict pockets 

**--infolder** argument below should be the output folder from the **extract_pdbs** subcommand.

```bash
python pockethunter.py predict_pockets --infolder path/to/input_dir --output path/to/output_dir --numthreads 4
```

### Cluster pockets 

**--infolder** argument below should be the output folder from the **predict_pockets** subcommand.

```bash
python pockethunter.py cluster_pockets --infolder path/to/input_dir --outfolder path/to/output_dir --method DBSCAN --min_prob 0.7
```

### Visualize clustered pockets

After clustering, you can visualize the residue composition of each pocket cluster using `plot.py`:

```bash
python pockethunter.py plot_clustermap /path/to/clustered/results_folder
```

This will open an interactive heatmap showing residue composition per cluster. Hover over points to see details.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
