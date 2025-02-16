# CLI Tool for Trajectory Analysis and Pocket Detection

PocketHunter is a command-line tool for analysis of protein molecular simulation trajectories prior to a docking-based virtual screening study. 

The tool can be used to:

* characterize druggable small-molecule binding pockets using p2rank from protein molecular simulation trajectories in XTC format. 
* identify potential cryptic binding pockets
* select ligand-receptive pockets for (ensemble) docking-based virtual screening based on their performance in separating actives from decoys/inactives in a given input molecule set.

The tool first identifies all small-molecule binding pockets in all input conformations. Then, performs a DBSCAN clustering analysis based on amino acids forming each pocket, using amino acid memberships as binary features for each pocket. Final pocket clusters identified denote various conformations of major binding pockets observed throughout the simulation.

The user can optionally select one of the clusters for active/inactive separation performance analysis. In this case, an annotated list of SMILES strings (indicating whether the molecule is active/inactive) should be provided. The tool can take this list as input, perform a docking screen of representatives of the selected cluster pockets, and report their separation performance in terms of AUROC values.


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Pocket Detection](#pocket-detection)
  - [Convert to PDB](#convert-to-pdb)
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

Pocket Hunter can be used in two main steps, namely identify and select pockets, respectively.

### Identify pockets

First, pockets in molecular simulation trajectory is detected, and pocket clusters are identifies.

```bash
python main_cli.py full_pipeline --xtc path/to/your.xtc --topology path/to/your_topology.pdb --numthreads 4 --output path/to/output_dir --min_prob 0.7 --stride 10
```

### Select pockets

Then, the inhibitor identification performance of pockets in a selected cluster can be tested on a user-provided list of actives/inactives. The cluster analyzed here is often a well-characterized active site of the protein, with previously available ligand-binding information.

```bash
python main_cli.py command pending....
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
