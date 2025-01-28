# CLI Tool for Trajectory Analysis and Pocket Detection

This repository contains a set of scripts to convert trajectory files to PDB format, generate binding site maps, and perform hierarchical clustering to identify representative pockets. The tool leverages P2Rank for pocket detection and provides visualization capabilities.

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
    git clone https://github.com/bogrum/PocketCLI.git
    cd PocketCLI
    ```

2. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Pocket Detection

Run the full pocket detection pipeline, which includes converting XTC files to PDB format, running P2Rank for pocket detection, and generating a heatmap of representative pockets at the specified clustering depth.

```bash
python main_cli.py pocket_detection --xtc path/to/your.xtc --topology path/to/your_topology.pdb --depth 5 --threads 4 --output path/to/output_dir
```

### Convert to PDB

Convert XTC files to PDB files.

```bash
python main_cli.py convert_to_pdb --xtc path/to/your.xtc --topology path/to/your_topology.pdb --output path/to/output_dir
```

## Configuration

The configuration settings are defined in the `get_config()` function in `main_cli.py`. You can modify the paths and directories as needed.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
