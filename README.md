# CLI Tool for Trajectory Analysis and Pocket Detection

This repository contains a set of scripts to convert trajectory files to PDB format, generate binding site maps, and perform hierarchical clustering to identify representative pockets. The tool leverages P2Rank for pocket detection and provides visualization capabilities.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Convert and Rank](#convert-and-rank)
  - [Pocket Detection](#pocket-detection)
- [Scripts Overview](#scripts-overview)
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

### Convert and Rank

Convert XTC files to PDB format and run P2Rank on the generated PDB files.

```bash
python main_cli.py convert_and_rank --xtc path/to/your.xtc --topology path/to/your_topology.pdb
```

### Pocket Detection

Convert XTC files to PDB format, run P2Rank, and generate a heatmap of representative pockets at the specified clustering depth.

```bash
python main_cli.py pocket_detection --xtc path/to/your.xtc --topology path/to/your_topology.pdb --depth 5
```
## Configuration

The configuration settings are defined in the `get_config()` function in `main_cli.py`. You can modify the paths and directories as needed.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
