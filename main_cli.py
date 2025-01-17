import argparse
import os
import step1_xtc_handling  # Script for step 1: XTC handling and PDB conversion
import step2_pocketscsv    # Script for step 2: Pockets CSV generation
import step3_clustering_conf_creator  # Script for step 3: Clustering config creation


def get_config():
    """Placeholder for actual configuration loading logic."""
    config = {
        'inputs_folder': './inputs',
        'pdb_dir': './pdb',
        'processed_dir': './processed',
        'p2rank_processed_dir': './p2rank_processed'
    }

    # Create directories if they do not exist
    for key, directory in config.items():
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")

    return config


def pocket_detection(args, config):
    """Convert XTC to PDB, run P2Rank, and generate a heatmap of representative pockets at the specified depth.""" 
    xtc_file = args.xtc
    topology_file = args.topology
    depth = args.depth
    
    # Convert XTC to PDB
    pdb_files = step1_xtc_handling.xtc_to_pdb(xtc_file, topology_file, config['processed_dir'])

    # Write PDB list to a file
    pdb_list_file = os.path.join(config['processed_dir'], 'pdb_list.ds')
    step1_xtc_handling.write_pdb_list(config['processed_dir'], pdb_list_file)

    # Run P2Rank on the PDB list
    step1_xtc_handling.run_p2rank(pdb_list_file, config['processed_dir'])

    # Path to the PDB list file
    step2_pocketscsv.merge_to_csv(config['processed_dir'],pdb_list_file)

    # Create the hiearchical clustering heatmap with the representative pockets
    step3_clustering_conf_creator.main(config['processed_dir'], os.path.join(config['processed_dir'], 'new_pockets.csv'), config['processed_dir'], depth)

def convert_and_rank(args, config):
    """Convert XTC to PDB and run P2Rank on the generated PDB files."""
    xtc_file = args.xtc
    topology_file = args.topology

    # Convert XTC to PDB
    pdb_files = step1_xtc_handling.xtc_to_pdb(xtc_file, topology_file, config['processed_dir'])

    # Write PDB list to a file
    pdb_list_file = os.path.join(config['processed_dir'], 'pdb_list.ds')
    step1_xtc_handling.write_pdb_list(config['processed_dir'], pdb_list_file)

    # Run P2Rank on the PDB list
    step1_xtc_handling.run_p2rank(pdb_list_file, config['processed_dir'])

    # Path to the PDB list file
    step2_pocketscsv.merge_to_csv(config['processed_dir'],pdb_list_file)

    print(f"P2Rank execution completed. In {config['processed_dir']} directory, you can find the output CSV file with pocket data.")

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="CLI tool for converting trajectory files to PDB format, generating binding site maps, and performing ensemble docking. Execute each step using dedicated commands."
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Command for :::pocket_detection:::  XTC to PDB conversion, P2Rank execution and pocket detection
    parser_pocket_detection = subparsers.add_parser("pocket_detection", help="Converts XTC to PDB, runs P2Rank and detects pockets. Outputs a heatmap of representative pockets.")
    parser_pocket_detection.add_argument("--xtc", type=str, help="Path to XTC file", required=True)
    parser_pocket_detection.add_argument("--topology", type=str, help="Path to topology PDB file", required=True)
    parser_pocket_detection.add_argument("--depth", type=int, help="Clustering depth", required=True)


    # Command for :::convert and rank::: XTC to PDB conversion and P2Rank execution
    parser_convert_and_rank = subparsers.add_parser("convert_and_rank", help="Converts XTC to PDB, runs P2Rank and outputs a summary CSV file of pockets.")
    parser_convert_and_rank.add_argument("--xtc", type=str, help="Path to XTC file", required=True)
    parser_convert_and_rank.add_argument("--topology", type=str, help="Path to topology PDB file", required=True)

    # Parse arguments
    args = parser.parse_args()

    # Load configuration
    config = get_config()

    # Execute commands based on parsed arguments
    if args.command == "pocket_detection":
        pocket_detection(args, config)
    elif args.command == "convert_and_rank":
        convert_and_rank(args, config)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
