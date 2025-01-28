import argparse
import os
import uuid
from datetime import datetime
import step1_xtc_handling  # Script for step 1: XTC handling and PDB conversion
import step2_pocketscsv    # Script for step 2: Pockets CSV generation
import step3_clustering_conf_creator  # Script for step 3: Clustering config creation


def get_config(output_dir):
    """Generate configuration based on the specified output directory."""

    os.makedirs(output_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") # Timestamp for job ID, actually we did not need it, but could have caused problems if we ran it multiple times.
    job_id = f"{timestamp}_{str(uuid.uuid4())[:4]}"
    
    config = {
        'processed_dir': os.path.join(output_dir, job_id)
    }

    # Create directories if they do not exist
    for key, directory in config.items():
        os.makedirs(directory, exist_ok=True)
        print(f"Created directory: {directory}") if not os.path.exists(directory) else None
    
    return config


def pocket_detection(args, config):
    """Run the full pocket detection pipeline."""
    xtc_file = args.xtc
    topology_file = args.topology
    depth = args.depth
    threads = args.threads
    
    # Convert XTC to PDB
    pdb_files = step1_xtc_handling.xtc_to_pdb(
        xtc_file, 
        topology_file, 
        config['processed_dir']
    )

    # Write PDB list to a file
    pdb_list_file = os.path.join(config['processed_dir'], 'pdb_list.ds')
    step1_xtc_handling.write_pdb_list(config['processed_dir'], pdb_list_file)

    # Run P2Rank on the PDB list
    step1_xtc_handling.run_p2rank(
        pdb_list_file, 
        config['processed_dir'],
        threads=threads
    )

    # Merge results to CSV
    step2_pocketscsv.merge_to_csv(config['processed_dir'], pdb_list_file)

    # Generate clustering heatmap
    step3_clustering_conf_creator.main(
        config['processed_dir'], 
        os.path.join(config['processed_dir'], 'new_pockets.csv'), 
        config['processed_dir'], 
        depth
    )

def convert_to_pdb(args, config):
    """Convert XTC to PDB files with optional threading."""
    xtc_file = args.xtc
    topology_file = args.topology
    threads = args.threads

    # Convert XTC to PDB
    pdb_files = step1_xtc_handling.xtc_to_pdb(
        xtc_file, 
        topology_file, 
        config['processed_dir']
    )


def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="CLI tool for converting trajectory files to PDB format, generating binding site maps, and performing ensemble docking."
    )
    parser.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Pocket detection command
    parser_pocket = subparsers.add_parser("pocket_detection", help="Full pocket detection pipeline")
    parser_pocket.add_argument("--xtc", required=True, help="Input XTC trajectory file")
    parser_pocket.add_argument("--topology", required=True, help="Topology PDB file")
    parser_pocket.add_argument("--depth", type=int, required=True, help="Clustering depth")
    parser_pocket.add_argument("--threads", type=int, default=1, help="Number of parallel threads (default: 4)")
    parser_pocket.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')

    # Convert to PDB command
    parser_convert = subparsers.add_parser("convert_to_pdb", help="Convert XTC to PDB files")
    parser_convert.add_argument("--xtc", required=True, help="Input XTC trajectory file")
    parser_convert.add_argument("--topology", required=True, help="Topology PDB file")
    parser_convert.add_argument("--threads", type=int, default=1, help="Number of parallel threads (default: 1)")
    parser_convert.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')

    args = parser.parse_args()
    
    # Load configuration with user-specified output directory
    config = get_config(args.output)

    # Execute commands
    if args.command == "pocket_detection":
        pocket_detection(args, config)
    elif args.command == "convert_to_pdb":
        convert_to_pdb(args, config)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()