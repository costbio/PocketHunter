import argparse
import os
import uuid
from datetime import datetime
import step1_xtc_handling  # Script for step 1: XTC handling and PDB conversion
import step2_pocketscsv    # Script for step 2: Pockets CSV generation
import step3_clustering_conf_creator  # Script for step 3: Clustering config creation
from logging_configuration import logger

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
        if key == 'log_file':
            continue  # Skip log file since it's handled by the logger
        if not os.path.exists(directory):
            os.makedirs(directory)
            logger.info(f"Created directory: {directory}")

    return config


def pocket_detection(args, config):
    try:

        """Run the full pocket detection pipeline."""
        logger.info("Starting the full pocket detection pipeline.")
        xtc_file = args.xtc
        topology_file = args.topology
        clustering_depth = args.clustering_depth
        numthreads = args.numthreads
    
        # Extract frames from XTC to PDB files
        logger.info(f"Extracting frames from XTC file: {xtc_file}")
        pdb_files = step1_xtc_handling.xtc_to_pdb(
            xtc_file, 
            topology_file, 
            config['processed_dir']
        )

        # Write PDB list to a file
        logger.info("Writing PDB list file")
        pdb_list_file = os.path.join(config['processed_dir'], 'pdb_list.ds')
        step1_xtc_handling.write_pdb_list(config['processed_dir'], pdb_list_file)

        # Run P2Rank on the PDB list
        logger.info(f"Running P2Rank with {numthreads} threads")
        step1_xtc_handling.run_p2rank(
            pdb_list_file, 
            config['processed_dir'],
            numthreads=numthreads
        )

        # Merge results to CSV
        logger.info("Merging pocket data to CSV")
        step2_pocketscsv.merge_to_csv(config['processed_dir'], pdb_list_file)

        # Generate clustering heatmap
        logger.info("Creating hierarchical clustering heatmap")
        step3_clustering_conf_creator.main(
            config['processed_dir'], 
            os.path.join(config['processed_dir'], 'new_pockets.csv'), 
            config['processed_dir'], 
            clustering_depth
        )

        logger.info("Full pocket detection pipeline completed successfully.")

    except Exception as e:
        logger.error(f"Error in pocket detection process: {str(e)}")
        raise

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
        description="CLI tool to identify representative small-molecular binding pockets "
         " in protein molecular simulation trajectories using p2rank and via a simple clustering analysis"
    )
    parser.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Pocket detection command
    parser_pocket = subparsers.add_parser("pocket_detection", help="Full pocket detection pipeline")
    parser_pocket.add_argument("--xtc", required=True, help="Input XTC trajectory file")
    parser_pocket.add_argument("--topology", required=True, help="Topology PDB file")
    parser_pocket.add_argument("--clustering_depth", type=int, required=True, help="Clustering depth")
    parser_pocket.add_argument("--numthreads", type=int, default=1, help="Number of parallel threads (default: 4)")
    parser_pocket.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')

    # Convert to PDB command
    parser_convert = subparsers.add_parser("convert_to_pdb", help="Convert XTC to PDB files")
    parser_convert.add_argument("--xtc", required=True, help="Input XTC trajectory file")
    parser_convert.add_argument("--topology", required=True, help="Topology PDB file")
    parser_convert.add_argument("--numthreads", type=int, default=1, help="Number of parallel threads (default: 1)")
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