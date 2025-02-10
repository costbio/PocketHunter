import argparse
import os
import uuid
from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
import step1_xtc_handling  
import step2_pocketscsv    
import step3_clustering_conf_creator  

def setup_logging(log_file):
    """
    Configure and return a logger that writes to both console and file.
    Implements rotating file handler and timestamped log file names.
    """
    logger = logging.getLogger("main_logger")
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
def get_config(output_dir):
    """Generate configuration based on the specified output directory."""

    os.makedirs(output_dir, exist_ok=False)

    # Timestamp for job ID
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") 
    job_id = f"{timestamp}_{str(uuid.uuid4())[:4]}"

    logger = setup_logging(os.path.join(output_dir,f"log_{timestamp}_{job_id}"))
    
    config = {'logger': logger}

    return config


def full_pipeline(args, config):
    try:

        """Run the full pipeline."""
        logger = config['logger']
        logger.info("Starting the full pocket detection pipeline.")
        xtc_file = args.xtc
        topology_file = args.topology
        clustering_depth = args.clustering_depth
        numthreads = args.numthreads
        stride = args.stride
        outfolder = args.outfolder
    
        # Extract frames from XTC to PDB files
        extract_to_pdb(args, config)

        # Detect pockets on PDB files
        detect_pockets(args, config)

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

def extract_to_pdb(args, config):
    """Extract frames from XTC to PDB files."""
    xtc_file = args.xtc
    topology_file = args.topology
    stride = args.stride
    outfolder = args.outfolder
    logger = config['logger']
    logger.info(f"Extracting frames from XTC file: {xtc_file}")

    # Convert XTC to PDB
    pdb_files = step1_xtc_handling.xtc_to_pdb(
        xtc_file, 
        topology_file,
        stride, 
        outfolder,
        config
    )

    # Write PDB list to a file
    logger.info("Writing PDB list file")
    pdb_list_file = os.path.join(outfolder, 'pdb_list.ds')
    step1_xtc_handling.write_pdb_list(outfolder, pdb_list_file, config)

def detect_pockets(args, config):
    """Detect pockets on PDB files in input folder using p2rank."""
    infolder = os.path.abspath(args.infolder)
    outfolder = args.outfolder
    numthreads = args.numthreads
    logger = config['logger']
    pdb_list_file = os.path.join(infolder,'pdb_list.ds')

    # Run P2Rank on the PDB list
    logger.info(f"Running P2Rank with {numthreads} threads")
    step1_xtc_handling.run_p2rank(
        pdb_list_file, 
        outfolder,
        numthreads,
        config,
        )

    # Merge results to CSV
    logger.info("Merging pocket data to CSV")
    step2_pocketscsv.merge_to_csv(outfolder, pdb_list_file, config)

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="CLI tool to identify representative small-molecular binding pockets "
         " in protein molecular simulation trajectories using p2rank and via a simple clustering analysis"
    )
    parser.add_argument('--output', type=str, default='.', help='Base output directory (default: current directory)')
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Full pipeline
    parser_pocket = subparsers.add_parser("full_pipeline", help="Full pipeline for pocket prediction")
    parser_pocket.add_argument("--xtc", required=True, help="Input XTC trajectory file")
    parser_pocket.add_argument("--topology", required=True, help="Topology PDB/GRO file")
    parser_pocket.add_argument("--clustering_depth", type=int, required=True, help="Clustering depth")
    parser_pocket.add_argument("--numthreads", type=int, default=4, help="Number of parallel threads (default: 4)"),
    parser_pocket.add_argument("--stride", type=int, default=1, help="Stride applied for frame extraction"),
    parser_pocket.add_argument('--outfolder', type=str, default='.', help='Output folder')

    # Extract frames from XTC
    parser_convert = subparsers.add_parser("extract_to_pdb", help="Extract frames from XTC to PDB files")
    parser_convert.add_argument("--xtc", required=True, help="Input XTC trajectory file")
    parser_convert.add_argument("--topology", required=True, help="Topology PDB/GRO file")
    parser_convert.add_argument("--stride", type=int, default=1, help="Stride applied for frame extraction")
    parser_convert.add_argument('--outfolder', type=str, help='Output folder')

    # Predict pockets
    parser_convert = subparsers.add_parser("detect_pockets", help="Detect pockets on PDB files in input folder using p2rank.")
    parser_convert.add_argument("--numthreads", type=int, default=4, help="Number of parallel threads (default: 4)")
    parser_convert.add_argument('--infolder', type=str, help='Input folder containing PDB files (should be generated by extract_to_pdb)')
    parser_convert.add_argument('--outfolder', type=str, help='Output folder')

    args = parser.parse_args()

        if args.command is None:
        parser.print_help()
        return  # Exit without running anything
    
    # Load configuration with user-specified output directory
    config = get_config(args.outfolder)

    # Execute commands
    if args.command == "full_pipeline":
        full_pipeline(args, config)
    elif args.command == "extract_to_pdb":
        extract_to_pdb(args, config)
    elif args.command == "detect_pockets":
        detect_pockets(args, config)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()