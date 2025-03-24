import argparse
import os
import uuid
from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
import extract_predict  
import cluster

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

def get_config(output_dir):
    """Generate configuration based on the specified output directory."""

    os.makedirs(output_dir, exist_ok=True)

    # Timestamp for job ID
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") 
    job_id = f"{timestamp}_{str(uuid.uuid4())[:4]}"

    logger = setup_logging(os.path.join(output_dir,f"log_{job_id}.log"))
    
    config = {'logger': logger}

    return config

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
    pdb_files = extract_predict.xtc_to_pdb(xtc_file, topology_file,
        stride, outfolder, overwrite, config)

    # Write PDB list to a file
    logger.info("Writing PDB list file")
    pdb_list_file = os.path.join(outfolder, 'pdb_list.ds')
    extract_predict.write_pdb_list(outfolder, pdb_list_file, config)

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
    extract_predict.run_p2rank(pdb_list_file, outfolder, numthreads, novis, compress, overwrite, config)

def cluster_pockets(args, config):
    """Identifies representative pockets via clustering."""
    infile = os.path.abspath(args.infile)
    outfolder = args.outfolder
    method = args.method
    depth = args.depth
    min_prob = args.min_prob

    cluster.cluster_pockets(infile, outfolder, method, depth, min_prob, config)

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
    parser_full.add_argument("--method", type=str, choices=["hierarchical","dbscan"], 
                                help='Clustering method. Choose between hierarchical and dbscan.')
    parser_full.add_argument("--min_prob",type=float, default=0.5, help="Minimum ligand-binding probability for clustering (default: 0.5).")
    parser_full.add_argument('--novis', action='store_true', help='Do not keep visualizations generated by p2rank.')
    parser_full.add_argument('--compress', action='store_true', help='Compress p2rank output directory.')
    parser_full.add_argument("--numthreads", type=int, default=4, help="Number of parallel threads (default: 4)."),
    parser_full.add_argument("--stride", type=int, default=1, help="Stride applied for frame extraction."),
    parser_full.add_argument("--overwrite", action='store_true', help="Overwrite existing output directory.")
    parser_full.add_argument('--outfolder', type=str, default='.', help='Output folder.')

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
    parser_cluster.add_argument("--method", type=str, choices=["hierarchical","dbscan"], 
                                help='Clustering method. Choose between hierarchical and dbscan.')
    parser_cluster.add_argument("--depth", type=int, default=4, 
                                help='Clustering depth parameter for identification of representatives pockets from '
                                'hierarchical clustering.')
    parser_cluster.add_argument("--min_prob",type=float, default=0.5, help="Minimum ligand-binding probability for clustering (default: 0.5).")

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
    elif args.command == "cluster_pockets":
        cluster_pockets(args, config)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()